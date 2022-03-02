from dataStructures import *
import soil
import waterBalance
import rectangularMesh
import criteria3D
#import visual3D
import assimilation
import os
from PenmanMonteith import computeHourlyET0
from transmissivity import computeNormTransmissivity
import exportUtils
import importUtils
import pandas as pd
import numpy as np
import crop
import time

from hyperopt import fmin, tpe, hp, STATUS_OK, SparkTrials, Trials, pyll
from hyperopt.base import miscs_update_idxs_vals
from hyperopt.pyll.base import dfs, as_apply
from hyperopt.pyll.stochastic import implicit_stochastic_symbols
from sklearn.metrics import mean_squared_error, mean_absolute_error

iteration = 0

def objective(params):
    global iteration
    # print(os.getcwd())
    dataPath = os.path.join("data", "errano")
    settingsFolder = os.path.join(dataPath, "settings")

    # print("Building rectangle mesh...")
    rectangularMesh.rectangularMeshCreation()
    # print("Nr. of rectangles:", C3DStructure.nrRectangles)
    # print("Total area [m^2]:", C3DStructure.totalArea)

    rectangularMesh.header = rectangularMesh.getHeader(rectangularMesh.C3DRM)

    # SOIL
    # print("Load soil...")
    soilFile = "soil.txt"
    soilPath = os.path.join(settingsFolder, soilFile)
    soil.readHorizon(soilPath, 1)
    soil.horizon.VG_alpha = params["VG_alpha"]
    soil.horizon.VG_n = params["VG_n"]
    soil.horizon.thetaS = params["thetaS"]
    soil.horizon.Ks = params["Ks"]
    totalDepth = soil.horizon.lowerDepth
    # print("Soil depth [m]:", totalDepth)

    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(totalDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.maxThicknessAt)
    # print("Nr. of layers:", C3DStructure.nrLayers)

    # Initialize memory
    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    # print("Nr. of cells: ", C3DStructure.nrCells)

    # print("Set cell properties...")
    for i in range(C3DStructure.nrRectangles):
        [x, y, z] = rectangularMesh.C3DRM[i].centroid
        for layer in range(C3DStructure.nrLayers):
            index = i + C3DStructure.nrRectangles * layer
            elevation = z - soil.depth[layer]
            volume = float(rectangularMesh.C3DRM[i].area * soil.thickness[layer])
            criteria3D.setCellGeometry(index, x, y,
                                       elevation, volume, rectangularMesh.C3DRM[i].area)
            if layer == 0:
                # surface
                if rectangularMesh.C3DRM[i].isBoundary and C3DParameters.isSurfaceRunoff:
                    criteria3D.setCellProperties(index, True, BOUNDARY_RUNOFF)
                    criteria3D.setBoundaryProperties(index,
                                                     rectangularMesh.C3DRM[i].boundarySide,
                                                     rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    criteria3D.setCellProperties(index, True, BOUNDARY_NONE)

                criteria3D.setMatricPotential(index, 0.0)

            elif layer == (C3DStructure.nrLayers - 1):
                # last layer
                if C3DParameters.isWaterTable:
                    criteria3D.setCellProperties(index, False, BOUNDARY_PRESCRIBEDTOTALPOTENTIAL)
                elif C3DParameters.isFreeDrainage:
                    criteria3D.setCellProperties(index, False, BOUNDARY_FREEDRAINAGE)
                else:
                    criteria3D.setCellProperties(index, False, BOUNDARY_NONE)

                criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)

            else:
                if rectangularMesh.C3DRM[i].isBoundary and C3DParameters.isFreeLateralDrainage:
                    criteria3D.setCellProperties(index, False, BOUNDARY_FREELATERALDRAINAGE)
                    criteria3D.setBoundaryProperties(index,
                                                     rectangularMesh.C3DRM[i].boundarySide * soil.thickness[layer],
                                                     rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    criteria3D.setCellProperties(index, False, BOUNDARY_NONE)

                criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)

    # print("Set links...")
    for i in range(C3DStructure.nrRectangles):
        # UP
        for layer in range(1, C3DStructure.nrLayers):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i
            linkIndex = index - C3DStructure.nrRectangles
            criteria3D.SetCellLink(index, linkIndex, UP, exchangeArea)
        # LATERAL
        for neighbour in rectangularMesh.C3DRM[i].neighbours:
            if neighbour != NOLINK:
                linkSide = rectangularMesh.getAdjacentSide(i, neighbour)
                for layer in range(C3DStructure.nrLayers):
                    if layer == 0:
                        # surface: boundary length [m]
                        exchangeArea = linkSide
                    else:
                        # sub-surface: boundary area [m2]
                        exchangeArea = soil.thickness[layer] * linkSide
                    index = C3DStructure.nrRectangles * layer + i
                    linkIndex = C3DStructure.nrRectangles * layer + neighbour
                    criteria3D.SetCellLink(index, linkIndex, LATERAL, exchangeArea)
        # DOWN
        for layer in range(C3DStructure.nrLayers - 1):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i
            linkIndex = index + C3DStructure.nrRectangles
            criteria3D.SetCellLink(index, linkIndex, DOWN, exchangeArea)

    # initial state
    stateFolder = os.path.join(dataPath, "state")
    initialState = pd.read_csv(os.path.join(stateFolder, "1630447200.csv"))
    assimilation.assimilate(initialState)
    waterBalance.initializeBalance()
    # print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    # print("Read drip position...")
    irrigationConfigurations = pd.read_csv(os.path.join(settingsFolder, "dripper.csv"))
    criteria3D.setDripIrrigationPositions(irrigationConfigurations)

    # print("Read plant position...")
    plantConfiguration = pd.read_csv(os.path.join(settingsFolder, "plant.csv"))
    criteria3D.setPlantPositions(plantConfiguration)
    crop.initializeCrop(plantConfiguration, params['rootDepthMax'], params['rootXDeformation'], params['rootZDeformation'], params['kcMax'])

    # print("Read weather data...")
    weatherDataFolder = "meteo"
    weatherDataPath = os.path.join(dataPath, weatherDataFolder)
    stationInfo, weatherData = importUtils.readMeteoData(weatherDataPath)
    height = stationInfo.iloc[0]["Height"]

    # print("Read irrigation data...")
    waterFolder = "water"
    waterPath = os.path.join(dataPath, waterFolder)
    waterData = importUtils.readWaterData(waterPath, weatherData.iloc[0]["timestamp"],
                                          weatherData.iloc[-1]["timestamp"])
    weatherData.set_index(["timestamp"])
    waterData.set_index(["timestamp"])
    weatherData, waterData = importUtils.transformDates(weatherData, waterData)

    # TIME LENGTH
    weatherTimeLength = (weatherData.iloc[1]["timestamp"] - weatherData.iloc[0]["timestamp"])  # [s]
    waterTimeLength = (waterData.iloc[1]["timestamp"] - waterData.iloc[0]["timestamp"])  # [s]
    # print("Weather data time length [s]:", weatherTimeLength)
    # print("Water data time length [s]:", waterTimeLength)
    if (weatherTimeLength % waterTimeLength) != 0:
        raise Exception("Water time length is not a divider of weather data time length")
    else:
        nrWaterEventsInTimeLength = int(weatherTimeLength / waterTimeLength)
    # print("Total simulation time [hours]:", len(weatherData) * weatherTimeLength / 3600)

    # initialize export
    outputPath = os.path.join(dataPath, "output")
    outputFileName = f"output_{str(iteration)}.csv"
    outputFilePath = os.path.join(outputPath, outputFileName)
    outputPointsPath = os.path.join(outputPath, "output_points.csv")
    exportUtils.createExportFile(outputPointsPath, outputFilePath)

    latitude = stationInfo.iloc[0]["Latitude"]
    longitude = stationInfo.iloc[0]["Longitude"]

    #visual3D.initialize(1280)
    #visual3D.isPause = True
    # wait for start
    #while visual3D.isPause:
        #time.sleep(0.00001)

    # main cycle
    weatherIndex = 0
    while weatherIndex < len(weatherData):
        obsWeather = weatherData.loc[weatherIndex]
        currentDateTime = pd.to_datetime(obsWeather["timestamp"], unit='s')

        # waterTable
        # for i in range(len(waterTableDepth)):
        #    if currentDateTime > waterTableDate[i]:
        #        C3DParameters.waterTableDepth = waterTableDepth[i]

        if not (np.isnan(obsWeather["air_temperature"])):
            airTemperature = obsWeather["air_temperature"]
        if not (np.isnan(obsWeather["solar_radiation"])):
            globalSWRadiation = obsWeather["solar_radiation"]
        if not (np.isnan(obsWeather["air_humidity"])):
            airRelHumidity = obsWeather["air_humidity"]
        if not (np.isnan(obsWeather["wind_speed"])):
            windSpeed_10m = obsWeather["wind_speed"]
        else:
            print("Missed")

        # evapotranspiration
        normTransmissivity = computeNormTransmissivity(weatherData, weatherIndex, latitude, longitude)
        ET0 = computeHourlyET0(height, airTemperature, globalSWRadiation, airRelHumidity, windSpeed_10m,
                               normTransmissivity)  # mm m^-2
        # print(currentDateTime, "ET0:", format(ET0, ".2f"))
        criteria3D.initializeSinkSource(ALL)
        crop.setEvapotranspiration(currentDateTime, ET0)

        for i in range(nrWaterEventsInTimeLength):
            waterIndex = weatherIndex * nrWaterEventsInTimeLength + i
            waterEvent = waterData.loc[waterIndex]
            precipitation = 0 if np.isnan(waterEvent["precipitation"]) else waterEvent["precipitation"]
            irrigation = 0 if np.isnan(waterEvent["irrigation"]) else waterEvent["irrigation"]

            criteria3D.initializeSinkSource(ONLY_SURFACE)
            waterBalance.currentPrec = precipitation / waterTimeLength * 3600.  # [mm m-2 hour-1]
            criteria3D.setRainfall(precipitation, waterTimeLength)

            if C3DParameters.assignIrrigation:
                waterBalance.currentIrr = irrigation / waterTimeLength * 3600.  # [l hour-1]
                criteria3D.setDripIrrigation(irrigation, waterTimeLength)

            if (waterBalance.currentIrr > 0) or (waterBalance.currentPrec > 0):
                C3DParameters.deltaT_max = 300
                C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT, C3DParameters.deltaT_max)
            else:
                C3DParameters.deltaT_max = waterTimeLength

            criteria3D.compute(waterTimeLength, False)

        exportUtils.takeScreenshot(obsWeather["timestamp"], outputFilePath)
        weatherIndex += 1

    #visual3D.isPause = True
    # print("\nEnd simulation.\n")

    simulated_data = pd.read_csv(outputFilePath)
    simulated_data = simulated_data.set_index('timestamp')
    simulated_data *= -1
    simulated_data[simulated_data < 20] = 20
    simulated_data = simulated_data.apply(lambda x: np.log(x))

    
    original_data = pd.read_csv(os.path.join(dataPath, 'ground_truth.csv'))
    original_data = original_data.set_index('timestamp')
    original_data = original_data.loc[simulated_data.index]
    original_data *= -1
    original_data[original_data < 20] = 20
    original_data = original_data.apply(lambda x: np.log(x))

    total_rmse = mean_squared_error(simulated_data, original_data, squared=False)
    
    iteration += 1
    return {'loss': total_rmse, 'status': STATUS_OK}

def main():
    print('Start')
    start_time = time.time()
    dataPath = os.path.join("data", "errano")
    outputPath = os.path.join(dataPath, "output")
    space = {
        'VG_alpha': hp.uniform('VG_alpha', 1.0, 3.0),
        'VG_n': hp.uniform('VG_n', 1.1, 1.4),
        'thetaS': hp.uniform('thetaS', 0.3, 0.5),
        'Ks': hp.loguniform('Ks', -16.1, -12.4),     
        #'water_table': hp.uniform('water_table', 1.5, 3.5),
        #'LAI': hp.uniform('LAI', 2, 4),
        'rootDepthMax': hp.uniform('rootDepthMax', 0.6, 1.0),
        'rootXDeformation': hp.uniform('rootXDeformation', 0.0, 1.0),
        'rootZDeformation': hp.uniform('rootZDeformation', 0.0, 1.0),
        'kcMax': hp.uniform('kcMax', 1.5, 2.5),
        #'root_max_distance': hp.uniform('root_max_distance', 1.5, 2.5)
        }
    trials = Trials()
    best = fmin(fn=objective,
                space=space,
                algo=tpe.suggest,
                max_evals=150,
                trials=trials,
                show_progressbar=True,
                rstate=np.random.RandomState(42))

    best_parameters = pd.DataFrame.from_records([best])
    best_parameters.to_csv(os.path.join(outputPath, 'best_parameters.csv'), index=False)

    all_trials = pd.DataFrame(trials.trials)
    all_trials.to_json(os.path.join(outputPath, 'all_trials.json'), indent=True)

    best_trial = pd.DataFrame(trials.best_trial)
    best_trial.to_json(os.path.join(outputPath,'best_trial.json'), indent=True)

    objective(best)
    end_time = time.time()
    print("Time of execution:", end_time-start_time)

main()
