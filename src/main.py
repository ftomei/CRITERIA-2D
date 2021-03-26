import numpy as np
from dataStructures import *
from readDataFile import readDataFile
from fileUtilities import loadState
import soil
import waterBalance
import rectangularMesh
import criteria3D
#import visual3D
import os
from PenmanMonteith import computeHourlyET0
from transmissivity import computeNormTransmissivity
import exportUtils
import importUtils
import pandas as pd
import datetime
import crop

from hyperopt import fmin, tpe, hp, STATUS_OK, SparkTrials, Trials, pyll
from hyperopt.base import miscs_update_idxs_vals
from hyperopt.pyll.base import dfs, as_apply
from hyperopt.pyll.stochastic import implicit_stochastic_symbols
from sklearn.metrics import mean_squared_error, mean_absolute_error

SIMULATED_TOP_DEPTH = -0.174
SIMULATED_CENTER_DEPTH = -0.373
SIMULATED_BOTTOM_DEPTH = -0.566
REAL_TOP_DEPTH = -20
REAL_CENTER_DEPTH = -40
REAL_BOTTOM_DEPTH = -60

SIMULATED_ZERO_DISTANCE = 1.05
SIMULATED_LEFT_DISTANCE = 1.35
SIMULATED_CENTER_DISTANCE = 1.65
SIMULATED_RIGHT_DISTANCE = 1.95
REAL_ZERO_DISTANCE = 0
REAL_LEFT_DISTANCE = 30
REAL_CENTER_DISTANCE = 60
REAL_RIGHT_DISTANCE = 90

depths = {
    SIMULATED_TOP_DEPTH: REAL_TOP_DEPTH, 
    SIMULATED_CENTER_DEPTH: REAL_CENTER_DEPTH, 
    SIMULATED_BOTTOM_DEPTH: REAL_BOTTOM_DEPTH
}
distances = {
    SIMULATED_ZERO_DISTANCE: REAL_ZERO_DISTANCE, 
    SIMULATED_LEFT_DISTANCE: REAL_LEFT_DISTANCE, 
    SIMULATED_CENTER_DISTANCE: REAL_CENTER_DISTANCE, 
    SIMULATED_RIGHT_DISTANCE: REAL_RIGHT_DISTANCE
}


def objective(params):
    C3DParameters.waterTableDepth = -params["water_table"]
    #print (os.getcwd())
    dataPath = os.path.join("..", "data", "fondo_1_tuning_4")

    #print("Building rectangle mesh...")
    rectangularMesh.rectangularMeshCreation()
    #print ("Nr. of rectangles:", C3DStructure.nrRectangles)
    #print ("Total area [m^2]:", C3DStructure.totalArea)

    rectangularMesh.header = rectangularMesh.getHeader(rectangularMesh.C3DRM)

    # SOIL
    #print ("Load soil...")
    soilFolder = "soil"
    soilFile = "soil_fitting.txt"
    soilPath = os.path.join(dataPath, soilFolder, soilFile)
    #soil.C3DSoil = soil.readHorizon(soilPath, 1, 1.4e-6, 0.43, 1.6, 1.25)
    soil.C3DSoil = soil.readHorizon(soilPath, 1, params["k_sat"], params["theta_sat"], params["alpha"], params["n"])
    totalDepth = soil.C3DSoil.lowerDepth
    #print("Soil depth [m]:", totalDepth)

    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(totalDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.geometricFactor)
    #print("Nr. of layers:", C3DStructure.nrLayers)

    # Initialize memory
    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    #print("Nr. of cells: ", C3DStructure.nrCells)

    #print("Set cell properties...")
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

    #print("Set links...")
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

    waterBalance.initializeBalance()
    #print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    #print("Read weather data...")
    weatherDataFolder = "arpae"
    weatherDataPath = os.path.join(dataPath, weatherDataFolder)
    stationInfo, weatherData = importUtils.readArpaeData(weatherDataPath)
    height = stationInfo.iloc[0]["Altezza (Metri sul livello del mare)"]

    #print("Read irrigations data...")
    waterFolder = "water"
    waterPath = os.path.join(dataPath, waterFolder)
    irrigationConfigurations, waterData = importUtils.readWaterData(waterPath, weatherData.iloc[0]["start"],
                                                                     weatherData.iloc[-1]["start"])
    criteria3D.setDripIrrigationPositions(irrigationConfigurations)

    # weatherData, waterData = importUtils.transformDates(weatherData, waterData)
    soilPath = os.path.join(dataPath, "soil")
    plantConfiguration = pd.read_csv(os.path.join(soilPath, "plant_configuration.csv"))
    crop.initializeCrop(plantConfiguration, irrigationConfigurations, params['root_max_distance'], params['LAI'], params['roots_depth'], params['roots_deformation'], params['kc_max'])

    # TIME LENGTH
    weatherTimeLength = (weatherData.iloc[0]["end"] - weatherData.iloc[0]["start"])  # [s]
    waterTimeLength = (waterData.iloc[0]["end"] - waterData.iloc[0]["start"])  # [s]
    #print("Weather relevations time lenght [s]:", weatherTimeLength)
    #print("Water relevations time lenght [s]:", waterTimeLength)
    if (weatherTimeLength % waterTimeLength) != 0:
        raise Exception("Water time lenght is not a divider of weather data time lenght")
    else:
        nrWaterEventsInWeatherTimeLength = int(weatherTimeLength / waterTimeLength)
    #print("Total simulation time [hours]:", len(weatherData) * weatherTimeLength / 3600)

    #visual3D.initialize(1280)
    #visual3D.isPause = True

    # initialize export
    outputPath = os.path.join(dataPath, "output")
    exportUtils.createExportFile(outputPath)

    latitude = stationInfo.iloc[0]["Latitudine (Gradi Centesimali)"]
    longitude = stationInfo.iloc[0]["Longitudine (Gradi Centesimali)"]

    # main cycle
    extendedWeatherData, extendedWaterData = importUtils.setDataIndeces(weatherData, waterData)
    weatherData, waterData = extendedWeatherData.iloc[12:-12], extendedWaterData.iloc[12:-12]
    minTimestamp, maxTimestamp = weatherData["end"].min(), weatherData["end"].max()
    dailyET0 = 0
    simulated_data = pd.DataFrame()
    
    waterTableDate, waterTableDepth = importUtils.readWaterTable(waterPath)

    for weatherIndex, obsWeather in weatherData.iterrows():
        currentDateTime = pd.to_datetime(obsWeather["end"], unit='s')
        # for i in range(len(waterTableDepth)):
        #    if currentDateTime > waterTableDate[i]:
        #        C3DParameters.waterTableDepth = waterTableDepth[i]

        airTemperature = obsWeather["temperature"]
        globalSWRadiation = obsWeather["radiations"]
        airRelHumidity = obsWeather["humidity"]
        windSpeed_10m = obsWeather["wind"]

        # evapotranspiration
        normTransmissivity = computeNormTransmissivity(extendedWeatherData, currentDateTime, latitude, longitude)
        ET0 = computeHourlyET0(height, airTemperature, globalSWRadiation, airRelHumidity, windSpeed_10m,
                               normTransmissivity)  # mm m^-2
        #print(currentDateTime, "ET0:", format(ET0, ".2f"))

        crop.setEvapotranspiration(ET0)

        for i in range(nrWaterEventsInWeatherTimeLength):

            waterIndex = weatherIndex + (i * waterTimeLength)
            waterEvent = waterData.loc[waterIndex]
            try:
                current_timestamp = waterEvent["end"]
            except:
                current_timestamp = waterIndex + waterTimeLength
            
            #print('{:.2f}'.format(((current_timestamp - minTimestamp)/(maxTimestamp - minTimestamp))*100))

            waterBalance.currentPrec = waterEvent["precipitations"] / waterTimeLength * 3600  # [mm m-2 hour-1]
            criteria3D.setRainfall(waterEvent["precipitations"], waterTimeLength)

            if C3DParameters.assignIrrigation:
                waterBalance.currentIrr = (len(criteria3D.irrigationIndeces) * waterEvent[
                    "irrigations"]) / waterTimeLength * 3600  # [l hour-1]
                criteria3D.setDripIrrigation(waterEvent["irrigations"], waterTimeLength)

            if (waterBalance.currentIrr > 0) or (waterBalance.currentPrec > 0):
                C3DParameters.deltaT_max = 300
                C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT, C3DParameters.deltaT_max)
            else:
                C3DParameters.deltaT_max = waterTimeLength

            criteria3D.compute(waterTimeLength)

            exportUtils.takeScreenshot(waterEvent["end"])
            for elem in exportUtils.getScreenshot(waterEvent["end"]):
                simulated_data = simulated_data.append(elem, ignore_index=True)

    print("\nEnd simulation.")


    #print ("\nEnd simulation.")
    original_data = pd.read_csv(os.path.join(dataPath, 'ground_truth.csv'))
    original_data = original_data.pivot(index='timestamp', columns=['z', 'x'], values='H')
    original_data.columns = ["_".join((str(z), str(x))) for z, x in original_data.columns]
    original_data = original_data.reset_index()

    #simulated_data = pd.read_csv(os.path.join(os.path.join(dataPath, 'output'), 'output.csv'))
    #simulated_data = simulated_data.drop(['y', 'Se'], axis=1)
    simulated_data = simulated_data[(simulated_data["z"] == SIMULATED_TOP_DEPTH) | (simulated_data["z"] == SIMULATED_CENTER_DEPTH) | (simulated_data["z"] == SIMULATED_BOTTOM_DEPTH)]
    simulated_data = simulated_data[(simulated_data["x"] == SIMULATED_ZERO_DISTANCE) | (simulated_data["x"] == SIMULATED_LEFT_DISTANCE) | (simulated_data["x"] == SIMULATED_CENTER_DISTANCE) | (simulated_data["x"] == SIMULATED_RIGHT_DISTANCE)]
    simulated_data = simulated_data.pivot(index='timestamp', columns=['z', 'x'], values='H')
    simulated_data.columns = ["_".join((str(z), str(x))) for z, x in simulated_data.columns]
    simulated_data = simulated_data.reset_index()
    for z in depths.keys():
        for x in distances.keys():
            simulated_data = simulated_data.rename(columns={"_".join([str(z), str(x)]): "_".join([str(depths[z]), str(distances[x])])})

    original_data_columns = original_data.columns.tolist()
    original_data_columns.sort()
    original_data = original_data[original_data_columns]

    simulated_data_columns = simulated_data.columns.tolist()
    simulated_data_columns.sort()
    simulated_data = simulated_data[simulated_data_columns]

    #print(original_data)
    #print(simulated_data)

    original_data = original_data.iloc[:, :-1].to_numpy() * -1
    simulated_data = simulated_data.iloc[:, :-1].to_numpy() * -1
    
    original_data = np.nan_to_num(np.log(original_data))
    simulated_data = np.nan_to_num(np.log(simulated_data))

    total_rmse = mean_squared_error(simulated_data, original_data, squared=False)
    
    return {'loss': total_rmse, 'status': STATUS_OK}

def main():
    dataPath = os.path.join("..", "data", "fondo_1_tuning_4")
    space = {
        'k_sat': hp.loguniform('k_sat', -14.5, -10),
        'theta_sat': hp.uniform('theta_sat', 0.3, 0.5),
        'alpha': hp.uniform('alpha', 1, 3),
        'n': hp.uniform('n', 1.01, 1.5),
        'water_table': hp.uniform('water_table', 1.5, 3.5),
        'LAI': hp.uniform('LAI', 2, 4),
        'roots_depth': hp.uniform('roots_depth', 0.5, 1.2),
        'roots_deformation': hp.uniform('roots_deformation', 0, 2),
        'kc_max': hp.uniform('kc_max', 1.5, 2.5),
        'root_max_distance': hp.uniform('root_max_distance', 1.5, 2.5)
        }
    trials = SparkTrials()
    best = fmin(fn=objective,
                space=space,
                algo=tpe.suggest,
                max_evals=150,
                trials=trials,
                show_progressbar=True)

    best_parameters = pd.DataFrame.from_records([best])
    best_parameters.to_csv(os.path.join(dataPath, 'best_parameters.csv'), index=False)

    all_trials = pd.DataFrame(trials.trials)
    all_trials.to_json(os.path.join(dataPath, 'all_trials.json'), indent=True)

    best_trial = pd.DataFrame(trials.best_trial)
    best_trial.to_json(os.path.join(dataPath,'best_trial.json'), indent=True)

    objective(best)

main()
