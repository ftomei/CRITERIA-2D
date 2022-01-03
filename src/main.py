from dataStructures import *
import soil
import waterBalance
import rectangularMesh
import criteria3D
import visual3D
import os
from PenmanMonteith import computeHourlyET0
from transmissivity import computeNormTransmissivity
import exportUtils
import importUtils
import pandas as pd
import crop


def main():
    print(os.getcwd())
    dataPath = os.path.join("data", "errano")
    settingsFolder = os.path.join(dataPath, "settings")

    print("Building rectangle mesh...")
    rectangularMesh.rectangularMeshCreation()
    print("Nr. of rectangles:", C3DStructure.nrRectangles)
    print("Total area [m^2]:", C3DStructure.totalArea)

    rectangularMesh.header = rectangularMesh.getHeader(rectangularMesh.C3DRM)

    # SOIL
    print("Load soil...")
    soilFile = "soil.txt"
    soilPath = os.path.join(settingsFolder, soilFile)
    soil.C3DSoil = soil.readHorizon(soilPath, 1)
    totalDepth = soil.C3DSoil.lowerDepth
    print("Soil depth [m]:", totalDepth)

    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(totalDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.geometricFactor)
    print("Nr. of layers:", C3DStructure.nrLayers)

    # Initialize memory
    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    print("Nr. of cells: ", C3DStructure.nrCells)

    print("Set cell properties...")
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

    print("Set links...")
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
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    print("Read drip position...")
    irrigationConfigurations = pd.read_csv(os.path.join(settingsFolder, "drippers.csv"))
    criteria3D.setDripIrrigationPositions(rectangularMesh, irrigationConfigurations)

    print("Read plant position...")
    plantConfiguration = pd.read_csv(os.path.join(settingsFolder, "plant.csv"))
    crop.initializeCrop(plantConfiguration, irrigationConfigurations)

    print("Read weather data...")
    weatherDataFolder = "arpae"
    weatherDataPath = os.path.join(dataPath, weatherDataFolder)
    stationInfo, weatherData = importUtils.readArpaeData(weatherDataPath)
    height = stationInfo.iloc[0]["Altezza (Metri sul livello del mare)"]

    print("Read irrigation data...")
    waterFolder = "water"
    waterPath = os.path.join(dataPath, waterFolder)
    waterData = importUtils.readWaterData(waterPath, weatherData.iloc[0]["start"], weatherData.iloc[-1]["start"])

    # weatherData, waterData = importUtils.transformDates(weatherData, waterData)

    # TIME LENGTH
    weatherTimeLength = (weatherData.iloc[0]["end"] - weatherData.iloc[0]["start"])  # [s]
    waterTimeLength = (waterData.iloc[0]["end"] - waterData.iloc[0]["start"])  # [s]
    print("Weather relevations time lenght [s]:", weatherTimeLength)
    print("Water relevations time lenght [s]:", waterTimeLength)
    if (weatherTimeLength % waterTimeLength) != 0:
        raise Exception("Water time lenght is not a divider of weather data time lenght")
    else:
        nrWaterEventsInWeatherTimeLength = int(weatherTimeLength / waterTimeLength)
    print("Total simulation time [hours]:", len(weatherData) * weatherTimeLength / 3600)

    visual3D.initialize(1280)
    visual3D.isPause = True

    # initialize export
    outputPath = os.path.join(dataPath, "output")
    exportUtils.createExportFile(outputPath)

    latitude = stationInfo.iloc[0]["Latitudine (Gradi Centesimali)"]
    longitude = stationInfo.iloc[0]["Longitudine (Gradi Centesimali)"]

    # main cycle
    extendedWeatherData, extendedWaterData = importUtils.setDataIndeces(weatherData, waterData)
    weatherData, waterData = extendedWeatherData.iloc[12:-12], extendedWaterData.iloc[12:-12]
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
        print(currentDateTime, "ET0:", format(ET0, ".2f"))

        crop.setEvapotranspiration(ET0)

        for i in range(nrWaterEventsInWeatherTimeLength):

            waterIndex = weatherIndex + (i * waterTimeLength)
            waterEvent = waterData.loc[waterIndex]

            waterBalance.currentPrec = waterEvent["precipitations"] / waterTimeLength * 3600  # [mm m-2 hour-1]
            criteria3D.setRainfall(waterEvent["precipitations"], waterTimeLength)

            if C3DParameters.assignIrrigation:
                waterBalance.currentIrr = (len(criteria3D.irrigationIndices) * waterEvent[
                    "irrigations"]) / waterTimeLength * 3600  # [l hour-1]
                criteria3D.setDripIrrigation(waterEvent["irrigations"], waterTimeLength)

            if (waterBalance.currentIrr > 0) or (waterBalance.currentPrec > 0):
                C3DParameters.deltaT_max = 300
                C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT, C3DParameters.deltaT_max)
            else:
                C3DParameters.deltaT_max = waterTimeLength

            criteria3D.compute(waterTimeLength)

            exportUtils.takeScreenshot(waterEvent["end"])

    print("\nEnd simulation.")


main()
