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
import numpy as np
import crop
import time


def main():
    print(os.getcwd())
    dataPath = os.path.join("data", "errano")
    settingsFolder = os.path.join(dataPath, "settings")
    fieldSettings = os.path.join(settingsFolder, "field.ini")

    if not importUtils.setField(fieldSettings):
        return

    # Soil
    print("Load soil...")
    soilFile = "soil.txt"
    soilPath = os.path.join(settingsFolder, soilFile)
    soil.readHorizon(soilPath)
    print("Soil depth [m]:", soil.horizon.lowerDepth)

    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(soil.horizon.lowerDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.maxThicknessAt)
    print("Nr. of layers:", C3DStructure.nrLayers)

    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    print("Nr. of cells: ", C3DStructure.nrCells)

    # todo modificare
    plantConfiguration = pd.read_csv(os.path.join(settingsFolder, "plant.csv"))
    crop.initializeCrop(plantConfiguration)

    criteria3D.initializeMesh()

    # initialize balance
    waterBalance.initializeBalance()
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    print("Read weather data...")
    weatherDataFolder = "meteo"
    weatherDataPath = os.path.join(dataPath, weatherDataFolder)
    weatherData = importUtils.readMeteoData(weatherDataPath)

    print("Read irrigation data...")
    waterFolder = "water"
    waterPath = os.path.join(dataPath, waterFolder)
    waterData = importUtils.readWaterData(waterPath, weatherData.iloc[0]["timestamp"],
                                          weatherData.iloc[-1]["timestamp"])
    weatherData.set_index(["timestamp"])
    waterData.set_index(["timestamp"])
    weatherData, waterData = importUtils.transformDates(weatherData, waterData)
    print("Total simulation time [hours]:", len(weatherData))

    # initialize export
    outputPath = os.path.join(dataPath, "output")
    exportUtils.createExportFile(outputPath)
    stateFolder = os.path.join(dataPath, "state")
    modelStateFileName = os.path.join(stateFolder, "modelState.bin")

    visual3D.initialize(1200)

    if C3DParameters.isAssimilation:
        print("Load obs water potential...")
        obsPath = os.path.join(dataPath, "obs_data")
        obsWaterPotential = pd.read_csv(os.path.join(obsPath, "waterPotential.csv"))
        obsStateFileName = os.path.join(stateFolder, "obsState.csv")

        for weatherIndex in range(1, 25):
            obsWeather = weatherData.loc[weatherIndex]
            waterEvent = waterData.loc[weatherIndex]
            currentDateTime = pd.to_datetime(obsWeather["timestamp"], unit='s')
            normTransmissivity = computeNormTransmissivity(weatherData, weatherIndex, C3DStructure.latitude,
                                                           C3DStructure.longitude)

            criteria3D.computeOneHour(obsWeather, waterEvent, normTransmissivity, currentDateTime)

            importUtils.writeObsState(obsStateFileName, obsWaterPotential, obsWeather["timestamp"])
            importUtils.loadObsState(obsStateFileName)
            exportUtils.takeScreenshot(obsWeather["timestamp"])

    # wait for start (press 'r')
    visual3D.isPause = True
    while visual3D.isPause:
        time.sleep(0.00001)

    # main cycle
    currentIndex = 1
    weatherIndex = 25
    restartIndex = 1
    isFirstRun = True
    while weatherIndex < len(weatherData):
        obsWeather = weatherData.loc[weatherIndex]
        waterEvent = waterData.loc[weatherIndex]
        currentDateTime = pd.to_datetime(obsWeather["timestamp"], unit='s')
        normTransmissivity = computeNormTransmissivity(weatherData, weatherIndex, C3DStructure.latitude, C3DStructure.longitude)

        criteria3D.computeOneHour(obsWeather, waterEvent, normTransmissivity, currentDateTime)

        # save model state
        if C3DParameters.isAssimilation and currentIndex == C3DParameters.assimilationInterval:
            importUtils.saveCurrentModelState(modelStateFileName)
            restartIndex = weatherIndex

        # save output
        if not C3DParameters.isForecast or isFirstRun:
            exportUtils.takeScreenshot(obsWeather["timestamp"])
        else:
            if currentIndex > (C3DParameters.forecastPeriod - C3DParameters.assimilationInterval):
                exportUtils.takeScreenshot(obsWeather["timestamp"])

        # restart
        if C3DParameters.isForecast and currentIndex == C3DParameters.forecastPeriod:
            importUtils.loadModelState(modelStateFileName)
            # assimilation
            if C3DParameters.isAssimilation:
                obsWeather = weatherData.loc[restartIndex]
                importUtils.writeObsState(obsStateFileName, obsWaterPotential, obsWeather["timestamp"])
                importUtils.loadObsState(obsStateFileName)
            # redraw
            waterBalance.totalTime = restartIndex * 3600
            visual3D.redraw()
            # re-initialize index
            weatherIndex = restartIndex + 1
            currentIndex = 1
            isFirstRun = False
        else:
            weatherIndex += 1
            currentIndex += 1

    visual3D.isPause = True
    print("\nEnd simulation.\n")

main()
