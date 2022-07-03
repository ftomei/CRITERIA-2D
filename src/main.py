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
    # path
    print(os.getcwd())
    projectPath = os.path.join("data", "errano")
    settingsFolder = os.path.join(projectPath, "settings")
    waterFolder = os.path.join(projectPath, "water")
    weatherFolder = os.path.join(projectPath, "meteo")
    obsDataFolder = os.path.join(projectPath, "obs_data")
    stateFolder = os.path.join(projectPath, "state")
    outputFolder = os.path.join(projectPath, "output")

    fieldSettings = os.path.join(settingsFolder, "field.ini")
    if not importUtils.setField(fieldSettings):
        return

    # Soil
    print("Load soil...")
    soilFile = "soil.txt"
    soilPath = os.path.join(settingsFolder, soilFile)
    soil.readHorizon(soilPath)

    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(C3DStructure.gridDepth,
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
    waterBalance.initializeBalance()
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    print("Read weather data...")
    weatherData = importUtils.readMeteoData(weatherFolder)

    print("Read precipitation and irrigation data...")
    waterData = importUtils.readWaterData(waterFolder, weatherData.iloc[0]["timestamp"],
                                          weatherData.iloc[-1]["timestamp"])
    weatherData.set_index(["timestamp"])
    waterData.set_index(["timestamp"])
    weatherData, waterData = importUtils.transformDates(weatherData, waterData)
    print("Total simulation time [hours]:", len(weatherData))

    # initialize export
    exportUtils.createExportFile(outputFolder)
    modelStateFileName = os.path.join(stateFolder, "modelState.bin")

    weatherIndex = 1
    if C3DParameters.isFirstAssimilation:
        print("Assimilate observed water potential...")
        obsWaterPotential = pd.read_csv(os.path.join(obsDataFolder, "waterPotential.csv"))
        obsFileName = os.path.join(stateFolder, "obsState.csv")
        obsWeather = weatherData.loc[weatherIndex]
        importUtils.writeObsData(obsFileName, obsWaterPotential, obsWeather["timestamp"])
        importUtils.loadObsData(obsFileName)

    visual3D.initialize(1200)
    # wait for start (press 'r')
    visual3D.isPause = True
    while visual3D.isPause:
        time.sleep(0.00001)

    # main cycle
    currentIndex = 1
    restartIndex = 1
    isFirstRun = True
    while weatherIndex < len(weatherData):
        obsWeather = weatherData.loc[weatherIndex]
        waterEvent = waterData.loc[weatherIndex]
        currentDateTime = pd.to_datetime(obsWeather["timestamp"], unit='s')
        normTransmissivity = computeNormTransmissivity(weatherData, weatherIndex, C3DStructure.latitude, C3DStructure.longitude)

        criteria3D.computeOneHour(obsWeather, waterEvent, normTransmissivity, currentDateTime)

        # assimilation
        if C3DParameters.isPeriodicAssimilation and not C3DParameters.isForecast:
            if (currentIndex % C3DParameters.assimilationInterval) == 0:
                importUtils.writeObsData(obsFileName, obsWaterPotential, obsWeather["timestamp"])
                importUtils.loadObsData(obsFileName)

        # save model state
        if C3DParameters.isForecast and currentIndex == C3DParameters.assimilationInterval:
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
            obsWeather = weatherData.loc[restartIndex]
            importUtils.writeObsData(obsFileName, obsWaterPotential, obsWeather["timestamp"])
            importUtils.loadObsData(obsFileName)
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
