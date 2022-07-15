import pandas as pd
import os
import time

from dataStructures import *
import soil
import waterBalance
import criteria3D
import visual3D
import exportUtils
import importUtils
import crop


def main():
    # path
    print(os.getcwd())
    projectPath = os.path.join("data", "errano")
    settingsFolder = os.path.join(projectPath, "settings")
    weatherFolder = os.path.join(projectPath, "meteo")
    waterFolder = os.path.join(projectPath, "water")
    obsDataFolder = os.path.join(projectPath, "obs_data")
    stateFolder = os.path.join(projectPath, "state")
    outputFolder = os.path.join(projectPath, "output")

    print("Read field settings...")
    fieldSettings = os.path.join(settingsFolder, "field.ini")
    if not importUtils.readFieldParameters(fieldSettings):
        return

    # todo modificare
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

    print("Read crop settings...")
    cropSettings = os.path.join(settingsFolder, "crop.ini")
    if not importUtils.readCropParameters(cropSettings):
        return
    crop.initializeCrop()

    print("Initialize mesh...")
    criteria3D.initializeMesh()

    waterBalance.initializeBalance()
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    print("Read weather and irrigation data...")
    weatherData = importUtils.readMeteoData(weatherFolder)
    waterData = importUtils.readWaterData(waterFolder, weatherData.iloc[0]["timestamp"],
                                          weatherData.iloc[-1]["timestamp"])
    weatherData.set_index(["timestamp"])
    waterData.set_index(["timestamp"])
    criteria3D.weatherData, criteria3D.waterData = importUtils.transformDates(weatherData, waterData)
    print("Total simulation time [hours]:", len(weatherData))

    # initialize export
    exportUtils.createExportFile(outputFolder)
    modelStateFileName = os.path.join(stateFolder, "modelState.bin")

    if C3DParameters.isPeriodicAssimilation or C3DParameters.isFirstAssimilation:
        print("Read observed water potential...")
        obsWaterPotential = pd.read_csv(os.path.join(obsDataFolder, "waterPotential.csv"))
        obsFileName = os.path.join(stateFolder, "obsState.csv")

    weatherIndex = 1
    if C3DParameters.isFirstAssimilation:
        print("Assimilate observed water potential...")
        obsWeather = weatherData.loc[weatherIndex]
        importUtils.writeObsData(obsFileName, obsWaterPotential, obsWeather["timestamp"])
        importUtils.loadObsData(obsFileName)

    visual3D.initialize(1200)
    visual3D.isPause = True
    # wait for start (press 'r')
    while visual3D.isPause:
        time.sleep(0.00001)

    # main cycle
    currentIndex = 1
    restartIndex = 1
    isFirstRun = True
    isRedraw = True
    while weatherIndex < len(weatherData):
        criteria3D.computeOneHour(weatherIndex, isRedraw)
        obsWeather = criteria3D.weatherData.loc[weatherIndex]

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
            obsWeather = criteria3D.weatherData.loc[restartIndex]
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
