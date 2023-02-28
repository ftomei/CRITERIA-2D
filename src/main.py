import pandas as pd
import os
import time
import datetime

from dataStructures import *
import soil
import waterBalance
import criteria3D
import visual3D
import exportUtils
import importUtils
import crop
import json
import argparse


def main(args):
    # path
    print(f"pwd: {os.getcwd()}")
    print(f"args: {args}")
    settingsFolder = os.path.join(args.path, "settings")
    weatherFolder = os.path.join(args.path, "meteo")
    waterFolder = os.path.join(args.path, "water")
    obsDataFolder = os.path.join(args.path, "obs_data")

    stateFolder = os.path.join(args.path, "state")
    if not os.path.exists(stateFolder):
        os.makedirs(stateFolder)

    outputFolder = os.path.join(args.path, "output")
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    try:
        params = json.load(open(os.path.join(outputFolder, f"input_{args.iteration}.json")))
    except:
        params = {}
    iterations_str = f"_{args.iteration}" if args.iteration != -1 else ""
    params["iteration"] = args.iteration
    print(f"tuning params: {params}")

    print("Read model settings...")
    modelSettings = os.path.join(settingsFolder, "settings.ini")
    if not importUtils.readModelParameters(modelSettings, params):
        return

    print("Read field settings...")
    fieldSettings = os.path.join(settingsFolder, "field.ini")
    if not importUtils.readFieldParameters(fieldSettings, params):
        print("Error in read field.ini")
        return

    print("read soil properties...")
    soilSettings = os.path.join(settingsFolder, "soil.csv")
    soil.readHorizon(soilSettings, params)
    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(C3DStructure.gridDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.maxThicknessAt)
    print("Nr. of layers:", C3DStructure.nrLayers)

    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    print("Nr. of cells: ", C3DStructure.nrCells)

    print("Initialize mesh...")
    criteria3D.initializeMesh()

    print("Read crop settings...")
    cropSettings = os.path.join(settingsFolder, "crop.ini")
    if not importUtils.readCropParameters(cropSettings, params):
        return

    if C3DParameters.computeTranspiration:
        crop.initializeCrop()
    else:
        crop.maxRootFactor = 1.0

    if C3DParameters.computeTranspiration or C3DParameters.computeEvaporation:
        print("Read weather and irrigation data...")
        criteria3D.weatherData = importUtils.readMeteoData(weatherFolder)
        criteria3D.weatherData.set_index(["timestamp"])
        criteria3D.waterData = importUtils.readWaterData(waterFolder, criteria3D.weatherData.iloc[0]["timestamp"],
                                              criteria3D.weatherData.iloc[-1]["timestamp"])
    else:
        print("Read water data...")
        criteria3D.waterData = importUtils.readWaterData(waterFolder, NODATA, NODATA)

    criteria3D.waterData.set_index(["timestamp"])
    # criteria3D.weatherData, criteria3D.waterData = importUtils.transformDates(weatherData, waterData)
    print("Total simulation time [hours]:", len(criteria3D.waterData))

    # initialize export and state
    exportUtils.createExportFile(outputFolder, settingsFolder, params["iteration"])
    modelState = os.path.join(stateFolder, f"modelState{iterations_str}.bin")
    obsState = os.path.join(stateFolder, f"obsState{iterations_str}.csv")

    if C3DParameters.isPeriodicAssimilation or C3DParameters.isFirstAssimilation:
        print("Read observed water potential...")
        obsWaterPotential = pd.read_csv(os.path.join(obsDataFolder, f"waterPotential.csv"))

    # first assimilation
    weatherIndex = 0
    if C3DParameters.isFirstAssimilation:
        print("Assimilate observed water potential (first day)...")
        timeStamp = criteria3D.waterData.loc[weatherIndex]["timestamp"]
        importUtils.writeObsData(obsState, obsWaterPotential, timeStamp)
        importUtils.loadObsData(obsState)

        waterBalance.initializeBalance()
        for i in range(24):
            criteria3D.computeOneHour(weatherIndex+i, False)
        importUtils.loadObsData(obsState)

    waterBalance.initializeBalance()
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))

    if C3DParameters.isVisual:
        visual3D.initialize(1200)
        visual3D.isPause = True
        # wait for start (press 'r')
        while visual3D.isPause:
            time.sleep(0.00001)

    # main cycle
    print("Start...")
    currentIndex = 1
    isFirstRun = True
    restartIndex = weatherIndex
    while weatherIndex < len(criteria3D.waterData):
        criteria3D.computeOneHour(weatherIndex, C3DParameters.isVisual)

        currentTimeStamp = criteria3D.waterData.loc[weatherIndex]["timestamp"]
        currentDateTime = pd.to_datetime(currentTimeStamp, unit='s')

        # assimilation
        if C3DParameters.isPeriodicAssimilation and not C3DParameters.isForecast:
            if (currentIndex % C3DParameters.assimilationInterval) == 0:
                print("Assimilation:", currentDateTime)
                importUtils.writeObsData(obsState, obsWaterPotential, currentTimeStamp)
                importUtils.loadObsData(obsState)

        # save model state
        if C3DParameters.isForecast and currentIndex == C3DParameters.assimilationInterval:
            importUtils.saveCurrentModelState(modelState)
            restartIndex = weatherIndex

        # save output
        if not C3DParameters.isForecast: # or isFirstRun:
            exportUtils.takeScreenshot(currentTimeStamp)
        else:
            if currentIndex > (C3DParameters.forecastPeriod - C3DParameters.assimilationInterval):
                exportUtils.takeScreenshot(currentTimeStamp)

        # daily water balance
        if currentDateTime.hour == 0:
            if currentIndex >= 24:
                lastDate = currentDateTime.date() - datetime.timedelta(days=1)
                exportUtils.saveDailyBalance(lastDate)
            waterBalance.dailyBalance.initialize()

        # restart
        if C3DParameters.isForecast and currentIndex == C3DParameters.forecastPeriod:
            print("Restart...")
            importUtils.loadModelState(modelState)
            # assimilation
            obsWeather = criteria3D.weatherData.loc[restartIndex]
            importUtils.writeObsData(obsState, obsWaterPotential, obsWeather["timestamp"])
            importUtils.loadObsData(obsState)
            # redraw
            waterBalance.totalTime = restartIndex * 3600
            if C3DParameters.isVisual:
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


def parse_args():
    parser = argparse.ArgumentParser(description="CRITERIA")

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=False,
        default=os.path.join("data", "errano_all"),
        help="path to working directory",
    )
    parser.add_argument(
        "-it",
        "--iteration",
        type=int,
        required=False,
        default=-1,
        help="tuning iteration index",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    main(args)