import pandas as pd
import numpy as np
import os
from PenmanMonteith import computeHourlyET0
from dataStructures import *
from transmissivity import computeNormTransmissivity
from datetime import datetime
import criteria3D
from exportUtils import takeInterpolatedSlice

UNINITIALIZED = 0
UNIXDAY = 86400
UNIXHOUR = 3600
TIMEZONE = 2


# Scheduler params
DRIPPER_CAPACITY = 4  # l/h
MAX_IRRIGATION = 15
first_run = True
startingTimestamp = 1718874000

# Field parms
IRRIGATION_HOUR = 11
IRRIGATION_DAYS_INTERVAL = 7

# Log vars
irrigationDataFrame = pd.DataFrame(columns=["timestamp", "irrigation"])


def predictEt0(weatherData, weatherIndex, timespan=(UNIXDAY * 7)):
    et0Params = weatherData.iloc[weatherIndex - timespan : weatherIndex]
    et0 = 0
    for _, row in et0Params.iterrows():
        et0 += computeHourlyET0(
            C3DStructure.z,
            row["air_temperature"],
            row["solar_radiation"],
            row["air_humidity"],
            row["wind_speed"],
            computeNormTransmissivity(
                weatherData, weatherIndex, C3DStructure.latitude, C3DStructure.longitude
            ),
        )
        weatherIndex += 1
    return et0


def cbar_to_meters(cbar_value):
    return cbar_value * 10.2 / 100


def computeHumidityBins(interpolated_data, settingsFolder):
    # print(len(interpolated_data))
    bins = {
        "-1000000:-10000": 0,
        "-10000:-1500": 0,
        "-1500:-300": 0,
        "-300:-100": 0,
        "-100:-30": 0,
        "-30:0": 0,
        "nCells": 0,
    }
    for cell_index in takeInterpolatedSlice(settingsFolder):
        cellMoisture = interpolated_data[cell_index].H
        if cellMoisture >= cbar_to_meters(-1000000) and cellMoisture < cbar_to_meters(
            -10000
        ):
            bins["-1000000:-10000"] = bins["-1000000:-10000"] + 1
        elif cellMoisture >= cbar_to_meters(-10000) and cellMoisture < cbar_to_meters(
            -1500
        ):
            bins["-10000:-1500"] = bins["-10000:-1500"] + 1
        elif cellMoisture >= cbar_to_meters(-1500) and cellMoisture < cbar_to_meters(
            -300
        ):
            bins["-1500:-300"] = bins["-1500:-300"] + 1
        elif cellMoisture >= cbar_to_meters(-300) and cellMoisture < cbar_to_meters(
            -100
        ):
            bins["-300:-100"] = bins["-300:-100"] + 1
        elif cellMoisture >= cbar_to_meters(-100) and cellMoisture < cbar_to_meters(
            -30
        ):
            bins["-100:-30"] = bins["-100:-30"] + 1
        elif cellMoisture >= cbar_to_meters(-30) and cellMoisture < cbar_to_meters(0):
            bins["-30:0"] = bins["-30:0"] + 1
        bins["nCells"] = bins["nCells"] + 1
    return bins


def applyWateringRules(et0, humidityBins):
    irrigationQuantity = et0
    # if irrigationQuantity > MAX_IRRIGATION:
    #     irrigationQuantity = MAX_IRRIGATION

    if (humidityBins["-30:0"] + humidityBins["-100:-30"]) / humidityBins[
        "nCells"
    ] < 0.7 and humidityBins["-30:0"] / humidityBins["nCells"] < 0.50:
        return irrigationQuantity
    else:
        return 0


def computeDailyIrrigation(cumulative_et0, C3DCells, settingsFolder):
    irrigation = []

    # humidityBins = computeHumidityBins(C3DCells, settingsFolder)
    # irrigationQuantity = applyWateringRules(cumulative_et0, humidityBins)

    irrigationQuantity = cumulative_et0

    print("Next week will be irrigating " + str(irrigationQuantity) + "l")
    for _ in range(0, int(irrigationQuantity / DRIPPER_CAPACITY)):
        irrigation.append(DRIPPER_CAPACITY)
    irrigation.append(
        round(irrigationQuantity - (len(irrigation) * DRIPPER_CAPACITY), 3)
    )
    return irrigation


def scheduleIrrigation(cumulative_et0, weatherIndex, C3DCells, waterFolder, settingsFolder):
    global startingTimestamp
    global irrigationDataFrame
    currentTimeStamp = criteria3D.weatherData.loc[weatherIndex]["timestamp"]
    if (currentTimeStamp - startingTimestamp) < (UNIXDAY * IRRIGATION_DAYS_INTERVAL):
        return False
    startingTimestamp = currentTimeStamp

    plannedIrrigations = pd.DataFrame(columns=["timestamp", "irrigation"])
    plannedIrrigations["timestamp"] = criteria3D.weatherData.iloc[
        weatherIndex : weatherIndex + (24 * IRRIGATION_DAYS_INTERVAL)
    ][
        "timestamp"
    ]  # .apply(lambda x : x - UNIXDAY)
    irrigation = computeDailyIrrigation(cumulative_et0, C3DCells, settingsFolder)
    irrigationHours = True
    irrigationDataFrame = pd.concat(
        [
            irrigationDataFrame,
            pd.DataFrame(
                {"timestamp": [currentTimeStamp], "irrigation": [sum(irrigation)]}
            ),
        ],
        ignore_index=True,
    )
    irrigationDataFrame.to_csv(
        os.path.join(waterFolder, "planned_irrigations.csv"), index=False
    )

    for index, row in plannedIrrigations.iterrows():
        if irrigationHours:
            hourIrrigation = irrigation.pop(0)
            plannedIrrigations.loc[index, "irrigation"] = hourIrrigation
            if len(irrigation) == 0:
                irrigationHours = False
        else:
            plannedIrrigations.loc[index, "irrigation"] = 0

    with open(os.path.join(waterFolder, "irrigation.csv"), "a") as fw:
        for row in plannedIrrigations.values:
            criteria3D.waterData.loc[
                criteria3D.waterData["timestamp"] == int(row[0]), "irrigation"
            ] = float(row[1])
            fw.write(str(row[0]) + "," + str(row[1]) + "\n")
    return True
