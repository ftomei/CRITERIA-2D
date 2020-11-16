import pandas as pd
import os
from enum import Enum 

TIME_LENGHT = 60

def readArpaeData(arpaePath):
    stationInfoFile = "station_info.csv" 
    stationInfo = pd.read_csv(os.path.join(arpaePath, stationInfoFile))

    humidityFile = "humidity.csv"
    humidity = pd.read_csv(os.path.join(arpaePath, humidityFile))

    precipitationsFile= "precipitations.csv"
    precipitations = pd.read_csv(os.path.join(arpaePath, precipitationsFile))

    radiationsFile= "radiations.csv" 
    radiations = pd.read_csv(os.path.join(arpaePath, radiationsFile))

    temperatureFile = "temperature.csv"
    temperature = pd.read_csv(os.path.join(arpaePath, temperatureFile))

    windFile = "wind.csv"
    wind = pd.read_csv(os.path.join(arpaePath, windFile))


    arpaeData = [humidity, precipitations, radiations, temperature, wind]

    for i in range(len(arpaeData) - 1):
        if arpaeData[i].iloc[0]["start"] != arpaeData[i + 1].iloc[0]["start"]:
            raise Exception("Arpae files have different time spans")

        if arpaeData[i].iloc[0]["end"] != arpaeData[i + 1].iloc[0]["end"]:
            raise Exception("Arpae files haa different time steps")

        if arpaeData[i].iloc[-1]["start"] != arpaeData[i + 1].iloc[-1]["start"]:
            raise Exception("Arpae files have different time spans")

    for i in range(len(arpaeData)):
        arpaeData[i] = arpaeData[i].set_index(['start', 'end'])

    mergedDf = arpaeData[0]
    for i in range(1, len(arpaeData)):
        mergedDf = mergedDf.merge(arpaeData[i], left_index=True, right_index=True)

    return stationInfo, mergedDf.reset_index()

def readIrrigationsData(irrigationsPath, arpae_start, arpae_end):
    irrigationsConfigurationsFile = "irrigations_configurations.csv" 
    irrigationsConfigurations = pd.read_csv(os.path.join(irrigationsPath, irrigationsConfigurationsFile))

    irrigationsFile = "irrigations.csv"
    irrigations = pd.read_csv(os.path.join(irrigationsPath, irrigationsFile))

    if arpae_start != irrigations.iloc[0]["start"] or arpae_end != irrigations.iloc[-1]["start"]:
        raise Exception("Irrigations file has a different time span")

    return irrigationsConfigurations, irrigations