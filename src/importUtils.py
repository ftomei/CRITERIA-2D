import pandas as pd
import os
from enum import Enum 

def readArpaeData(arpaePath):
    stationInfoFile = "station_info.csv" 
    stationInfo = pd.read_csv(os.path.join(arpaePath, stationInfoFile))

    humidityFile = "humidity.csv"
    humidity = pd.read_csv(os.path.join(arpaePath, humidityFile))

    radiationsFile= "radiations.csv" 
    radiations = pd.read_csv(os.path.join(arpaePath, radiationsFile))

    temperatureFile = "temperature.csv"
    temperature = pd.read_csv(os.path.join(arpaePath, temperatureFile))

    windFile = "wind.csv"
    wind = pd.read_csv(os.path.join(arpaePath, windFile))

    arpaeData = [humidity, radiations, temperature, wind]

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

def readWaterData(waterPath, arpae_start, arpae_end):
    irrigationsConfigurationsFile = "irrigations_configurations.csv" 
    irrigationsConfigurations = pd.read_csv(os.path.join(waterPath, irrigationsConfigurationsFile))

    irrigationsFile = "irrigations.csv"
    irrigations = pd.read_csv(os.path.join(waterPath, irrigationsFile))

    precipitationsFile= "precipitations.csv"
    precipitations = pd.read_csv(os.path.join(waterPath, precipitationsFile))

    if irrigations.iloc[0]["start"] != precipitations.iloc[0]["start"]:
        raise Exception("Arpae files have different time spans")

    if irrigations.iloc[0]["end"] != precipitations.iloc[0]["end"]:
        raise Exception("Arpae files have different time steps")

    if irrigations.iloc[-1]["start"] != precipitations.iloc[-1]["start"]:
        raise Exception("Arpae files have different time spans")

    irrigations = irrigations.set_index(['start', 'end'])
    precipitations = precipitations.set_index(['start', 'end'])

    mergedDf = irrigations.merge(precipitations, left_index=True, right_index=True)
    mergedDf = mergedDf.reset_index()

    if arpae_start != mergedDf.iloc[0]["start"] or arpae_end != mergedDf.iloc[-1]["start"]:
        raise Exception("Irrigations file has a different time span")

    return irrigationsConfigurations, mergedDf

def transformDates(arpaeData, waterData):
    arpaeData["start"] = pd.to_datetime(arpaeData["start"], infer_datetime_format=True)
    arpaeData["end"] = pd.to_datetime(arpaeData["end"], infer_datetime_format=True)

    waterData["start"] = pd.to_datetime(waterData["start"], infer_datetime_format=True)
    waterData["end"] = pd.to_datetime(waterData["end"], infer_datetime_format=True)
    
    return arpaeData, waterData

def setDataIndeces(arpaeData, waterData):
    arpaeData = arpaeData.set_index(["start"])

    waterData = waterData.set_index(["start"])
    
    return arpaeData, waterData
    