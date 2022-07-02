import pandas as pd
import os
import criteria3D
import assimilation
import waterBalance
import tkinter.filedialog
from dataStructures import *
import rectangularMesh
from array import array
from configparser import ConfigParser


def setField(settingsFilename):
    config = ConfigParser()
    config.read(settingsFilename)

    # size
    print("Read field settings...")
    try:
        width = config.getfloat('size', 'width')
    except:
        print("ERROR! Missing size.width in field.ini")
        return False
    height = config.getfloat('size', 'height')
    gridStep = config.getfloat('size', 'gridStep')
    C3DStructure.z0 = config.getfloat('size', 'z0')
    initialize3DStructure(width, height, gridStep)

    # slope
    C3DStructure.slopeX = config.getfloat('slope', 'slopeX')
    C3DStructure.slopeY = config.getfloat('slope', 'slopeY')
    C3DStructure.plantSlope = config.getfloat('slope', 'plantSlope')
    C3DStructure.plantSlopeWidth = config.getfloat('slope', 'plantSlopeWidth')

    # Build rectangular mesh
    print("Building rectangle mesh...")
    rectangularMesh.rectangularMeshCreation()
    print("Nr. of rectangles:", C3DStructure.nrRectangles)
    print("Total area [m^2]:", C3DStructure.totalArea)

    # set plant
    print("set plant positions...")
    xStr = config.get('plant', 'x').split(',')
    x = [float(each) for each in xStr]
    yStr = config.get('plant', 'y').split(',')
    y = [float(each) for each in yStr]
    if len(x) != len(y):
        print("ERROR! Different number of plant.x,y in field.ini")
        return False
    plantIndices.clear()
    for i in range(len(x)):
        surfaceIndex = rectangularMesh.getSurfaceIndex(x[i], y[i])
        if surfaceIndex != NODATA:
            plantIndices.append(surfaceIndex)

    # set dripper
    print("set dripper positions...")
    xStr = config.get('dripper', 'x').split(',')
    x = [float(each) for each in xStr]
    yStr = config.get('dripper', 'y').split(',')
    y = [float(each) for each in yStr]
    if len(x) != len(y):
        print("ERROR: different number of dripper x,y positions.")
        return False
    dripperIndices.clear()
    for i in range(len(x)):
        surfaceIndex = rectangularMesh.getSurfaceIndex(x[i], y[i])
        if surfaceIndex != NODATA:
            dripperIndices.append(surfaceIndex)

    return True


def getStateFileName(isSave):
    root = tkinter.Tk()
    options = {'defaultextension': ".csv", 'filetypes': [("Comma separated values", ".csv")], 'initialdir': "data"}
    if isSave:
        fileName = tkinter.filedialog.asksaveasfilename(**options)
    else:
        fileName = tkinter.filedialog.askopenfilename(**options)
    root.destroy()
    return fileName


def readWaterTable(waterPath):
    date = []
    depth = []
    waterTableData = pd.read_csv(os.path.join(waterPath, "watertable.csv"))
    for _, waterTable in waterTableData.iterrows():
        date.append(pd.to_datetime(waterTable['date']))
        depth.append(waterTable['depth'])
    return date, depth


def readMeteoData(meteoPath):
    stationInfoFile = "station_info.csv"
    stationInfo = pd.read_csv(os.path.join(meteoPath, stationInfoFile))

    humidityFile = "air_humidity.csv"
    humidity = pd.read_csv(os.path.join(meteoPath, humidityFile))

    radiationsFile = "solar_radiation.csv"
    radiations = pd.read_csv(os.path.join(meteoPath, radiationsFile))

    temperatureFile = "air_temperature.csv"
    temperature = pd.read_csv(os.path.join(meteoPath, temperatureFile))

    windFile = "wind_speed.csv"
    wind = pd.read_csv(os.path.join(meteoPath, windFile))

    meteoData = [humidity, radiations, temperature, wind]

    for i in range(len(meteoData) - 1):
        if meteoData[i].iloc[0]["timestamp"] != meteoData[i + 1].iloc[0]["timestamp"]:
            raise Exception("Meteo files have different time spans")

    for i in range(len(meteoData)):
        meteoData[i] = meteoData[i].set_index(['timestamp'])

    mergedDf = meteoData[0]
    for i in range(1, len(meteoData)):
        mergedDf = mergedDf.merge(meteoData[i], how='outer', left_index=True, right_index=True)

    return stationInfo, mergedDf.reset_index()


def readWaterData(waterPath, meteo_start, meteo_end):
    irrigationFile = "irrigation.csv"
    irrigation = pd.read_csv(os.path.join(waterPath, irrigationFile))

    precipitationsFile = "precipitation.csv"
    precipitations = pd.read_csv(os.path.join(waterPath, precipitationsFile))

    if irrigation.iloc[0]["timestamp"] != precipitations.iloc[0]["timestamp"]:
        raise Exception("Water files have different time spans")

    if irrigation.iloc[-1]["timestamp"] != precipitations.iloc[-1]["timestamp"]:
        raise Exception("Water files have different time spans")

    irrigation = irrigation.set_index(['timestamp'])
    precipitations = precipitations.set_index(['timestamp'])

    mergedDf = irrigation.merge(precipitations, how='outer', left_index=True, right_index=True)
    mergedDf = mergedDf.reset_index()

    if meteo_start != mergedDf.iloc[0]["timestamp"] or meteo_end != mergedDf.iloc[-1]["timestamp"]:
        raise Exception("Irrigation file has a different time span")

    return mergedDf


def transformDates(meteoData, waterData):
    meteoData["time"] = pd.to_datetime(meteoData["timestamp"], infer_datetime_format=True)
    waterData["time"] = pd.to_datetime(waterData["timestamp"], infer_datetime_format=True)
    return meteoData, waterData


def loadObsState(fileName):
    obsState = pd.read_csv(fileName)
    assimilation.assimilate(obsState)
    waterBalance.updateStorage()
    return True


def writeObsState(stateFileName, obsData, timeStamp):
    header = "x,y,z,value\n"
    f = open(stateFileName, "w")
    f.write(header)

    df = obsData[obsData["timestamp"] == timeStamp]
    if not df.empty:
        for column in [column for column in list(df.columns) if column != "timestamp"]:
            splitted_column = column.split("_")
            f.write(
                "{:.2f}".format(float(splitted_column[2][1:]) / 100) + ","  # x
                + "{:.1f}".format(float(splitted_column[1][1:]) / 100) + ","  # y
                + "{:.1f}".format(float(splitted_column[0][1:]) / 100) + ","  # z
                + "{:.1f}".format(df[column].values[0]) + "\n"  # psi
            )
    return not df.empty


def saveCurrentModelState(stateFileName):
    float_array = array('d')
    for i in range(len(C3DCells)):
        float_array.append(C3DCells[i].H)

    f = open(stateFileName, "wb")
    float_array.tofile(f)
    f.close()


def loadModelState(stateFileName):
    f = open(stateFileName, "rb")
    float_array = array('d')
    float_array.fromfile(f, len(C3DCells))
    f.close()

    for i in range(len(C3DCells)):
        H = float_array[i]
        signPsi = H - C3DCells[i].z
        criteria3D.setMatricPotential(i, signPsi)

    waterBalance.updateStorage()
