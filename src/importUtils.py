import pandas as pd
import os
import criteria3D
import assimilation
import waterBalance
import tkinter
import tkinter.filedialog


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


def loadState(fileName):
    currentState = pd.read_csv(fileName)
    assimilation.assimilate(currentState)
    waterBalance.updateStorage()
    return True


def loadState_old(fileName):
    if not os.path.isfile(fileName):
        return False

    state = pd.read_csv(fileName)
    pos = []
    potential = []
    for _, position in state.iterrows():
        x = position['x']                   # [m]
        y = position['y']                   # [m]
        depth = position['z']               # [m]
        psi = position['value'] / 9.81      # water potential - from [kPa] to [m]
        pos.append([x, y, depth])
        potential.append(psi)
        # symmetric values
        if x != 0:
            pos.append([-x, y, depth])
            potential.append(psi)

    criteria3D.setModelState(pos, potential)
    return True


def writeState(stateFileName, obsData, timeStamp):
    header = "x,y,z,value\n"
    f = open(stateFileName, "w")
    f.write(header)

    for i in range(len(obsData)):
        currentTime = obsData.iloc[i].at["timestamp"]
        if currentTime == timeStamp:
            y = 0
            z = 0.6
            x = 0
            for j in range(1, len(obsData.columns)):
                psi = obsData.iloc[i].at[obsData.columns[j]]
                row = '{:.2f}'.format(x) + "," + '{:.1f}'.format(y) + "," + '{:.1f}'.format(z) + ","
                row += '{:.1f}'.format(psi) + "\n"
                f.write(row)
                x += 0.25
                if (x == 0.75):
                    x = 0.8
                if (x > 1):
                    x = 0
                    z -= 0.2
            return True
    return False
    