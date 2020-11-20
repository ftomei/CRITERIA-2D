# PenmanMonteith.py

# Penman Monteith equation for potential evapotranspiration
# from Allen et al., 1998 
# Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements

import math
import transmissivity
import datetime
import pandas as pd
from copy import copy

# [W m-2 K-4] Stefan-Boltzmann constant
STEFAN_BOLTZMANN = 5.670373E-8      
# [Pa] standard atmospheric pressure at sea level
P0 = 101300.
# [K m-1] constant lapse rate of moist air
LAPSE_RATE_MOIST_AIR = 0.0065
# [J kg-1 K-1] specific gas constant for dry air
R_DRY_AIR = 287.058
# [J kg-1 K-1] specific heat at constant pressure
CP = 1013.  
# [m s-2] gravity acceleration
GRAVITY = 9.80665
# [K]
ZEROCELSIUS = 273.15
# [K] temperature at reference pressure level (P0)
TEMP_P0 = 20. + ZEROCELSIUS
# [-] albedo of reference crop (grass)
ALBEDO_CROP_REFERENCE = 0.23     
# [-] ratio molecular weight of water vapour/dry air
RATIO_WATER_VD = 0.622

TRASMISSIVITY_THREESHOLD = 300

# pressure [Pa]
# height [m]
def pressureFromAltitude(height):
    return P0 * math.pow(1. + height * LAPSE_RATE_MOIST_AIR / TEMP_P0, 
                         - GRAVITY / (LAPSE_RATE_MOIST_AIR * R_DRY_AIR))

# saturation vapor pressure [Pa]
# airTemperature [degC] 
def saturationVaporPressure(airTemperature):
    return 611 * math.exp(17.502 * airTemperature / (airTemperature + 240.97))

def computeNormTransmissivity(arpaeRelevation, arpaeData, latitude, longitude):
    currentDate = copy(arpaeRelevation["end"])
    potentialRad = 0
    observedRad = 0
    nrHoursAhead = -1
    
    while potentialRad < TRASMISSIVITY_THREESHOLD:
        nrHoursAhead += 1
        date = datetime.date(currentDate.year, currentDate.month, currentDate.day)
        hour = currentDate.hour
        potentialRad += transmissivity.clearSkyRad(date, hour, latitude, longitude)
        observedRad += arpaeData.loc[currentDate - pd.Timedelta('1 hour')]["radiations"]
        currentDate += pd.Timedelta('1 hour')

    currentDate = copy(arpaeRelevation["end"])
    for i in range(nrHoursAhead):
        currentDate -= pd.Timedelta('1 hour')
        date = datetime.date(currentDate.year, currentDate.month, currentDate.day)
        hour = currentDate.hour
        potentialRad += transmissivity.clearSkyRad(date, hour, latitude, longitude)
        observedRad += arpaeData.loc[currentDate - pd.Timedelta('1 hour')]["radiations"]
    
    return observedRad / potentialRad

# return reference evapotranspiration (mm)
# height               elevation above mean sea level (meters)
# airTemperature       air temperature (C)
# globalSWRadiation    global Short Wave radiation (W m-2)
# airRelHumidity       air relative humidity (%)
# windSpeed_10m        wind speed at 10 meters (m s-1)
# normTransmissivity   normalized transmissivity [0-1]
def computeHourlyET0(height, airTemperature, globalSWRadiation, airRelHumidity, 
                     windSpeed_10m, normTransmissivity):   
    # air temperature [Kelvin]  
    airTempKelvin = airTemperature + ZEROCELSIUS
    # wind speed at 2 meters [m s-1]
    windSpeed_2m = windSpeed_10m * 0.748
    # barometric pressure [kPa]
    pressure = pressureFromAltitude(height) / 1000.
    # saturation vapor pressure [kPa] 
    satVapPressure = saturationVaporPressure(airTemperature) / 1000.
    # current vapor pressure [kPa] 
    vaporPressure = satVapPressure * (airRelHumidity / 100.)
    # net emissivity of the surface [-]
    emissivity = 0.34 - 0.14 * math.sqrt(vaporPressure)
    # cloudiness factor for long wave radiation [-]
    cloudFactor = max(0, 1.35 * min(normTransmissivity, 1) - 0.35)
    # net longwave radiation [J m-2 s-1]
    netLWRadiation = cloudFactor * emissivity * STEFAN_BOLTZMANN * math.pow(airTempKelvin, 4.)
    # net radiation [J m-2 s-1] 
    netRadiation = (1. - ALBEDO_CROP_REFERENCE) * globalSWRadiation - netLWRadiation
    
    # from [W m-2] to [J m-2 h-1]
    netRadiation = netRadiation * 3600.
    
    # values for grass
    # g   soil heat flux density [J m-2 h-1]  
    # Cd  bulk surface resistance and aerodynamic resistance coefficient  
    if (netRadiation > 0):
        g = 0.1 * netRadiation
        Cd = 0.24
    else:
        g = 0.5 * netRadiation
        Cd = 0.96
              
    # slope of saturation vapor pressure curve [kPa K]
    slope = 4098. * satVapPressure / (airTempKelvin * airTempKelvin)
    # latent heat of vaporization [J kg-1]
    latentHeatVap = 2501000. - 2369.2 * airTemperature;
    # psychrometric instrument constant [kPa K-1] 
    psychro = CP * pressure / (RATIO_WATER_VD * latentHeatVap);

    denominator = slope + psychro * (1. + Cd * windSpeed_2m)
    firstTerm = slope * (netRadiation - g) / (latentHeatVap * denominator)
    secondTerm = (psychro * (37. / airTempKelvin) * windSpeed_2m * (satVapPressure - vaporPressure)) / denominator

    return max(firstTerm + secondTerm, 0.)



