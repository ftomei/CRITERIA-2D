#transmissivity.py

import math
import datetime
import pandas as pd

SOLAR_CONSTANT = 1367           # [W/m²]
TRASMISSIVITY_THRESHOLD = 300   # [W/m²]
TIMEZONE = 1
MAXIMUM_TRANSMISSIVITY = 0.75


def degreeToRadians(degree):
    return degree * (math.pi / 180.)

def dateToDOY(year, month, day):
    myDate = datetime.date(year, month, day)
    januaryFirst = datetime.date(year, 1, 1)
    return (myDate - januaryFirst).days + 1
       
#return clear sky solar radiation 
def clearSkyRad(myDate, finalHourUTC, latDegrees, lonDegrees):
    latRad = degreeToRadians(latDegrees)
    longitudeCorrection = ((TIMEZONE * 15.) - lonDegrees) / 15.
    
    solarTime = finalHourUTC - 0.5 + TIMEZONE
    if (solarTime < 0):
        solarTime += 24
        myDate = myDate - 1
    if (solarTime > 24):
        solarTime -= 24
        myDate = myDate -1
    
    doy = dateToDOY(myDate.year, myDate.month, myDate.day)
    timeAdjustment = degreeToRadians(279.575 + 0.986 * doy)
    
    timeEq = (-104.7 * math.sin(timeAdjustment) + 596.2 * math.sin(2. * timeAdjustment) 
              + 4.3 * math.sin(3. * timeAdjustment) - 12.7 * math.sin(4 * timeAdjustment) 
              - 429.3 * math.cos(timeAdjustment) - 2. * math.cos(2. * timeAdjustment) 
              + 19.3 * math.cos(3. * timeAdjustment)) / 3600.
              
    solarNoon = 12. + longitudeCorrection - timeEq
    
    solarDeclination = 0.4102 * math.sin(2. * math.pi / 365. * (doy - 80.))
    
    solarAngle = math.asin(math.sin(latRad) * math.sin(solarDeclination)
                        + math.cos(latRad) * math.cos(solarDeclination)
                        * math.cos(math.pi / 12 * (solarTime - solarNoon)))
    
    return max(0, SOLAR_CONSTANT * math.sin(solarAngle))


def computeNormTransmissivity(obsData, currentDateTime, latitude, longitude):
    potentialRad = 0
    observedRad = 0
    nrHoursAhead = 0
    
    myDateTime = currentDateTime
    while potentialRad < TRASMISSIVITY_THRESHOLD:
        date = datetime.date(myDateTime.year, myDateTime.month, myDateTime.day)
        hour = myDateTime.hour
        potentialRad += clearSkyRad(date, hour, latitude, longitude)
        observedRad += obsData.loc[myDateTime - pd.Timedelta('1 hour')]["radiations"]
        myDateTime += pd.Timedelta('1 hour')
        if (potentialRad >= TRASMISSIVITY_THRESHOLD): 
            nrHoursAhead += 1

    myDateTime = currentDateTime
    for i in range(nrHoursAhead):
        myDateTime -= pd.Timedelta('1 hour')
        date = datetime.date(myDateTime.year, myDateTime.month, myDateTime.day)
        hour = myDateTime.hour
        potentialRad += clearSkyRad(date, hour, latitude, longitude)
        observedRad += obsData.loc[myDateTime - pd.Timedelta('1 hour')]["radiations"]
    
    return min(1, observedRad / (potentialRad * MAXIMUM_TRANSMISSIVITY))
    
    
    