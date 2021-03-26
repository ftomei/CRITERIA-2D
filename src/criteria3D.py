# criteria3D.py

from math import fabs
from dataStructures import *
from rectangularMesh import distance3D
import waterBalance
#import visual3D
import soil
import time
import crop

CYTHON = True
if CYTHON:
    import solverCython as solver
else:
    import solver as solver

irrigationIndeces = []


def memoryAllocation(nrLayers, nrRectangles):
    C3DStructure.nrRectangles = nrRectangles
    C3DStructure.nrLayers = nrLayers
    C3DStructure.nrCells = nrLayers * nrRectangles
    solver.setCriteria3DArrays(C3DStructure.nrCells, C3DStructure.nrMaxLinks)
    for i in range(C3DStructure.nrCells):
        C3DCells.append(Ccell())


def setCellGeometry(i, x, y, z, volume, area):
    C3DCells[i].x = x
    C3DCells[i].y = y
    C3DCells[i].z = z
    C3DCells[i].volume = volume
    C3DCells[i].area = area


def setCellProperties(i, isSurface, boundaryType):
    C3DCells[i].isSurface = isSurface
    C3DCells[i].boundary.type = boundaryType


def setBoundaryProperties(i, area, slope):
    C3DCells[i].boundary.area = area
    C3DCells[i].boundary.slope = slope


def setDripIrrigationPositions(irrigationConfigurations):
    for _, position in irrigationConfigurations.iterrows():
        # x
        xOffset = position['x'] / C3DStructure.gridStep
        if xOffset < 0:
            xOffset = 0
        if xOffset >= C3DStructure.nrRectanglesInXAxis:
            xOffset = C3DStructure.nrRectanglesInXAxis - 1
        # y
        yOffset = position['y'] / C3DStructure.gridStep
        if yOffset < 0:
            yOffset = 0
        if yOffset >= C3DStructure.nrRectanglesInYAxis:
            yOffset = C3DStructure.nrRectanglesInYAxis - 1
        # index
        index = int(C3DStructure.nrRectanglesInXAxis * yOffset + xOffset)
        irrigationIndeces.append(index)


def getCellDistance(i, j):
    v1 = [C3DCells[i].x, C3DCells[i].y, C3DCells[i].z]
    v2 = [C3DCells[j].x, C3DCells[j].y, C3DCells[j].z]
    return distance3D(v1, v2)


# -----------------------------------------------------------
# direction:         UP, DOWN, LATERAL
# interfaceArea      [m^2]            
# -----------------------------------------------------------
def SetCellLink(i, linkIndex, direction, interfaceArea):
    if direction == UP:
        C3DCells[i].upLink.index = linkIndex
        C3DCells[i].upLink.area = interfaceArea
        C3DCells[i].upLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return OK
    elif direction == DOWN:
        C3DCells[i].downLink.index = linkIndex
        C3DCells[i].downLink.area = interfaceArea
        C3DCells[i].downLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return OK
    elif direction == LATERAL:
        for j in range(C3DStructure.nrLateralLinks):
            if C3DCells[i].lateralLink[j].index == NOLINK:
                C3DCells[i].lateralLink[j].index = linkIndex
                C3DCells[i].lateralLink[j].area = interfaceArea
                C3DCells[i].lateralLink[j].distance = getCellDistance(i, linkIndex)
                return OK
    else:
        return LINK_ERROR


def setMatricPotential(i, signPsi):
    if C3DCells[i].isSurface:
        C3DCells[i].H = C3DCells[i].z + max(signPsi, 0.0)
        C3DCells[i].Se = 1.
        C3DCells[i].k = soil.C3DSoil.Ks
    else:
        C3DCells[i].H = C3DCells[i].z + signPsi
        C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        C3DCells[i].k = soil.getHydraulicConductivity(i)
    C3DCells[i].H0 = C3DCells[i].H
    return OK


# -----------------------------------------------------------
# set uniform rainfall rate
# rain            [mm]
# duration        [s]            
# -----------------------------------------------------------
def setRainfall(rain, duration):
    rate = (rain * 0.001) / duration  # [m s^-1]
    for i in range(C3DStructure.nrRectangles):
        area = C3DCells[i].area  # [m^2]
        C3DCells[i].sinkSource += rate * area  # [m^3 s^-1]


# -----------------------------------------------------------
# set drip irrigation
# irrigation      [l]
# duration        [s]            
# -----------------------------------------------------------
def setDripIrrigation(irrigation, duration):
    rate = irrigation / duration  # [l s^-1]
    for index in irrigationIndeces:
        C3DCells[index].sinkSource += rate * 0.001  # [m^3 s^-1]


# timeLength        [s]          
def compute(timeLength):
    currentTime = 0
    while currentTime < timeLength:
        residualTime = timeLength - currentTime
        acceptedStep = False
        while not acceptedStep:
            #if visual3D.isPause:
            #    print("\nPress 'r' to run")
            #while visual3D.isPause:
            #    time.sleep(0.00001)

            deltaT = min(C3DParameters.currentDeltaT, residualTime)
            # print("\ntime step [s]: ", deltaT)
            # print("sink/source [l]:", format(waterBalance.sumSinkSource(deltaT) * 1000., ".5f"))

            acceptedStep = solver.computeStep(deltaT)
            if not acceptedStep:
                # restoreWater
                for i in range(C3DStructure.nrCells):
                    C3DCells[i].H = C3DCells[i].H0

        #visual3D.redraw()
        currentTime += deltaT
