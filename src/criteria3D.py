# criteria3D.py

from math import fabs
from dataStructures import *
import rectangularMesh
import waterBalance
import visual3D
import soil
import time

if CYTHON:
    import solverCython as solver
else:
    import solver as solver

irrigationIndices = []
plantIndices = []


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
        x = position['x']
        y = position['y']
        surfaceIndex = rectangularMesh.getSurfaceIndex(x, y)
        if surfaceIndex != NODATA:
            irrigationIndices.append(surfaceIndex)


def setPlantPositions(plantConfigurations):
    for _, position in plantConfigurations.iterrows():
        x = position['plant_x']
        y = position['plant_y']
        surfaceIndex = rectangularMesh.getSurfaceIndex(x, y)
        if surfaceIndex != NODATA:
            plantIndices.append(surfaceIndex)


def getCellDistance(i, j):
    v1 = [C3DCells[i].x, C3DCells[i].y, C3DCells[i].z]
    v2 = [C3DCells[j].x, C3DCells[j].y, C3DCells[j].z]
    return rectangularMesh.distance3D(v1, v2)


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
        C3DCells[i].k = soil.horizon.Ks
    else:
        C3DCells[i].H = C3DCells[i].z + signPsi
        C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        C3DCells[i].k = soil.getHydraulicConductivity(i)
    C3DCells[i].H0 = C3DCells[i].H
    return OK


def initializeSinkSource(cellType):
    if cellType == ONLY_SURFACE:
        for i in range(C3DStructure.nrRectangles):
            C3DCells[i].sinkSource = 0
    else:
        for i in range(C3DStructure.nrCells):
            C3DCells[i].sinkSource = 0


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
    if rate > 0:
        C3DParameters.pond = C3DParameters.pondIrrigation
    else:
        C3DParameters.pond = C3DParameters.pondRainfall

    for index in irrigationIndices:
        C3DCells[index].sinkSource += rate * 0.001  # [m^3 s^-1]


def setModelState(position, psiValues):
    for surfaceIndex in range(C3DStructure.nrRectangles):
        for layer in range(C3DStructure.nrLayers):
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            # surface
            if layer == 0:
                setMatricPotential(i, 0)
            # subsurface
            else:
                p0 = rectangularMesh.getXYDepth(i)
                # inverse distance weight
                sumWeights = 0
                sumPsi = 0
                for index in range(len(position)):
                    distance = rectangularMesh.distance3D(position[index], p0)
                    distance = max(distance, 0.001)
                    weight = 1. / (distance * distance * distance)
                    sumWeights += weight
                    sumPsi += psiValues[index] * weight

                psi = (sumPsi / sumWeights)
                setMatricPotential(i, psi)

    waterBalance.initializeBalance()


# timeLength        [s]          
def compute(timeLength, isRedraw):
    currentTime = 0
    while currentTime < timeLength:
        residualTime = timeLength - currentTime
        deltaT = min(C3DParameters.currentDeltaT, residualTime)
        acceptedStep = False

        while not acceptedStep:
            if isRedraw:
                if visual3D.isPause:
                    # print("\nPress 'r' to run")
                    while visual3D.isPause:
                        time.sleep(0.00001)

            deltaT = min(C3DParameters.currentDeltaT, residualTime)
            # print("time step [s]: ", deltaT)
            # print("sink/source [l]:", format(waterBalance.sumSinkSource(deltaT) * 1000., ".5f"))

            acceptedStep = solver.computeStep(deltaT)
            if not acceptedStep:
                # restoreWater
                for i in range(C3DStructure.nrCells):
                    C3DCells[i].H = C3DCells[i].H0

        if isRedraw:
            visual3D.redraw()
        currentTime += deltaT
