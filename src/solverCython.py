#solverCython.py

from dataStructures import *
from waterProcesses import runoff, infiltration, redistribution
import boundaryConditions
import soil
import waterBalance
from solverC import setArraysC, set_x, set_C, set_A, set_indices, get_x, arrangeMatrix, GaussSeidel


def setCriteria3DArrays(nrCells, nrLinks):
    setArraysC(nrCells, nrLinks)
    
def computeStep(deltaT):
    global x, C, indices
    #initialize
    approximation = 1
    isValidStep = False
    for i in range(C3DStructure.nrCells):
        C3DCells[i].H0 = C3DCells[i].H
        set_x(i, C3DCells[i].H)
        if  (C3DCells[i].isSurface): 
            set_C(i, C3DCells[i].area)
    
    while ((not isValidStep) 
    and (approximation <= C3DParameters.maxApproximationsNr)):
        isFirstApprox = (approximation == 1)
        waterBalance.maxCourant = 0.0
        for i in range(C3DStructure.nrCells):
            if (not C3DCells[i].isSurface):
                C3DCells[i].Se = soil.getDegreeOfSaturation(i)
                C3DCells[i].k = soil.getHydraulicConductivity(i)
                dTheta_dH = soil.getdTheta_dH(i)
                set_C(i, C3DCells[i].volume * dTheta_dH)
        boundaryConditions.updateBoundary(deltaT)
        
        print ("approximation nr:", approximation)
        print ("Sum flows (abs) [l]:", format(waterBalance.sumWaterFlow(deltaT, True) * 1000.,".5f"))  
        
        for i in range(C3DStructure.nrCells):
            k = 0
            if (newMatrixElement(i, C3DCells[i].upLink, k, 
                                 False, deltaT, isFirstApprox)): 
                k += 1
            for l in range(C3DStructure.nrLateralLinks):
                if (newMatrixElement(i, C3DCells[i].lateralLink[l], k,
                                 True, deltaT, isFirstApprox)): 
                    k += 1
            if (newMatrixElement(i, C3DCells[i].downLink, k,
                                 False, deltaT, isFirstApprox)): 
                k += 1
            if (k < C3DStructure.nrMaxLinks): 
                set_indices(i, k, NOLINK)
            
            arrangeMatrix(i, deltaT, C3DCells[i].H0, C3DCells[i].flow)
            
        if ((waterBalance.maxCourant > 1.0) 
        and (deltaT > C3DParameters.deltaT_min)):
            print ("Courant too high:", waterBalance.maxCourant)
            print ("Decrease time step")
            while (waterBalance.maxCourant > 1.0):
                waterBalance.halveTimeStep()
                waterBalance.maxCourant *= 0.5
            return False

        if not solveMatrix(approximation):
            waterBalance.halveTimeStep()
            print("System not convergent.")
            return False
       
        # new hydraulic head
        for i in range(0, C3DStructure.nrCells):
            C3DCells[i].H = get_x(i)
            # check surface error
            if (C3DCells[i].isSurface):
                if (C3DCells[i].H < C3DCells[i].z):
                    C3DCells[i].H = C3DCells[i].z
            C3DCells[i].Se = soil.getDegreeOfSaturation(i)
            
        # waterBalance
        isValidStep = waterBalance.waterBalance(deltaT, approximation)
        if (waterBalance.forceExit): return False
        approximation += 1
        
    if isValidStep:
        for i in range(C3DStructure.nrCells):
            C3DCells[i].Hprev = C3DCells[i].H0
    else: 
        # restoreWater
        for i in range(C3DStructure.nrCells):
            C3DCells[i].H = C3DCells[i].H0
        
    return isValidStep


def  newMatrixElement(i, link, k, isLateral, deltaT, isFirstApprox):
    global A, indices
    j = link.index
    if (j == NOLINK): return False
    
    value = 0.0
    if C3DCells[i].isSurface:
        if C3DCells[j].isSurface:
            if (C3DParameters.computeSurfaceFlux):
                value = runoff(i, link, deltaT, isFirstApprox)
            else:
                value = 0.0
        else:
            value = infiltration(i, j, link, deltaT, isFirstApprox)
    else:
        if C3DCells[j].isSurface:
            value = infiltration(j, i, link, deltaT, isFirstApprox)
        else:
            value = redistribution(i, link, isLateral, deltaT)
    if (value == 0.0): return False
    
    set_indices(i, k, j)
    set_A(i, k, value)
    return True

def solveMatrix(approximation):
    ratio = (C3DParameters.maxIterationsNr / C3DParameters.maxApproximationsNr)
    maxIterationsNr = max(10, ratio * approximation)
    
    iteration = 0
    norm = 1000.
    bestNorm = norm
    while ((norm > C3DParameters.residualTolerance) 
    and (iteration < maxIterationsNr)): 
        norm = GaussSeidel()
        if norm > (bestNorm * 10.): return False
        bestNorm = min(norm, bestNorm)
        iteration += 1
    return True
