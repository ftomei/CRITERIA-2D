#waterProcesses.py

from math import fabs, sqrt
from dataStructures import *
import waterBalance
import soil

CYTHON = True
if CYTHON:
    from solverC import meanK
else:
    from soil import meanK


def redistribution(i, link, isLateral, deltaT):
    j = link.index
    k = meanK(C3DParameters.meanType, C3DCells[i].k, C3DCells[j].k)
    if (isLateral):
        k *= C3DParameters.conductivityHVRatio
        
    return (k * link.area) / link.distance


def infiltration(surf, sub, link, deltaT, isFirstApprox):
    if (C3DCells[surf].z > C3DCells[sub].H):
        #unsaturated
        Havg = (C3DCells[surf].H + C3DCells[surf].H0) * 0.5
        psi = Havg - C3DCells[surf].z
        if isFirstApprox:
            rain = (C3DCells[surf].sinkSource / C3DCells[surf].area) * (deltaT * 0.5)
            psi += rain
        if (psi < EPSILON): return 0.0
        
        interfaceK = meanK(C3DParameters.meanType, C3DCells[sub].k, soil.C3DSoil.Ks)
        dH = C3DCells[surf].H - C3DCells[sub].H
        maxK = (psi / deltaT) * (link.distance / dH)
        k = min(interfaceK , maxK)
    else:
        #saturated
        k = soil.C3DSoil.Ks
    
    return (k  * link.area) / link.distance


def runoff(i, link, deltaT, isFirstApprox):
    j = link.index
    zmax = max(C3DCells[i].z, C3DCells[j].z)
    Hmax = max((C3DCells[i].H + C3DCells[i].H0)/ 2.0, 
               (C3DCells[j].H + C3DCells[j].H0)/ 2.0)
    Hs = Hmax - (zmax + C3DParameters.pond) 
    if (Hs <= EPSILON_METER): return 0.
    
    dH = fabs(C3DCells[i].H - C3DCells[j].H)
    if (dH < EPSILON_METER): return 0.
    
    # pond
    Hs = min(Hs, dH)
    
    # [m/s] Manning equation
    v = (pow(Hs, 2.0 / 3.0) * sqrt(dH/link.distance)) / C3DParameters.roughness
    Courant = v * deltaT / link.distance
    waterBalance.maxCourant = max(waterBalance.maxCourant, Courant)

    # on surface: link.area = side length [m]
    area = link.area * Hs 
    return (v / dH) * area

