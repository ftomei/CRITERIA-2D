#boundaryConditions.py

from math import sqrt
from dataStructures import *
import soil

def updateBoundary(deltaT):
    retentionCurve = C3DParameters.waterRetentionCurve 
    for i in range(C3DStructure.nrCells):
        # Initialize
        C3DCells[i].flow = 0.0
        C3DCells[i].boundary.flow = 0.0
        
        # sink/source: precipitation, irrigation, evapotranspiration
        if (C3DCells[i].sinkSource != NODATA) and (C3DCells[i].sinkSource != 0):
            #[m3 s-1]
            C3DCells[i].flow = C3DCells[i].sinkSource
  
        if (C3DCells[i].boundary.type != BOUNDARY_NONE):
            slope = C3DCells[i].boundary.slope 
            if (slope == 0.0):
                    slope = 0.001
            meanH = (C3DCells[i].H + C3DCells[i].H0) * 0.5
                 
            if (C3DCells[i].boundary.type == BOUNDARY_RUNOFF):
                if (slope > 0.0):
                    Hs = meanH - (C3DCells[i].z + C3DParameters.pond)
                    if (Hs > EPSILON_METER):
                        boundaryArea = C3DCells[i].boundary.area * Hs
                        maxFlow = (Hs * C3DCells[i].area) / deltaT
                        # Manning equation [m3 s-1]
                        flow = ((boundaryArea / C3DParameters.roughness) * (Hs**(2./3.)) * sqrt(slope))
                        C3DCells[i].boundary.flow = -min(flow, maxFlow)
                    
            elif (C3DCells[i].boundary.type == BOUNDARY_FREELATERALDRAINAGE): 
                k = C3DCells[i].k * C3DParameters.conductivityHVRatio             
                C3DCells[i].boundary.flow = -k * C3DCells[i].boundary.area * slope 
                                      
            elif (C3DCells[i].boundary.type == BOUNDARY_FREEDRAINAGE):
                C3DCells[i].boundary.flow = -C3DCells[i].k * C3DCells[i].upLink.area
            
            elif (C3DCells[i].boundary.type == BOUNDARY_PRESCRIBEDTOTALPOTENTIAL):
                curve = C3DParameters.waterRetentionCurve
                prescribedH = C3DParameters.waterTableDepth
                boundaryPsi = prescribedH - C3DCells[i].z
                boundarySe = soil.degreeOfSaturation(curve, boundaryPsi)
                boundaryK = soil.hydraulicConductivity(curve, boundarySe)
                k = soil.meanK(C3DParameters.meanType, C3DCells[i].k, boundaryK)
                dH = prescribedH - C3DCells[i].H
                C3DCells[i].boundary.flow = k * dH * C3DCells[i].area
                
            C3DCells[i].flow += C3DCells[i].boundary.flow
            
            