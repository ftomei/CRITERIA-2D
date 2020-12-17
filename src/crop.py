#crop.py
from dataStructures import  *
import math
import soil
import numpy as np


#default: kiwifruit
class Ccrop:
    laiMin = 1.0                        # [m2 m-2]
    laiMax = 3.5                        # [m2 m-2]
    rootDepthZero = 0.1                 # [m]
    rootDepthMax = 0.9                  # [m]
    rootDeformation = 1.                # [-]
    kcMax = 1.1                         # [-]
    fRAW = 0.55                         # [-]
    currentLAI = NODATA              
    currentRootDepth = NODATA
    currentRootLenght = NODATA
  
  
def setGrass():
    grass = Ccrop()
    grass.laiMin = 0.5                   # [m2 m-2]
    grass.laiMax = 2.0                   # [m2 m-2]
    grass.rootDepthZero = 0.02           # [m]
    grass.rootDepthMax = 0.65            # [m]
    grass.rootDeformation = 1.0          # [-]
    grass.kcMax = 0.8                    # [-]
    grass.fRAW = 0.65                    # [-]
    return grass 


def setMaxValues(crop): 
    crop.currentLAI = crop.laiMax     
    crop.currentRootDepth = crop.rootDepthMax
    crop.currentRootLenght = crop.currentRootDepth - crop.rootDepthZero   

    
def getCropSurfaceCover(currentLAI):
    k = 0.6      # [-] light extinction coefficient
    
    if (currentLAI <= 0):
        return 0.
    else:
        return 1. - math.exp(-k * currentLAI)


def getMaxEvaporation(currentLAI, ET0):
    maxEvapRatio = 0.66
    
    cropSurfaceCover = getCropSurfaceCover(currentLAI)
    return ET0 * maxEvapRatio * (1. - cropSurfaceCover)


def getMaxTranspiration(currentLAI, kcMax, ET0):
    if (currentLAI <= 0):
        return 0.

    cropSurfaceCover = getCropSurfaceCover(currentLAI)
    kcMaxFactor = 1. + (kcMax - 1.) * cropSurfaceCover
    kc = cropSurfaceCover * kcMaxFactor
    return ET0 * kc 


def cardioidDistribution(deformationFactor, nrLayersWithRoot):
    # initialize
    halfCircle = np.zeros(nrLayersWithRoot)
    cardiod = np.zeros(nrLayersWithRoot*2)
    
    for i in range(nrLayersWithRoot):
        sinAlfa = 1.0 - float(i+1.0) / float(nrLayersWithRoot)
        cosAlfa = max(math.sqrt(1.0 - math.pow(sinAlfa, 2.0)), EPSILON)
        alfa = math.atan(sinAlfa / cosAlfa)
        halfCircle[i] = ((math.pi/2.0) - alfa - sinAlfa*cosAlfa) / math.pi

    lastLayer = nrLayersWithRoot * 2 - 1
    cardiod[0] = halfCircle[0]
    cardiod[lastLayer] = cardiod[0]
    for i in range(1, nrLayersWithRoot):
        cardiod[i] = halfCircle[i] - halfCircle[i-1]
        cardiod[lastLayer-i] = cardiod[i]

    # cardioid deformation
    LiMin = -math.log(0.2) / float(nrLayersWithRoot)
    Limax = -math.log(0.05) / float(nrLayersWithRoot)
    k = LiMin + (Limax - LiMin) * (deformationFactor - 1.0)

    rootDensitySum = 0 ;
    for i in range(nrLayersWithRoot*2):
        cardiod[i] *= math.exp(-k*(i+0.5))
        rootDensitySum += cardiod[i]
        
    # normalize
    for i in range(nrLayersWithRoot*2):
        cardiod[i] /= rootDensitySum
    
    # assign layer density
    layerDensity = np.zeros(nrLayersWithRoot, np.float64)
    for i in range(nrLayersWithRoot):
        layerDensity[i] = cardiod[2*i] + cardiod[2*i+1];

    return layerDensity


def computeRootDensity(crop, nrLayers):
    # Initialize
    rootDensity = np.zeros(nrLayers, np.float64)
    if (crop.currentRootLenght <= 0):
        return rootDensity

    # smallest unit of computation (1 mm)
    atoms = np.zeros(nrLayers, np.int16)
    for i in range(nrLayers):
        atoms[i] = int(round(soil.thickness[i] * 1000))
              
    nrUnrootedAtoms = int(round(crop.rootDepthZero * 1000))
    nrRootedAtoms = int(round(crop.currentRootLenght * 1000))
    densityAtoms = cardioidDistribution(crop.rootDeformation, nrRootedAtoms)
    
    # assign root density
    counter = 0
    for layer in range(nrLayers):
        for i in range(atoms[layer]):
            if (counter >= nrUnrootedAtoms) and (counter < nrUnrootedAtoms + nrRootedAtoms):
                rootDensity[layer] += densityAtoms[counter - nrUnrootedAtoms]
            counter+=1

    # check (rootDensitySum == 1)
    rootDensitySum = 0.0
    for i in range(nrLayers):
        rootDensitySum += rootDensity[i]
    #print ("Sum of root density:", rootDensitySum)

    return rootDensity


def setTranspiration(index, crop, rootDensity, maxTranspiration, timeStep):
    if (maxTranspiration < EPSILON):
        return 
    
    FC = soil.getFieldCapacityWC()                  # [m3 m-3] water content at field capacity
    WP = soil.getWiltingPointWC()                   # [m3 m-3] water content at wilting point
    WSThreshold = FC - crop.fRAW * (FC - WP);       # [m3 m-3] water scarcity stress threshold
    
    # Initialize
    rootDensityWithoutStress = 0.0                  # [-]
    actualTranspiration = 0.0                       # [mm]
    
    nrLayers = len(rootDensity)
    isLayerStressed = np.zeros(nrLayers, dtype=bool)
    layerTranspiration = np.zeros(nrLayers, np.float64)
    for i in range(nrLayers):
        isLayerStressed[i] = False
        layerTranspiration[i] = 0
    
    for layer in range(nrLayers):
        if (rootDensity[layer] > 0):
            i = index + C3DStructure.nrRectangles * layer
            theta = soil.getVolumetricWaterContent(i)
            # TODO water surplus
            
            # WATER SCARSITY
            if (theta < WSThreshold):
                if (theta <= WP):
                    layerTranspiration[layer] = 0
                else:
                    layerTranspiration[layer] = maxTranspiration * rootDensity[layer] * ((theta - WP) / (WSThreshold - WP))
                isLayerStressed[layer] = True
            else:
                # normal conditions
                layerTranspiration[layer] = maxTranspiration * rootDensity[layer]
                
                # check stress
                theta_mm = theta * soil.thickness[layer] * 1000.
                WSThreshold_mm = WSThreshold * soil.thickness[layer] * 1000.
    
                if ((theta_mm - layerTranspiration[layer]) > WSThreshold_mm):
                    isLayerStressed[layer] = False;
                    rootDensityWithoutStress += rootDensity[layer]
                else:
                    isLayerStressed[layer] = True
                    
            actualTranspiration += layerTranspiration[layer]

    # WATER STRESS [-]
    waterStress = 1. - (actualTranspiration / maxTranspiration)

    # Hydraulic redistribution
    # the movement of water from moist to dry soil through plant roots
    # TODO add numerical process
    if (waterStress > EPSILON and rootDensityWithoutStress > EPSILON):
        redistribution = min(waterStress, rootDensityWithoutStress) * maxTranspiration
        
         # redistribution acts on not stressed roots
        for layer in range(nrLayers):
            if (rootDensity[layer] > 0) and (not isLayerStressed[layer]):
                addTransp = redistribution * (rootDensity[layer] / rootDensityWithoutStress)
                layerTranspiration[layer] += addTransp;
                actualTranspiration += addTransp;
                
    # Assign transpiration flux [m3 s-1]
    for layer in range(nrLayers):
        if (layerTranspiration[layer] > 0):
            i = index + C3DStructure.nrRectangles * layer
            rate = (layerTranspiration[layer] * 0.001) / timeStep       # [m s-1]
            C3DCells[i].sinkSource = -rate * C3DCells[i].area           # [m3 s-1]
            
    return actualTranspiration



'''
def setEvaporation(index, maxEvaporation):
{
   # TODO extend to geometric soilLayers

    # surface evaporation
    double surfaceEvaporation = MINVALUE(maxEvaporation, soilLayers[0].waterContent);
    soilLayers[0].waterContent -= surfaceEvaporation;

    double actualEvaporation = surfaceEvaporation;

    double residualEvaporation = maxEvaporation - surfaceEvaporation;
    if (residualEvaporation < EPSILON)
        return actualEvaporation;

    # soil evaporation
    unsigned int lastLayerEvap = unsigned(floor(MAX_EVAPORATION_DEPTH / soilLayers[1].thickness)) +1;
    double* coeffEvap = new double[lastLayerEvap];
    double layerDepth, coeffDepth;

    double sumCoeff = 0;
    double minDepth = soilLayers[1].depth + soilLayers[1].thickness / 2;
    for (unsigned int i=1; i <= lastLayerEvap; i++)
    {
        layerDepth = soilLayers[i].depth + soilLayers[i].thickness / 2.0;

        coeffDepth = MAXVALUE((layerDepth - minDepth) / (MAX_EVAPORATION_DEPTH - minDepth), 0);
        # coeffEvap: 1 at depthMin, ~0.1 at MAX_EVAPORATION_DEPTH
        coeffEvap[i-1] = exp(-2 * coeffDepth);

        coeffEvap[i-1] = MINVALUE(1.0, exp((-layerDepth * 2.0) / MAX_EVAPORATION_DEPTH));
        sumCoeff += coeffEvap[i-1];
    }

    bool isWaterSupply = true;
    double sumEvap, evapLayerThreshold, evapLayer;
    while ((residualEvaporation > EPSILON) && (isWaterSupply == true))
    {
        isWaterSupply = false;
        sumEvap = 0.0;

        for (unsigned int i=1; i<=lastLayerEvap; i++)
        {
            evapLayerThreshold = soilLayers[i].FC - coeffEvap[i-1] * (soilLayers[i].FC - soilLayers[i].HH);
            evapLayer = (coeffEvap[i-1] / sumCoeff) * residualEvaporation;

            if (soilLayers[i].waterContent > (evapLayerThreshold + evapLayer))
                isWaterSupply = true;
            else if (soilLayers[i].waterContent > evapLayerThreshold)
                evapLayer = soilLayers[i].waterContent - evapLayerThreshold;
            else
                evapLayer = 0.0;

            soilLayers[i].waterContent -= evapLayer;
            sumEvap += evapLayer;
        }

        residualEvaporation -= sumEvap;
        actualEvaporation  += sumEvap;
    }

    delete[] coeffEvap;

    return actualEvaporation;
}
'''



