#crop.py
from dataStructures import  *
import math
import soil
import rectangularMesh
import numpy as np

MAX_EVAPORATION_DEPTH = 0.15                # [m]

class Ccrop:
    laiMin = NODATA                         # [m2 m-2]
    laiMax = NODATA                         # [m2 m-2]
    rootDepthZero = NODATA                  # [m]
    rootDepthMax = NODATA                   # [m]
    rootDeformation = NODATA                # [-]
    kcMax = NODATA                          # [-]
    fRAW = NODATA                           # [-]
    currentLAI = NODATA              
    currentRootDepth = NODATA
    currentRootLenght = NODATA
    
    def setMaxValues(self):
        self.currentLAI = self.laiMax     
        self.currentRootDepth = self.rootDepthMax
        self.currentRootLenght = self.currentRootDepth - self.rootDepthZero
        
    def setKiwifruit(self):
        self.laiMin = 1.0                        # [m2 m-2]
        self.laiMax = 3.5                        # [m2 m-2]
        self.rootDepthZero = 0.1                 # [m]
        self.rootDepthMax = 1.0                  # [m]
        self.rootDeformation = 1.0               # [-]
        self.kcMax = 1.1                         # [-]
        self.fRAW = 0.55                         # [-]
        self.setMaxValues()

    def setGrass(self):
        self.laiMin = 0.5                       # [m2 m-2]
        self.laiMax = 2.0                       # [m2 m-2]
        self.rootDepthZero = 0.05               # [m]
        self.rootDepthMax = 0.8                 # [m]
        self.rootDeformation = 1.0              # [-]
        self.kcMax = 1.0                        # [-]
        self.fRAW = 0.65                        # [-]
        self.setMaxValues()


kiwi = Ccrop()
grass = Ccrop()


def initializeCrop():
    global rootDensityGrass, rootDensityKiwi, LAI_kiwi, LAI_grass, surfaceEvaporation
    # kiwifruit
    kiwi.setKiwifruit()
    kiwi.currentLAI = 2.0
    rootDensityKiwi = computeRootDensity(kiwi, C3DStructure.nrLayers)
    # grass
    grass.setGrass()
    grass.currentLAI = 1.0
    rootDensityGrass = computeRootDensity(grass, C3DStructure.nrLayers)
    # assign LAI
    LAI_kiwi = np.zeros(C3DStructure.nrRectangles)
    LAI_grass = np.zeros(C3DStructure.nrRectangles)
    for i in range(C3DStructure.nrRectangles):
        [x, y, z] = rectangularMesh.C3DRM[i].centroid
        # assign kiwi to whole  area
        LAI_kiwi[i] = kiwi.currentLAI
        # assign grass to right area
        if (x >= 0.8):
            LAI_grass[i] = grass.currentLAI
    # initialize surface evaporation
    surfaceEvaporation = np.zeros(C3DStructure.nrRectangles)
            
    
def getCropSurfaceCover(currentLAI):
    k = 0.6      # [-] light extinction coefficient
    if (currentLAI == NODATA) or (currentLAI <= 0):
        return 0.
    else:
        return 1. - math.exp(-k * currentLAI)


def getMaxEvaporation(currentLAI, ET0):
    maxEvapRatio = 0.66
    cropSurfaceCover = getCropSurfaceCover(currentLAI)
    return ET0 * maxEvapRatio * (1. - cropSurfaceCover)


def getMaxTranspiration(currentLAI, kcMax, ET0):
    if (currentLAI == NODATA) or (currentLAI <= 0):
        return 0.
    else:
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


# assign hourly transpiration
def setTranspiration(surfaceIndex, crop, rootDensity, maxTranspiration):
    if (maxTranspiration < EPSILON):
        return 0.0 
    
    FC = soil.getFieldCapacityWC()                  # [m3 m-3] water content at field capacity
    WP = soil.getWiltingPointWC()                   # [m3 m-3] water content at wilting point
    WSThreshold = FC - crop.fRAW * (FC - WP);       # [m3 m-3] water scarcity stress threshold
    
    # Initialize
    rootDensityWithoutStress = 0.0                  # [-]
    actualTranspiration = 0.0                       # [mm]
    
    nrLayers = len(rootDensity)
    isLayerStressed = np.zeros(nrLayers, dtype=bool)
    layerTranspiration = np.zeros(nrLayers, np.float64)
    for layer in range(nrLayers):
        isLayerStressed[layer] = False
        layerTranspiration[layer] = 0
    
    for layer in range(nrLayers):
        if (rootDensity[layer] > 0):
            i = surfaceIndex + C3DStructure.nrRectangles * layer
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
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            rate = (layerTranspiration[layer] * 0.001) / 3600.0         # [m s-1]
            if (C3DCells[i].sinkSource == NODATA):
                C3DCells[i].sinkSource = -rate * C3DCells[i].area       # [m3 s-1]
            else:
                C3DCells[i].sinkSource -= rate * C3DCells[i].area       # [m3 s-1]
            
    return actualTranspiration


def setEvaporation(surfaceIndex, maxEvaporation):
    # surface evaporation
    surfaceWater = C3DCells[surfaceIndex].H - C3DCells[surfaceIndex].z
    surfaceEvaporation = min(maxEvaporation, surfaceWater)
    rate = (surfaceEvaporation * 0.001) / 3600.0                                # [m s-1]
    C3DCells[surfaceIndex].sinkSource -= rate * C3DCells[surfaceIndex].area     # [m3 s-1]

    actualEvaporation = surfaceEvaporation
    residualEvaporation = maxEvaporation - surfaceEvaporation
    
    if (residualEvaporation < EPSILON):
        return actualEvaporation
    
    # soil evaporation
    '''
    lastLayerEvap = unsigned(floor(MAX_EVAPORATION_DEPTH / soilLayers[1].thickness)) +1
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
    '''
    return actualEvaporation


def setEvapotranspiration(ET0):
    global surfaceEvaporation
    
    # clean sinkSource and surface evaporation
    for i in range(C3DStructure.nrCells):
        C3DCells[i].sinkSource = 0 
    for i in range(C3DStructure.nrRectangles):
        surfaceEvaporation[i] = 0
        
    if C3DParameters.computeTranspiration:
        for i in range(C3DStructure.nrRectangles):
            # kiwifruit
            maxTranspKiwi = getMaxTranspiration(LAI_kiwi[i], kiwi.kcMax, ET0)
            setTranspiration(i, kiwi, rootDensityKiwi, maxTranspKiwi) 
            # grass
            maxTranspGrass = getMaxTranspiration(LAI_grass[i], grass.kcMax, ET0)
            setTranspiration(i, grass, rootDensityGrass, maxTranspGrass) 
                
    if C3DParameters.computeEvaporation:
        for i in range(C3DStructure.nrRectangles):
            LAIsum = LAI_grass[i] + LAI_kiwi[i]
            maxEvaporation = getMaxEvaporation(LAIsum, ET0)
            setEvaporation(i, maxEvaporation) 
            

