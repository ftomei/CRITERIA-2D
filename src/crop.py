# crop.py
from dataStructures import *
import math
import soil
import rectangularMesh
import numpy as np

MAX_EVAPORATION_DEPTH = 0.2     # [m]


class CCrop:
    laiMin = NODATA             # [m2 m-2]
    laiMax = NODATA             # [m2 m-2]
    rootDepthZero = NODATA      # [m]
    rootDepthMax = NODATA       # [m]
    rootDeformation = NODATA    # [-]
    kcMax = NODATA              # [-]
    fRAW = NODATA               # [-]
    currentLAI = NODATA
    currentRootDepth = NODATA
    currentRootLength = NODATA

    def setMaxValues(self):
        self.currentLAI = self.laiMax
        self.currentRootDepth = self.rootDepthMax
        self.currentRootLength = self.currentRootDepth - self.rootDepthZero

    def setKiwifruit(self):
        self.laiMin = 1.0           # [m2 m-2]
        self.laiMax = 4.0           # [m2 m-2]
        self.rootDepthZero = 0.05   # [m]
        self.rootDepthMax = 0.92    # [m]
        self.rootDeformation = 1.51 # [-] 0: symmetric 1: cardioid 2: cardioid more accentuated
        self.kcMax = 1.59           # [-]
        self.fRAW = 0.6             # [-]
        self.setMaxValues()
        self.currentLAI = 2.89      # [m2 m-2]


kiwi = CCrop()
rootDensityKiwi = []
k_root = np.array


def initializeCrop(plantConfiguration, irrigationConfigurations):
    global rootDensityKiwi, k_root

    # initialize kiwifruit
    kiwi.setKiwifruit()

    # initialize root factor
    k_root = np.zeros(C3DStructure.nrRectangles)
    if C3DParameters.computeTranspiration:
        # line plant-sprinkler
        max_distance = plantConfiguration.iloc[0]['max_distance']
        x1 = plantConfiguration.iloc[0]['plant_x']
        y1 = plantConfiguration.iloc[0]['plant_y']
        x2 = irrigationConfigurations.iloc[0]['x']
        y2 = irrigationConfigurations.iloc[0]['y']
        a = y2 - y1
        b = x1 - x2
        c = y1 * (x2 - x1) - x1 * (y2 - y1)
        denominator = math.sqrt(a * a + b * b)

        for i in range(C3DStructure.nrRectangles):
            [x, y, z] = rectangularMesh.C3DRM[i].centroid
            line_distance = abs(a * x + b * y + c) / denominator
            # plant_cell_distance = rectangularMesh.distance2D(plant_xy, rectangularMesh.C3DRM[i].centroid)

            if line_distance > max_distance:
                k_root[i] = 0.0
            else:
                k_root[i] = 1 - (line_distance / max_distance)
    else:
        kiwi.currentLAI = 0

    # set root density
    for i in range(C3DStructure.nrRectangles):
        rootDensity = computeRootDensity(kiwi, C3DStructure.nrLayers, k_root[i])
        rootDensityKiwi.append(rootDensity)


def getCropSurfaceCover(currentLAI):
    k = 0.6  # [-] light extinction coefficient
    if (currentLAI == NODATA) or (currentLAI <= 0):
        return 0.
    else:
        return 1. - math.exp(-k * currentLAI)


def getMaxEvaporation(currentLAI, ET0):
    maxEvaporationRatio = 1.0
    cropSurfaceCover = getCropSurfaceCover(currentLAI)
    return ET0 * maxEvaporationRatio * (1. - cropSurfaceCover)


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
    cardioid = np.zeros(nrLayersWithRoot * 2)

    for i in range(nrLayersWithRoot):
        sinAlfa = 1.0 - float(i + 1.0) / float(nrLayersWithRoot)
        cosAlfa = max(math.sqrt(1.0 - math.pow(sinAlfa, 2.0)), EPSILON)
        alfa = math.atan(sinAlfa / cosAlfa)
        halfCircle[i] = ((math.pi / 2.0) - alfa - sinAlfa * cosAlfa) / math.pi

    lastLayer = nrLayersWithRoot * 2 - 1
    cardioid[0] = halfCircle[0]
    cardioid[lastLayer] = cardioid[0]
    for i in range(1, nrLayersWithRoot):
        cardioid[i] = halfCircle[i] - halfCircle[i - 1]
        cardioid[lastLayer - i] = cardioid[i]

    # cardioid deformation
    LiMin = -math.log(0.2) / float(nrLayersWithRoot)
    LiMax = -math.log(0.05) / float(nrLayersWithRoot)
    k = LiMin + (LiMax - LiMin) * (deformationFactor - 1.0)

    rootDensitySum = 0
    for i in range(nrLayersWithRoot * 2):
        cardioid[i] *= math.exp(-k * (i + 0.5))
        rootDensitySum += cardioid[i]

    # normalize
    for i in range(nrLayersWithRoot * 2):
        cardioid[i] /= rootDensitySum

    # assign layer density
    layerDensity = np.zeros(nrLayersWithRoot, np.float64)
    for i in range(nrLayersWithRoot):
        layerDensity[i] = cardioid[2 * i] + cardioid[2 * i + 1]

    return layerDensity


def computeRootDensity(crop, nrLayers, rootFactor):
    # initialize
    rootDensity = np.zeros(nrLayers, np.float64)
    if crop.currentRootLength <= 0 or rootFactor == 0:
        return rootDensity

    rootLength = crop.currentRootLength * rootFactor

    # smallest unit of computation (1 mm)
    atoms = np.zeros(nrLayers, np.int16)
    for i in range(nrLayers):
        atoms[i] = int(round(soil.thickness[i] * 1000))

    nrUnrootedAtoms = int(round(crop.rootDepthZero * 1000))
    nrRootedAtoms = int(round(rootLength * 1000))
    densityAtoms = cardioidDistribution(crop.rootDeformation, nrRootedAtoms)

    # assign root density
    counter = 0
    for layer in range(nrLayers):
        for i in range(atoms[layer]):
            if (counter >= nrUnrootedAtoms) and (counter < nrUnrootedAtoms + nrRootedAtoms):
                rootDensity[layer] += densityAtoms[counter - nrUnrootedAtoms]
            counter += 1

    # check (rootDensitySum == 1)
    rootDensitySum = 0.0
    for i in range(nrLayers):
        rootDensitySum += rootDensity[i]
    if abs(rootDensitySum - 1.0) > EPSILON:
        print("WARNING! Sum of root density:", rootDensitySum)

    return rootDensity


# assign hourly transpiration
def setTranspiration(surfaceIndex, crop, rootDensity, maxTranspiration):
    if maxTranspiration < EPSILON:
        return 0.0

    FC = soil.getFieldCapacityWC()  # [m3 m-3] water content at field capacity
    WP = soil.getWiltingPointWC()  # [m3 m-3] water content at wilting point
    WSThreshold = FC - crop.fRAW * (FC - WP)  # [m3 m-3] water scarcity stress threshold

    # Initialize
    rootDensityWithoutStress = 0.0  # [-]
    actualTranspiration = 0.0  # [mm]

    nrLayers = len(rootDensity)
    isLayerStressed = np.zeros(nrLayers, dtype=bool)
    layerTranspiration = np.zeros(nrLayers, np.float64)
    for layer in range(nrLayers):
        isLayerStressed[layer] = False
        layerTranspiration[layer] = 0

    for layer in range(nrLayers):
        if rootDensity[layer] > 0:
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            theta = soil.getVolumetricWaterContent(i)
            # TODO water surplus

            # WATER SCARCITY
            if theta < WSThreshold:
                if theta <= WP:
                    layerTranspiration[layer] = 0
                else:
                    layerTranspiration[layer] = maxTranspiration * rootDensity[layer] * (
                                (theta - WP) / (WSThreshold - WP))
                isLayerStressed[layer] = True
            else:
                # normal conditions
                layerTranspiration[layer] = maxTranspiration * rootDensity[layer]

                # check stress
                theta_mm = theta * soil.thickness[layer] * 1000.
                WSThreshold_mm = WSThreshold * soil.thickness[layer] * 1000.

                if (theta_mm - layerTranspiration[layer]) > WSThreshold_mm:
                    isLayerStressed[layer] = False
                    rootDensityWithoutStress += rootDensity[layer]
                else:
                    isLayerStressed[layer] = True

            actualTranspiration += layerTranspiration[layer]

    # WATER STRESS [-]
    waterStress = 1. - (actualTranspiration / maxTranspiration)

    # Hydraulic redistribution
    # the movement of water from moist to dry soil through plant roots
    # TODO add numerical process
    if waterStress > EPSILON and rootDensityWithoutStress > EPSILON:
        redistribution = min(waterStress, rootDensityWithoutStress) * maxTranspiration

        # redistribution acts on not stressed roots
        for layer in range(nrLayers):
            if (rootDensity[layer] > 0) and (not isLayerStressed[layer]):
                addTranspiration = redistribution * (rootDensity[layer] / rootDensityWithoutStress)
                layerTranspiration[layer] += addTranspiration
                actualTranspiration += addTranspiration

    # Assign transpiration flux [m3 s-1]
    for layer in range(nrLayers):
        if layerTranspiration[layer] > 0:
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            rate = (layerTranspiration[layer] * 0.001) / 3600.0  # [m s-1]
            C3DCells[i].sinkSource -= rate * C3DCells[i].area  # [m3 s-1]

    return actualTranspiration


def setEvaporation(surfaceIndex, maxEvaporation):
    # TODO: enable surface evaporation - numerical problem
    """
    surfaceWater = (C3DCells[surfaceIndex].H - C3DCells[surfaceIndex].z)    # [m]
    surfaceEvaporation = min(maxEvaporation, surfaceWater * 1000.0)         # [mm]
    rate = (surfaceEvaporation * 0.001) / 3600.0                            # [m s-1]
    C3DCells[surfaceIndex].sinkSource -= rate * C3DCells[surfaceIndex].area     # [m3 s-1]
    """
    surfaceEvaporation = 0
    actualEvaporation = surfaceEvaporation
    residualEvaporation = maxEvaporation - surfaceEvaporation

    if residualEvaporation < EPSILON:
        return actualEvaporation

    # soil evaporation
    FC = soil.getFieldCapacityWC()  # [m3 m-3] water content at field capacity
    HH = soil.getHygroscopicWC()    # [m3 m-3] water content at Hygroscopic moisture
    half_FC = HH + (FC - HH) * 0.5
    lastIndex = 0
    while (lastIndex < len(soil.depth)) and (soil.depth[lastIndex] <= MAX_EVAPORATION_DEPTH):
        lastIndex += 1

    nrEvapLayers = lastIndex

    # evaporation coefficient: 1 at first layer, ~0.1 at MAX_EVAPORATION_DEPTH
    coeffEvap = np.zeros(nrEvapLayers, np.float64)
    sumCoefficient = 0
    for i in range(1, nrEvapLayers):
        coeffDepth = max((soil.depth[i] - soil.depth[1]) / (MAX_EVAPORATION_DEPTH - soil.depth[1]), 0)
        coeffEvap[i] = math.exp(-coeffDepth * math.e)
        sumCoefficient += (coeffEvap[i] * soil.thickness[i])

    isWaterSupply = True
    while (residualEvaporation > EPSILON) and (isWaterSupply == True):
        isWaterSupply = False
        sumEvaporation = 0.0

        for layer in range(1, nrEvapLayers):
            index = surfaceIndex + C3DStructure.nrRectangles * layer
            theta = soil.getVolumetricWaterContent(index)                       # [m3 m-3]
            evaporationThreshold = half_FC - coeffEvap[layer] * (half_FC - HH)  # [m3 m-3]
            evaporation = residualEvaporation * ((coeffEvap[layer]
                                                  * soil.thickness[layer]) / sumCoefficient)  # [mm]
            if theta < evaporationThreshold:
                evaporation = 0.0
            else:
                availableWC = (theta - evaporationThreshold) * soil.thickness[layer] * 1000.  # [mm]
                if availableWC <= evaporation:
                    evaporation = availableWC
                else:
                    isWaterSupply = True

            sumEvaporation += evaporation

            rate = (evaporation * 0.001) / 3600.  # [m s-1]
            C3DCells[index].sinkSource -= rate * C3DCells[index].area  # [m3 s-1]

        residualEvaporation -= sumEvaporation
        actualEvaporation += sumEvaporation

    return actualEvaporation


def setEvapotranspiration(ET0):

    # initialize sinkSource
    for i in range(C3DStructure.nrCells):
        C3DCells[i].sinkSource = 0

    if C3DParameters.computeTranspiration:
        for i in range(C3DStructure.nrRectangles):
            maxTrKiwi = getMaxTranspiration(kiwi.currentLAI, kiwi.kcMax, ET0)
            maxTranspiration = k_root[i] * maxTrKiwi
            setTranspiration(i, kiwi, rootDensityKiwi[i], maxTranspiration)

    if C3DParameters.computeEvaporation:
        maxEvaporation = getMaxEvaporation(kiwi.currentLAI, ET0)
        for i in range(C3DStructure.nrRectangles):
            setEvaporation(i, maxEvaporation)
