# dataStructures.py

from commonConst import *


class C3DStructure:
    nrDimensions = 3
    nrVerticesPerRectangle = 4
    z0 = 0.0            # [m] z
    slopePlant = 0.2    # [-] slope around plant (baulatura)
    slopeX = -0.01      # [-] slope on x axis
    slopeY = -0.025     # [-] slope on y axis
    gridWidth = 2.375   # [m] x axis
    gridHeight = 0.625  # [m] y axis
    gridStep = 0.125    # [m]
    slopeSide = 1.0     # [m] side of baulatura

    nrRectanglesInXAxis = int(gridWidth / gridStep)
    if abs(gridWidth - gridStep * nrRectanglesInXAxis) > EPSILON:
        nrRectanglesInXAxis += 1

    nrRectanglesInYAxis = int(gridHeight / gridStep)
    if abs(gridHeight - gridStep * nrRectanglesInYAxis) > EPSILON:
        nrRectanglesInYAxis += 1

    nrRectangles = nrRectanglesInXAxis * nrRectanglesInYAxis

    nrLayers = 0
    nrCells = 0
    nrLateralLinks = 4
    nrMaxLinks = 6
    totalArea = 0


class Clink:
    def __init__(self):
        self.index = NOLINK     # [-] index of linked cell
        self.area = NODATA      # [m^2] area of interface
        self.distance = NODATA  # [m]


class Cboundary:
    def __init__(self):
        self.type = BOUNDARY_NONE
        self.area = NODATA      # [m^2] area of interface
        self.slope = NODATA     # [-] slope
        self.flow = NODATA      # [m^3 s^-1] boundary water flow


class Ccell:
    def __init__(self):
        self.x = NODATA         # [m]
        self.y = NODATA         # [m]
        self.z = NODATA         # [m]
        self.volume = NODATA    # [m^3] volume
        self.area = NODATA      # [m^3] area (for surface cells)
        self.isSurface = False  # true if the cell is on surface
        self.Se = NODATA        # [-] degree of saturation
        self.H = NODATA         # [m] current total water potential
        self.H0 = NODATA        # [m] total water potential
        self.k = NODATA         # [m s^-1] hydraulic conductivity
        self.sinkSource = NODATA  # [m^3 s^-1] water sink/source
        self.flow = NODATA      # [m^3 s^-1] sink/source + boundary
        self.boundary = Cboundary()
        self.upLink = Clink()
        self.downLink = Clink()
        self.lateralLink = [Clink(), Clink(), Clink(), Clink()]


# model parameters
class C3DParameters:
    initialWaterPotential = -4.0        # [m]

    # water retention curve and conductivity
    waterRetentionCurve = IPPISCH_VG
    conductivityMean = LOGARITHMIC
    conductivityHVRatio = 1.0

    # infiltration
    computeInfiltration = True

    # soil layers thickness
    minThickness = 0.01                 # [m]
    maxThickness = 0.04                 # [m]
    maxThicknessAt = 0.2                # [m]

    # sink-source
    assignIrrigation = True
    computeEvaporation = True
    computeTranspiration = True

    # surface flow
    computeSurfaceFlow = True
    roughness = 0.24                    # [s m^0.33]
    pondIrrigation = 0.01               # [m]
    pondRainfall = 0.002                # [m]
    pond = pondRainfall                 # [m]

    # boundary
    isSurfaceRunoff = True
    isFreeLateralDrainage = True
    isFreeDrainage = True
    isWaterTable = False
    waterTableDepth = -6.0              # [m]

    # numerical solution parameters
    currentDeltaT = 60.0                # [s]
    deltaT_min = 2                      # [s]
    deltaT_max = 600.0                  # [s]
    maxIterationsNr = 100
    maxApproximationsNr = 10
    residualTolerance = 1E-12
    MBRThreshold = 1E-5


# global
C3DCells = []
