# dataStructures.py

from commonConst import *


class C3DStructure:
    nrDimensions = 3
    nrVerticesPerRectangle = 4
    nrLateralLinks = 4
    nrMaxLinks = 6

    z0 = NODATA                 # [m] topographic elevation (center grid)
    gridStep = NODATA           # [m] step on x and y axis
    gridWidth = NODATA          # [m] grid width (axis x)
    gridHeight = NODATA         # [m] grid height (axis y)

    slopeX = NODATA             # [-] slope on x axis
    slopeY = NODATA             # [-] slope on y axis

    plantSlope = NODATA         # [-] slope around plant (baulatura)
    plantSlopeWidth = NODATA    # [m] width of baulatura

    nrLayers = NODATA
    nrCells = NODATA
    totalArea = NODATA

    nrRectanglesInXAxis = NODATA
    nrRectanglesInYAxis = NODATA
    nrRectangles = NODATA


def initialize3DStructure(width, height, gridStep):
    C3DStructure.gridStep = gridStep

    C3DStructure.nrRectanglesInXAxis = int(width / gridStep)
    if (C3DStructure.nrRectanglesInXAxis % 2) == 0:
        C3DStructure.nrRectanglesInXAxis += 1

    C3DStructure.nrRectanglesInYAxis = int(height / gridStep)
    if (C3DStructure.nrRectanglesInYAxis % 2) == 0:
        C3DStructure.nrRectanglesInYAxis += 1

    C3DStructure.nrRectangles = C3DStructure.nrRectanglesInXAxis * C3DStructure.nrRectanglesInYAxis

    C3DStructure.gridWidth = C3DStructure.nrRectanglesInXAxis * gridStep
    C3DStructure.gridHeight = C3DStructure.nrRectanglesInYAxis * gridStep


class C3DStructure_old:
    nrDimensions = 3
    nrVerticesPerRectangle = 4

    z0 = 100.0          # [m] topographic elevation
    slopePlant = 0.2    # [-] slope around plant (baulatura)
    slopeX = -0.01      # [-] slope on x axis
    slopeY = -0.025     # [-] slope on y axis
    slopeSide = 1.0     # [m] side of baulatura

    gridStep = 0.125    # [m]
    DX = 2.5            # [m]
    DY = 0.66           # [m]

    nrRectanglesInXAxis = int(DX / gridStep)
    if (nrRectanglesInXAxis % 2) == 0:
        nrRectanglesInXAxis += 1

    nrRectanglesInYAxis = int(DY / gridStep)
    if (nrRectanglesInYAxis % 2) == 0:
        nrRectanglesInYAxis += 1

    gridWidth = nrRectanglesInXAxis * gridStep      # [m] x axis
    gridHeight = nrRectanglesInYAxis * gridStep     # [m] x axis

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
    computeSurfaceFlow = False
    roughness = 0.24                    # [s m^0.33]
    pondIrrigation = 0.005              # [m]
    pondRainfall = 0.002                # [m]
    pond = pondRainfall                 # [m]

    # boundary
    isSurfaceRunoff = False
    isFreeLateralDrainage = True
    isFreeDrainage = True
    isWaterTable = False
    waterTableDepth = -6.0              # [m]

    # numerical solution parameters
    currentDeltaT = 60.0                # [s]
    deltaT_min = 1                      # [s]
    deltaT_max = 600.0                  # [s]
    maxIterationsNr = 100
    maxApproximationsNr = 10
    residualTolerance = 1E-12
    MBRThreshold = 1E-5


# global
C3DCells = []
dripperIndices = []
plantIndices = []
