#dataStructures.py

from commonConst import  *

class C3DStructure:
    nrLayers = 0
    nrTriangles = 0
    nrCells = 0
    nrLateralLinks = 3
    nrMaxLinks = 5
    totalArea = 0

class Clink:
    def __init__(self):
        self.index = NOLINK               # [-] index of linked cell
        self.area = NODATA                # [m^2] area of interface
        self.distance = NODATA            # [m]
    
class Cboundary:
    def __init__(self):
        self.type = BOUNDARY_NONE
        self.area = NODATA                # [m^2] area of interface
        self.slope = NODATA               # [-] slope
        self.flow = NODATA                # [m^3 s^-1] boundary water flow 
    
class Ccell:
    def __init__(self):
        self.x = NODATA                   # [m]
        self.y = NODATA                   # [m]            
        self.z = NODATA                   # [m]
        self.volume = NODATA              # [m^3] volume
        self.isSurface = False            # true if the cell is on surface
        self.Se = NODATA                  # [-] degree of saturation
        self.H = NODATA                   # [m] total water potential
        self.H0 = NODATA                  # [m] previous total water potential
        self.k = NODATA                   # [m s^-1] hydraulic conductivity
        self.sinkSource = NODATA          # [m^3 s^-1] water sink/source
        self.flow = NODATA                # [m^3 s^-1] sink/source + boundary
        self.boundary = Cboundary() 
        self.upLink = Clink()
        self.downLink = Clink()
        self.lateralLink = [Clink(), Clink(), Clink()]
        
# user choices
class C3DParameters:
    waterRetentionCurve = IPPISCH_VG
    meanType = LOGARITHMIC
    initialWaterPotential = -2.0                # [m]
    computeOnlySurface = False
    isFreeDrainage = True
    minThickness = 0.01                         # [m]
    maxThickness = 0.1                          # [m]
    geometricFactor = 1.2
    roughness = 0.24                            # [s m^0.33]
    pond = 0.002                                # [m]
    currentDeltaT = 60.0                        # [s]
    deltaT_min = 6.0                            # [s]
    deltaT_max = 900.0                          # [s]
    maxIterationsNr = 100
    maxApproximationsNr = 10
    residualTolerance = 1E-12
    MBRThreshold = 1E-3
    conductivityHVRatio = 10.0

#global
C3DCells = []
