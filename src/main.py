import numpy as np
from dataStructures import *
from readDataFile import readDataFile
from fileUtilities import loadState
import soil
import waterBalance
import rectangularMesh
import criteria3D
import visual3D
import os
from PenmanMonteith import computeHourlyET0
import exportUtils
import importUtils

 
def main():
    print (os.getcwd())
    dataPath = "data"

    print("Building rectangle mesh...")
    rectangularMesh.rectangularMeshCreation()
    print ("Nr. of rectangles:", C3DStructure.nrRectangles)
    print ("Total area [m^2]:", C3DStructure.totalArea)
    
    rectangularMesh.header = rectangularMesh.getHeader(rectangularMesh.C3DRM)
   
    # SOIL
    print ("Load soil...")
    soilFolder = "soil"
    soilFile = "soil.txt"
    soilPath = os.path.join(dataPath, soilFolder, soilFile)
    soil.C3DSoil = soil.readHorizon(soilPath, 1)
    totalDepth = soil.C3DSoil.lowerDepth
    print("Soil depth [m]:", totalDepth)
    
    nrLayers, soil.depth, soil.thickness = soil.setLayers(totalDepth, 
                     C3DParameters.minThickness, C3DParameters.maxThickness, 
                     C3DParameters.geometricFactor) 
    print("Nr. of layers:", nrLayers)
    
    # Initialize memory
    criteria3D.memoryAllocation(nrLayers, C3DStructure.nrRectangles)
    print("Nr. of cells: ", C3DStructure.nrCells)
    
    print("Set cell properties...")   
    for i in range(C3DStructure.nrRectangles):
        for layer in range(C3DStructure.nrLayers): 
            [x, y, z] = rectangularMesh.C3DRM[i].centroid 
            index = i + C3DStructure.nrRectangles * layer
            elevation = z - soil.depth[layer]
            volume = float(rectangularMesh.C3DRM[i].area * soil.thickness[layer])
            criteria3D.setCellGeometry(index, x, y, 
                                elevation, volume, rectangularMesh.C3DRM[i].area)
            if (layer == 0):
                # surface 
                if rectangularMesh.C3DRM[i].isBoundary:
                    criteria3D.setCellProperties(index, True, BOUNDARY_RUNOFF)
                    criteria3D.setBoundaryProperties(index, 
                                  rectangularMesh.C3DRM[i].boundarySide, rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    criteria3D.setCellProperties(index, True, BOUNDARY_NONE)
                criteria3D.setMatricPotential(index, 0.0)
                
            elif (layer == (nrLayers-1)):
                # last layer
                if C3DParameters.isFreeDrainage:
                    criteria3D.setCellProperties(index, False, BOUNDARY_FREEDRAINAGE)
                else:
                    criteria3D.setCellProperties(index, False, BOUNDARY_NONE)
                    
                criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)
                
            else:
                if rectangularMesh.C3DRM[i].isBoundary: 
                    criteria3D.setCellProperties(index, False, BOUNDARY_FREELATERALDRAINAGE)
                    criteria3D.setBoundaryProperties(index, rectangularMesh.C3DRM[i].boundarySide * soil.thickness[layer], 
                                                     rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    criteria3D.setCellProperties(index, False, BOUNDARY_NONE)
                    
                criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)
                 
    print("Set links...")   
    for i in range(C3DStructure.nrRectangles): 
        # UP
        for layer in range(1, nrLayers):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i 
            linkIndex = index - C3DStructure.nrRectangles
            criteria3D.SetCellLink(index, linkIndex, UP, exchangeArea)   
        # LATERAL
        for neighbour in rectangularMesh.C3DRM[i].neighbours:
            if (neighbour != NOLINK):
                linkSide = rectangularMesh.getAdjacentSide(i, neighbour)
                for layer in range(nrLayers): 
                    if (layer == 0):
                        #surface: boundary length [m]
                        exchangeArea = linkSide
                    else:
                        #sub-surface: boundary area [m2]
                        exchangeArea = soil.thickness[layer] * linkSide
                    index = C3DStructure.nrRectangles * layer + i 
                    linkIndex = C3DStructure.nrRectangles * layer + neighbour
                    criteria3D.SetCellLink(index, linkIndex, LATERAL, exchangeArea)
        # DOWN
        for layer in range(nrLayers-1):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i 
            linkIndex = index + C3DStructure.nrRectangles
            criteria3D.SetCellLink(index, linkIndex, DOWN, exchangeArea)
            
    # LOAD INITIAL STATE - comment if you dont't have one
    #print ("Load initial state...")
    #loadState(dataPath + "state_0.csv")
    
    waterBalance.initializeBalance()
    print("Initial water storage [m^3]:", format(waterBalance.currentStep.waterStorage, ".3f"))
        
    print("Read arpae data...")
    arpaeFolder = "arpae"
    arpaePath = os.path.join(dataPath, arpaeFolder)
    stationInfo, arpaeData = importUtils.readArpaeData(arpaePath)

    print("Read irrigations data...")
    irrigationsFolder = "irrigations"
    irrigationsPath = os.path.join(dataPath, irrigationsFolder)
    irrigationsConfigurations, irrigationsData = importUtils.readIrrigationsData(irrigationsPath, arpaeData.iloc[0]["start"], arpaeData.iloc[-1]["start"])

    # TIME LENGHT 
    # change it if your observed data are different (ex: hourly)
    timeLength = importUtils.TIME_LENGHT * 60         # [s]
    C3DParameters.deltaT_max = timeLength
    print("Time lenght [s]:", timeLength)
    print("Total simulation time [s]:", len(arpaeData) * timeLength)
    
    visual3D.initialize(1280)
    visual3D.isPause = True
    
    # export inizialization 
    exportUtils.createExportFile()

    # main cycle
    for index, hourly_relevation in arpaeData.iterrows():
        criteria3D.cleanSurfaceSinkSource()
        waterBalance.currentPrec = hourly_relevation["precipitations"] / timeLength * 3600   #[l/hour]
        criteria3D.setRainfall(hourly_relevation["precipitations"], timeLength)
        #criteria3D.setDripIrrigation("irrigation", timeLength)
        exportUtils.takeScreenshot(hourly_relevation["start"])
        criteria3D.compute(timeLength)
    
    print ("\nEnd simulation.")   
main()
