#visual3D.py

import vpython as visual
from dataStructures import *
from color import *
import fileUtilities
import waterBalance
import rectangularMesh
import soil
import exportUtils

from copy import copy


# surfaceRectangles = []
sliceRectangles = []
subSurfaceRectangles = []
visualizedLayer = 0
visualizedSlice = int(C3DStructure.nrRectanglesInYAxis / 2)
nrColorLevels = 10
degreeMaximum = 1.0
degreeMinimum = 0.4
waterLevelMaximum = max(0.005, C3DParameters.pond)
isWaterPotential = False
isRootVisualization = False
isPause = False
  

def initialize(totalWidth):
    global sliceCanvas, soilCanvas, interface
    global sliceLabel, layerLabel, timeLabel, precLabel, irrLabel, stepLabel, storageLabel
    global flowLabel, totalFlowLabel, totalErrorLabel, colorScale
    global visualizedSlice, visualizedLayer
    
    setAllColorScale()
    
    #CANVAS DIMENSION
    interfaceWidth = int(totalWidth * 0.2)
    dx = int((totalWidth-interfaceWidth) / 2.0)
    dy = int(dx * 0.8)
    h  = int(dy / 30)
    
    #CENTER
    cX = (rectangularMesh.header.xMin + rectangularMesh.header.xMax) * 0.5
    cY = (rectangularMesh.header.yMin + rectangularMesh.header.yMax) * 0.5
    cZ = rectangularMesh.header.zMin * rectangularMesh.header.magnify
    #RANGE
    rx = rectangularMesh.header.xMax - rectangularMesh.header.xMin 
    ry = rectangularMesh.header.yMax - rectangularMesh.header.yMin
    rangeXY = max(rx, ry) * rectangularMesh.header.magnify
    lastLayer = C3DStructure.nrLayers-1
    rz = soil.depth[lastLayer] + soil.thickness[lastLayer]*0.5
    rangeZ = rz * rectangularMesh.header.magnify

    #INTERFACE CANVAS
    interface = visual.canvas(width = interfaceWidth, height = dy, align="left")
    interface.background = visual.color.white
    interface.center = visual.vector(1,0,0)
    interface.range = 5

    timeLabel = visual.label(height=h, pos=visual.vector(0,3,0), text = "", canvas = interface)
    stepLabel = visual.label(height=h, pos=visual.vector(0,2,0), text = "", canvas = interface)
    precLabel = visual.label(height=h,pos=visual.vector(0,1,0), text = "", canvas = interface)
    irrLabel = visual.label(height=h,pos=visual.vector(0,0,0), text = "", canvas = interface)
    storageLabel = visual.label(height=h,pos=visual.vector(0,-1,0), text = "", canvas = interface)
    flowLabel = visual.label(height=h, pos=visual.vector(0,-2,0), text = "", canvas = interface)
    totalFlowLabel = visual.label(height=h, pos=visual.vector(0,-3,0), text = "", canvas = interface)
    totalErrorLabel = visual.label(height=h, pos=visual.vector(0,-4,0), text = "", canvas = interface)
    
    #COLOR LEGEND
    stepY = 10 / nrColorLevels
    colorScale = []
    for i in range (nrColorLevels+1):
        label = visual.label(canvas = interface, pos=visual.vector(4, -5+(i*stepY), 0), height=h, 
                             background = visual.vec(0,0,0))
        colorScale.append(label)

    #SURFACE CANVAS
    soilCanvas = visual.canvas(width = dx, height = dy, align="left")
    soilCanvas.background = visual.color.white
    soilCanvas.center = visual.vector(cX, cY, cZ+(rangeXY*0.2))
    soilCanvas.ambient = visual.vector(0.5, 0.5, 0.5)
    soilCanvas.up = visual.vector(0,0,1)
    soilCanvas.forward = visual.vector(0, 0.01,-1.0)
    soilCanvas.range = rangeXY
    layerLabel = visual.label(canvas = soilCanvas, height = h, pos=visual.vector(cX, cY+rangeXY*0.8, cZ))
    
    drawColorScale()
    drawSubSurface(True)

    #SLICE CANVAS
    sliceCanvas = visual.canvas(width = dx, height = dy, align="left")
    sliceCanvas.background = visual.color.white
    sliceCanvas.center = visual.vector(cX, cY, cZ-(rangeZ*0.5))
    sliceCanvas.ambient = visual.vector(0.5, 0.5, 0.5)
    sliceCanvas.up = visual.vector(0, 0, 1)
    sliceCanvas.forward = visual.vector(0, 1, 0)
    sliceCanvas.range = rangeZ
    sliceLabel = visual.label(canvas = sliceCanvas, height = h, pos=visual.vector(cX, cY, cZ+(rangeZ*0.3)))
    
    sliceCanvas.caption = " *** COMMANDS ***\n\n 'r': run simulation \n 'p': pause "
    sliceCanvas.caption += "\n '^': move up (soil layer) \n 'v': move down (soil layer) "
    sliceCanvas.caption += "\n '<': move left (soil slice) \n '>': move right (soil slice) "
    sliceCanvas.caption += "\n 's': save state \n 'l': load state "
    sliceCanvas.caption += "\n 'c': colorscale range"
    sliceCanvas.caption += "\n 'w': switch water content/potential"
    sliceCanvas.caption += "\n 'q': switch to root factor visualization"
    
    drawSlice(True)
    updateInterface()
    interface.bind('keydown', keyInput)


def drawColorScale():
    if (isRootVisualization):
        for i in range (nrColorLevels):
            colorScale[i].visible = True 
        step = 1 / nrColorLevels
        for i in range (nrColorLevels+1):
            value = step * i
            c = getSEColor(value, 0, 1)
            colorScale[i].background = visual.vector(c[0], c[1], c[2])
            colorScale[i].text = format(value,".2f")
    elif (isWaterPotential):
        i = nrColorLevels
        c = getMatricPotentialColor(-1)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-10kPa"
        i -= 1
        c = getMatricPotentialColor(-2)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-20    "
        i -= 1
        c = getMatricPotentialColor(-3)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-30    "
        i -= 1
        c = getMatricPotentialColor(-5)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-50    "
        i -= 1
        c = getMatricPotentialColor(-10)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-100  "
        i -= 1
        c = getMatricPotentialColor(-30)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-300  "
        i -= 1
        c = getMatricPotentialColor(-100)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-1000"
        for j in range (i):
            colorScale[j].visible = False
    else:  
        for i in range (nrColorLevels):
            colorScale[i].visible = True 
        step = (degreeMaximum - degreeMinimum) / nrColorLevels
        for i in range (nrColorLevels+1):
            degree = degreeMinimum + step * i
            c = getSEColor(degree, degreeMinimum, degreeMaximum)
            colorScale[i].background = visual.vector(c[0], c[1], c[2])
            colorScale[i].text = format(degree,".2f")


def updateColorScale():
    global degreeMinimum, degreeMaximum
    
    rangeOk = False
    while not rangeOk:
        degreeMinimum = float(input("\nSet min. value (Degree of saturation):"))
        degreeMaximum = float(input("Set max. value (Degree of saturation):"))
        if degreeMaximum <= degreeMinimum:
            print("Wrong range!")
        else:
            rangeOk = True
            
    drawColorScale()
    drawSubSurface(False)
    drawSlice(False)


def updateLayer(s):
    global visualizedLayer
    
    if s == 'down':
        if visualizedLayer < C3DStructure.nrLayers-1:
            visualizedLayer += 1
    elif s == 'up':
        if visualizedLayer > 0:
            visualizedLayer -= 1
             
    updateInterface()
    drawSubSurface(False)


def updateSlice(s):
    global visualizedSlice
    
    if s == 'right':
        if visualizedSlice < C3DStructure.nrRectanglesInYAxis-1:
            visualizedSlice += 1
    elif s == 'left':
        if visualizedSlice > 0:
            visualizedSlice -= 1
             
    updateInterface()
    drawSlice(False)
          
                
def keyInput(evt):
    global isPause, isWaterPotential, isRootVisualization
    s = evt.key
    if s == 'r':
        isPause = False
    elif s == 'p':
        isPause = True
        print ("Pause...")
    elif s == 'up' or s == 'down':
        updateLayer(s)
    elif s == 'left' or s == 'right':
        updateSlice(s)
    elif s == 's':
        isPause = True
        print ("Save State...")
        fileUtilities.saveState()
    elif s == 'l':
        isPause = True
        print ("Load State...")
        if fileUtilities.loadState(""):
            waterBalance.initializeBalance()
            redraw()
    elif s == "c":
        if (isPause):
            updateColorScale()
    elif s == "w":
        isWaterPotential = not isWaterPotential
        drawColorScale()
        redraw()
    elif s == "q":
        isRootVisualization = not isRootVisualization
        drawColorScale()
        redraw()


def getNewRectangle(myColor, myCanvas, v):
    vert = []
    for i in range(C3DStructure.nrVerticesPerRectangle):
        vert.append(visual.vertex(pos = visual.vector(v[i,0], v[i,1], v[i,2] * rectangularMesh.header.magnify)))
        vert[i].color = myColor
                     
    newRectangle = visual.quad(canvas = myCanvas, vs=[vert[0],vert[1],vert[3],vert[2]])
    return newRectangle 


def drawSlice(isFirst):
    global sliceRectangles
    from crop import k_root, rootDensity

    firstIndex = visualizedSlice * C3DStructure.nrRectanglesInXAxis
    posY = C3DCells[firstIndex].y

    if isRootVisualization:
        var = "Root density"
    elif isWaterPotential:
        var = "Water potential"
    else:
        var = "Degree of saturation"
    sliceLabel.text = var + " - slice at " + format(posY*100, ".1f") + "cm"

    for z in range(C3DStructure.nrLayers):
        for x in range(C3DStructure.nrRectanglesInXAxis):
            index = firstIndex + x + (z * C3DStructure.nrRectangles)
            i = z * C3DStructure.nrRectanglesInXAxis + x
            if isRootVisualization:
                surfaceIndex = firstIndex + x
                c = getSEColor(2 * k_root[surfaceIndex]*rootDensity[surfaceIndex][z], 0, 1)
            elif isWaterPotential:
                c = getMatricPotentialColor(C3DCells[index].H - C3DCells[index].z)
            else:
                c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
            myColor = visual.vector(c[0], c[1], c[2])
            
            if isFirst:
                vertices = copy(rectangularMesh.C3DRM[firstIndex + x].v)
                for v in vertices[:2]:
                    v[2] = v[2] - soil.depth[z] + (soil.thickness[z] * 0.5)
                vertices[2] = vertices[0]
                vertices[3] = vertices[1]
                for v in vertices[2:]:
                    v[2] = v[2] - soil.thickness[z]
                newRectangle = getNewRectangle(myColor, sliceCanvas, vertices)
                sliceRectangles.append(newRectangle)
            else:
                sliceRectangles[i].v0.color = myColor
                sliceRectangles[i].v1.color = myColor
                sliceRectangles[i].v2.color = myColor  
                sliceRectangles[i].v3.color = myColor
 
 
def drawSubSurface(isFirst):
    global subSurfaceRectangles

    from crop import k_root
    maxWaterLevel = 0
    for i in range(C3DStructure.nrRectangles):
        index = visualizedLayer * C3DStructure.nrRectangles + i
        # color
        if visualizedLayer == 0:
            if isRootVisualization:
                c = getSEColor(k_root[i], 0, 1)
            else:
                waterLevel = max(C3DCells[i].H - C3DCells[i].z, 0.0)
                maxWaterLevel = max(waterLevel, maxWaterLevel)
                c = getSurfaceWaterColor(waterLevel, waterLevelMaximum)
        else:
            if isWaterPotential:
                c = getMatricPotentialColor(C3DCells[index].H - C3DCells[index].z)
            else:
                c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
            
        myColor = visual.vector(c[0], c[1], c[2])
        
        if isFirst:
            newRectangle = getNewRectangle(myColor, soilCanvas, rectangularMesh.C3DRM[i].v)
            subSurfaceRectangles.append(newRectangle)
        else:
            subSurfaceRectangles[i].v0.color = myColor
            subSurfaceRectangles[i].v1.color = myColor
            subSurfaceRectangles[i].v2.color = myColor  
            subSurfaceRectangles[i].v3.color = myColor 
    # label
    if visualizedLayer == 0:
        if isRootVisualization:
            layerLabel.text = "Root factor"
        else:
            layerLabel.text = "Surface water level - max:" + format(maxWaterLevel * 1000,".1f")+"mm"
    else:
        depth = soil.depth[visualizedLayer] * 100
        if isWaterPotential:
            var = "Water potential"
        else:
            var = "Degree of saturation"
        layerLabel.text = var + " - layer at " + format(depth, ".1f")+"cm"
    
    
def updateInterface():       
    timeLabel.text = "Time: " + format(waterBalance.totalTime / 3600.0, ".3f") + " [h]"
    precLabel.text = "Rainfall: " + format(waterBalance.currentPrec,".1f") + " [mm/hour]"
    irrLabel.text = "Irrigation: " + format(waterBalance.currentIrr,".3f") + " [l/hour]"
    storage = waterBalance.currentStep.waterStorage
    flow = waterBalance.currentStep.waterFlow
    timeStep = C3DParameters.currentDeltaT
    totalFlow = waterBalance.allSimulation.waterFlow
    totalError = waterBalance.allSimulation.MBE
    stepLabel.text = "Time step: " + str(timeStep) +" [s]"
    storageLabel.text = "Storage: " + format(storage,".4f") +" [m3]"
    flowLabel.text = "Flow: " + format(flow * 1000,".4f") + " [l]"
    totalFlowLabel.text = "Total flow: " + format(totalFlow * 1000,".4f") + " [l]"
    totalErrorLabel.text = "Total error: " + format(totalError * 1000,".4f") + " [l]"
    
    
def redraw():
    updateInterface()
    drawSlice(False)
    drawSubSurface(False)
