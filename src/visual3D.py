#visual3D.py

import vpython as visual
from dataStructures import *
from color import *
import waterBalance
import rectangularMesh
import soil
import exportUtils

from copy import copy


#surfaceRectangles = []
sliceRectangles = []
subSurfaceRectangles = []
visualizedLayer = 0
visualizedSlice = C3DStructure.nrRectanglesInXAxis * int(C3DStructure.nrRectanglesInYAxis / 2)
nrColorLevels = 10
degreeMaximum = 1
degreeMinimum = 0.5
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
    h = int(dy / 30)
    
    #CENTER
    cX = (rectangularMesh.header.xMin + rectangularMesh.header.xMax) * 0.5
    cY = (rectangularMesh.header.yMin + rectangularMesh.header.yMax) * 0.5
    cZ = rectangularMesh.header.zMin * rectangularMesh.header.magnify
    Zlabel = (rectangularMesh.header.zMax + rectangularMesh.header.dz*0.75 + 0.2) * rectangularMesh.header.magnify

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
    soilCanvas.center = visual.vector(cX, cY, cZ)
    soilCanvas.ambient = visual.vector(0.33, 0.33, 0.5)
    soilCanvas.up = visual.vector(0,0,1)
    soilCanvas.forward = visual.vector(0.33, -0.33, -0.15)
    soilCanvas.range = (rectangularMesh.header.xMax - rectangularMesh.header.xMin) * 0.55
    layerLabel = visual.label(canvas = soilCanvas, height = h, pos=visual.vector(cX, cY, Zlabel))
    
    drawColorScale()
    drawSubSurface(True)

    #SLICE CANVAS
    sliceCanvas = visual.canvas(width = dx, height = dy, align="left")
    sliceCanvas.background = visual.color.white
    sliceCanvas.center = visual.vector(cX, cY+0.5, cZ-0.25)
    sliceCanvas.ambient = visual.vector(0.33, 0.33, 0.5)
    sliceCanvas.up = visual.vector(0,0,1)
    sliceCanvas.forward = visual.vector(0, -0.5, 0)
    sliceCanvas.range = (rectangularMesh.header.xMax - rectangularMesh.header.xMin) * 0.55
    
    sliceCanvas.caption = " *** COMMANDS ***\n\n 'r': run simulation \n 'p': pause "
    sliceCanvas.caption += "\n 'w': move up (soil layer) \n 's': move down (soil layer) "
    sliceCanvas.caption += "\n 'a': move left (soil slice) \n 'd': move right (soil slice) "
    sliceCanvas.caption += "\n 'c': colorscale range"
    sliceLabel = visual.label(canvas = sliceCanvas, height = h, pos=visual.vector(cX, cY, Zlabel))
    
    drawSlice(True)
    updateInterface()
    interface.bind('keydown', keyInput)


def drawColorScale():   
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
    
    if s == 's':
        if (visualizedLayer < C3DStructure.nrLayers-1):
            visualizedLayer += 1
    elif s == 'w':
        if (visualizedLayer > 0):
            visualizedLayer -= 1
             
    updateInterface()
    drawSubSurface(False)

def updateSlice(s):
    global visualizedSlice
    
    if s == 'a':
        if (visualizedSlice < (C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis)):
            visualizedSlice += C3DStructure.nrRectanglesInXAxis
    elif s == 'd':
        if (visualizedSlice > 0):
            visualizedSlice -= C3DStructure.nrRectanglesInXAxis
             
    updateInterface()
    drawSlice(False)
          
                
def keyInput(evt):
    global isPause
    s = evt.key
    if s == 'r':
        isPause = False
    elif s == 'p':
        isPause = True
        print ("Pause...")
    elif s == 's' or s == 'w':
        updateLayer(s)
    elif s == 'a' or s == 'd':
        updateSlice(s)
    elif s == "c":
        if (isPause):
            updateColorScale()


def getNewRectangle(myColor, myCanvas, v):
    vert = []
    for i in range(C3DStructure.nrVerticesPerRectangle):
        vert.append(visual.vertex(pos = visual.vector(v[i,0], v[i,1], v[i,2] * rectangularMesh.header.magnify)))
        vert[i].color = myColor
                     
    newRectangle = visual.quad(canvas = myCanvas, vs=[vert[0],vert[1],vert[3],vert[2]])
    return newRectangle 

def drawSlice(isFirst):
    global sliceRectangles
    
    posY = (visualizedSlice / C3DStructure.nrRectanglesInXAxis) * C3DStructure.gridStep + C3DStructure.gridStep*0.5
    sliceLabel.text = "Degree of saturation, slice at " + format(posY*100,".1f")+"cm"
    
    for z in range(C3DStructure.nrLayers):
        for x in range(C3DStructure.nrRectanglesInXAxis):
            index = visualizedSlice + x + (z * C3DStructure.nrRectangles)
            i = z * C3DStructure.nrRectanglesInXAxis + x
            c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
            myColor = visual.vector(c[0], c[1], c[2])
            
            if (isFirst):
                vertices = copy(rectangularMesh.C3DRM[visualizedSlice + x].v)
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
    
    depth = soil.depth[visualizedLayer] * 100
    layerLabel.text = "Degree of saturation, layer at " + format(depth,".1f")+"cm"
        
    for i in range(C3DStructure.nrRectangles):
        index = visualizedLayer * C3DStructure.nrRectangles + i
        c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
        myColor = visual.vector(c[0], c[1], c[2])
        
        if (isFirst):
            newRectangle = getNewRectangle(myColor, soilCanvas, rectangularMesh.C3DRM[i].v)
            subSurfaceRectangles.append(newRectangle)
        else:
            subSurfaceRectangles[i].v0.color = myColor
            subSurfaceRectangles[i].v1.color = myColor
            subSurfaceRectangles[i].v2.color = myColor  
            subSurfaceRectangles[i].v3.color = myColor
    
def updateInterface():       
    timeLabel.text = "Time: " + str(int(waterBalance.totalTime)) + " [s]"
    precLabel.text = "Precipitation: " + str(waterBalance.currentPrec) + " [mm/hour]"
    irrLabel.text = "Irrigation: " + str(waterBalance.currentIrr) + " [l/hour]"
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
      
    