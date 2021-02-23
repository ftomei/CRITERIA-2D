#fileUtilities

from math import fabs
from readDataFile import *
from dataStructures import NODATA, C3DStructure, C3DCells
import soil
import criteria3D
    
    
def getStateFileName(isSave):
    root = tkinter.Tk()
    options = {}
    options['defaultextension'] = ".csv"
    options['filetypes'] = [("Comma separated values", ".csv")]
    options['initialdir'] = "data"
    if isSave:
        fileName = tkinter.filedialog.asksaveasfilename(**options)
    else:
        fileName = tkinter.filedialog.askopenfilename(**options)
    root.destroy()
    return fileName


def saveState():
    fileName = getStateFileName(True)
    if fileName != "":
        f = open(fileName, "w")
        lastLayer = C3DStructure.nrLayers-1
        #depth
        for layer in range(C3DStructure.nrLayers):
            f.write(str(soil.depth[layer]))
            if layer < lastLayer: f.write(",")
            else: f.write("\n")
        #water potential
        for i in range(C3DStructure.nrRectangles):
            for layer in range(C3DStructure.nrLayers):
                index = i + C3DStructure.nrRectangles * layer
                # matric potential [m]
                h = C3DCells[index].H - C3DCells[index].z
                if fabs(h) < 1E-12: h = 0.0
                f.write(str(h))
                if layer < lastLayer: f.write(",")
                else: f.write("\n")
                   
                   
def loadState(fileName):
    if (fileName == ""):
        fileName = getStateFileName(False)
    if (fileName == ""): return False
    state, isFileOk = readDataFile(fileName, 0, ",", False)
    if not isFileOk or (C3DStructure.nrRectangles + 1) > len(state): 
        print("*** Wrong state file!")
        return False
    #first row: depth
    depth = state[0]
    #surface
    for i in range(C3DStructure.nrRectangles):
        criteria3D.setMatricPotential(i, state[i+1][0])
    #subsurface
    for layer in range(1, C3DStructure.nrLayers):
        minDistance = 100.0
        l = NODATA
        for stateLayer in range(1, len(depth)):
            dz =  fabs(soil.depth[layer] - depth[stateLayer])
            if dz < minDistance:
                minDistance = dz
                l = stateLayer
        if (l != NODATA):
            for i in range(C3DStructure.nrRectangles):
                index = layer * C3DStructure.nrRectangles + i 
                criteria3D.setMatricPotential(index, state[i+1][l])
    return True

