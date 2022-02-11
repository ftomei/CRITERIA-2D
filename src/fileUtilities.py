# fileUtilities

from math import fabs
from readDataFile import *
from dataStructures import NODATA, C3DStructure, C3DCells
import soil
import tkinter
import tkinter.filedialog
import criteria3D


def getStateFileName(isSave):
    root = tkinter.Tk()
    options = {'defaultextension': ".csv", 'filetypes': [("Comma separated values", ".csv")], 'initialdir': "data"}
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
        lastLayer = C3DStructure.nrLayers - 1
        # depth
        for layer in range(C3DStructure.nrLayers):
            f.write(str(soil.depth[layer]))
            if layer < lastLayer:
                f.write(",")
            else:
                f.write("\n")
        # water potential
        for i in range(C3DStructure.nrRectangles):
            for layer in range(C3DStructure.nrLayers):
                index = i + C3DStructure.nrRectangles * layer
                # matric potential [m]
                h = C3DCells[index].H - C3DCells[index].z
                if fabs(h) < 1E-12: h = 0.0
                f.write(str(h))
                if layer < lastLayer:
                    f.write(",")
                else:
                    f.write("\n")

