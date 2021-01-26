import waterBalance
import os
from dataStructures import *

exportIndeces = []
outputDirectory = os.path.join("data", "fondo_1", "output")
outputFile = os.path.join(outputDirectory, "output.csv")
#nrDetections = -1
heightSlice = C3DStructure.gridHeight * 0.5
oneTimestampPerRow = False

def createExportFile():
    if heightSlice == None:
        takeAll()
    else:
        takeSlice()
    
    if oneTimestampPerRow:
        header = "timestamp," + ",".join(map(lambda index: str(index), exportIndeces)) + "\n"
    else:
        header = "timestamp,x,y,z,Se,H\n"

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    with open(outputFile, "w") as f:
        f.write(header)


def takeSlice():
    offset = C3DStructure.nrRectanglesInYAxis / (C3DStructure.gridHeight / heightSlice)
    i = offset * C3DStructure.nrRectanglesInXAxis
    for layer in range(C3DStructure.nrLayers):
        for j in range(C3DStructure.nrRectanglesInXAxis):
            z = C3DStructure.nrRectangles * layer
            index = i + j + z
            exportIndeces.append(int(index))

def takeAll():
    for index in range(C3DStructure.nrCells):
        exportIndeces.append(index)


def takeScreenshot(timestamp):
    #global nrDetections
    
    #print("_______________________________")
    #print(waterBalance.totalTime)
    #print(int(waterBalance.totalTime / 60), nrDetections)
    #print("_______________________________")
    
    #if int(waterBalance.totalTime / 60) > nrDetections:
        #nrDetections += 1
        #anche tutto il resto andrebbe sotto l'if, se lo scommento

    if oneTimestampPerRow:
        row = str(int(timestamp))
        for index in exportIndeces:
            row += "," + '{:.3f}'.format(C3DCells[index].Se)
        row += "\n"

        with open(outputFile, "a") as f:
            f.write(row)
    else:
        for index in exportIndeces:
            if C3DCells[index].z != 0.0 and C3DCells[index].x >= 1.0:
                row = str(int(timestamp))
                row += "," + '{:.3f}'.format(C3DCells[index].x)
                row += "," + '{:.3f}'.format(C3DCells[index].y)
                row += "," + '{:.3f}'.format(C3DCells[index].z)
                row += "," + '{:.3f}'.format(C3DCells[index].Se)
                psi = (C3DCells[index].H - C3DCells[index].z) * 9.81    # water potential [kPa] equivalent to [centibar]
                row += "," + '{:.3f}'.format(psi)
                row += "\n"

                with open(outputFile, "a") as f:
                    f.write(row)
