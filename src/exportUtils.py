import waterBalance
import os
from dataStructures import *

exportIndeces = []
outputFile = os.path.join("data", "fondo_1", "output", "output.csv")
#nrDetections = -1
heightSlice = 0.5
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

    if not os.path.exists(outputFile):
        open(outputFile, 'w').close()

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
        #row = str(waterBalance.totalTime)
        row = str(timestamp.value // 10 ** 9)
        for index in exportIndeces:
            row += "," + '{:.6f}'.format(C3DCells[index].Se)
        row += "\n"

        with open(outputFile, "a") as f:
            f.write(row)
    else:
        for index in exportIndeces:
            #row = str(waterBalance.totalTime)
            row = str(timestamp)
            row += "," + '{:.6f}'.format(C3DCells[index].x)
            row += "," + '{:.6f}'.format(C3DCells[index].y)
            row += "," + '{:.6f}'.format(C3DCells[index].z)
            row += "," + '{:.6f}'.format(C3DCells[index].Se)
            row += "," + '{:.6f}'.format(C3DCells[index].H)
            row += "\n"

            with open(outputFile, "a") as f:
                f.write(row)
