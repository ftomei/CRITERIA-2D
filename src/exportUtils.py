import waterBalance
from dataStructures import *

exportIndeces = []
outputFile = "./data/output/slice_60s_heatmap.csv"
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
        header = "timestamp,x,z,value\n"

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
        row = str(timestamp)
        for index in exportIndeces:
            row += "," + str(C3DCells[index].Se)
        row += "\n"

        with open(outputFile, "a") as f:
            f.write(row)
    else:
        for index in exportIndeces:
            #row = str(waterBalance.totalTime)
            row = str(timestamp)
            row += "," + str(C3DCells[index].x)
            row += "," + str(C3DCells[index].z)
            row += "," + str(C3DCells[index].Se)
            row += "\n"

            with open(outputFile, "a") as f:
                f.write(row)
