import os
from dataStructures import *
import rectangularMesh
import pandas as pd

outputIndices = []
outputSurfaceIndices = []
outputFile = ""
heightSlice = C3DStructure.gridHeight * 0.5
oneTimestampPerRow = True


def createExportFile(outputPath):
    global outputFile
    outputFile = os.path.join(outputPath, "output.csv")

    if oneTimestampPerRow:
        outputPoints = pd.read_csv(os.path.join(outputPath, "output_points.csv"))
        outputIndicesString = takeSelected(outputPoints)
        header = "timestamp," + outputIndicesString + "\n"
    else:
        if heightSlice == 0:
            takeAll()
        else:
            takeSlice()
        header = "timestamp,x,y,z,Se,H\n"

    if not os.path.exists(outputFile):
        open(outputFile, 'w').close()

    with open(outputFile, "w") as f:
        f.write(header)


def takeSelected(outputPoints):
    outputIndicesStrings = []
    outputSurfaceIndices.clear()
    outputIndices.clear()
    for _, position in outputPoints.iterrows():
        x = position['x']
        y = position['y']
        depth = position['z']
        surfaceIndex = rectangularMesh.getSurfaceIndex(x, y)
        if surfaceIndex != NODATA:
            index = rectangularMesh.getCellIndex(x, y, depth)
            if index != NODATA:
                outputSurfaceIndices.append(surfaceIndex)
                outputIndices.append(index)
                outputIndicesStrings.append(f'z{depth}_y{y}_x{x}')
    return ",".join(outputIndicesStrings)


def takeSlice():
    outputIndices.clear()
    offset = C3DStructure.nrRectanglesInYAxis / (C3DStructure.gridHeight / heightSlice)
    i = offset * C3DStructure.nrRectanglesInXAxis
    for layer in range(C3DStructure.nrLayers):
        for j in range(C3DStructure.nrRectanglesInXAxis):
            z = C3DStructure.nrRectangles * layer
            index = i + j + z
            outputIndices.append(int(index))


def takeAll():
    outputIndices.clear()
    for index in range(C3DStructure.nrCells):
        outputIndices.append(index)


def takeScreenshot(timestamp):
    if oneTimestampPerRow:
        row = str(int(timestamp))
        for index in outputIndices:
            psi = (C3DCells[index].H - C3DCells[index].z) * 9.81  # water potential [kPa] equivalent to [centibar]
            row += "," + '{:.3f}'.format(psi)
        row += "\n"

        with open(outputFile, "a") as f:
            f.write(row)
    else:
        for index in outputIndices:
            if C3DCells[index].z != 0.0 and C3DCells[index].x >= 1.0:
                row = str(int(timestamp))
                row += "," + '{:.3f}'.format(C3DCells[index].x)
                row += "," + '{:.3f}'.format(C3DCells[index].y)
                row += "," + '{:.3f}'.format(C3DCells[index].z)
                row += "," + '{:.3f}'.format(C3DCells[index].Se)
                psi = (C3DCells[index].H - C3DCells[index].z) * 9.81  # water potential [kPa] equivalent to [centibar]
                row += "," + '{:.3f}'.format(psi)
                row += "\n"

                with open(outputFile, "a") as f:
                    f.write(row)
