import config
import numpy as np
import rectangleMesh
from dataStructures import *

def rectangleMeshIndexTest():
    for i in range(C3DStructure.nrRectangles):
        assert i == rectangleMesh.C3DRM[i].index, "index of the {}th element should be {} instead of{}".format(i, rectangleMesh.C3DRM[i].index, i)

def rectangleMeshVertexTest():
    originVertex = np.array([
        [C3DStructure.gridOrigin, C3DStructure.gridOrigin, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin, C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin]])
    
    lastWidthVertex = np.array([
        [C3DStructure.gridWidth - C3DStructure.gridStep, C3DStructure.gridOrigin, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth, C3DStructure.gridOrigin, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth - C3DStructure.gridStep, C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth, C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridOrigin]])
    
    lastHeightVertex = np.array([
        [C3DStructure.gridOrigin, C3DStructure.gridHeight - C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridHeight - C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin, C3DStructure.gridHeight, C3DStructure.gridOrigin],
        [C3DStructure.gridOrigin + C3DStructure.gridStep, C3DStructure.gridHeight, C3DStructure.gridOrigin]])
    
    lastVertex = np.array([
        [C3DStructure.gridWidth - C3DStructure.gridStep, C3DStructure.gridHeight - C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth, C3DStructure.gridHeight - C3DStructure.gridStep, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth - C3DStructure.gridStep, C3DStructure.gridHeight, C3DStructure.gridOrigin],
        [C3DStructure.gridWidth, C3DStructure.gridHeight, C3DStructure.gridOrigin]])

    comparison = rectangleMesh.C3DRM[0].v == originVertex
    assert comparison.all(), "originVertex should be \n{}\ninstead of\n{}".format(originVertex, rectangleMesh.C3DRM[0].v)

    comparison = rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].v == lastWidthVertex
    assert comparison.all(), "lastWidthVertex should be \n{}\ninstead of\n{}".format(lastWidthVertex, rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].v)

    comparison = rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].v == lastHeightVertex
    assert comparison.all(), "lastHeightVertex should be \n{}\ninstead of\n{}".format(lastHeightVertex, rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].v)

    comparison = rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].v == lastVertex
    assert comparison.all(), "lastVertex should be \n{}\ninstead of\n{}".format(lastVertex, rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].v)

def rectangleMeshCentroidTest():
    originRectangleCentroid = np.array(
        [C3DStructure.gridOrigin + C3DStructure.gridStep / 2,  C3DStructure.gridOrigin + C3DStructure.gridStep / 2, C3DStructure.gridOrigin]
        )
    
    lastWidthRectangleCentroid = np.array(
        [C3DStructure.gridWidth - C3DStructure.gridStep / 2,  C3DStructure.gridOrigin + C3DStructure.gridStep / 2, C3DStructure.gridOrigin]
        )
    
    lastHeightRectangleCentroid = np.array(
        [C3DStructure.gridOrigin + C3DStructure.gridStep / 2,  C3DStructure.gridHeight - C3DStructure.gridStep / 2, C3DStructure.gridOrigin]
        )
    
    lastRectangleCentroid = np.array(
        [C3DStructure.gridWidth - C3DStructure.gridStep / 2,  C3DStructure.gridHeight - C3DStructure.gridStep / 2, C3DStructure.gridOrigin]
        )
    
    assert np.allclose(
        rectangleMesh.C3DRM[0].centroid, 
        originRectangleCentroid, 
        rtol = 1e-05), "originRectangleCentroid should be \n{}\ninstead of\n{}".format(originRectangleCentroid, rectangleMesh.C3DRM[0].centroid)

    assert np.allclose(
        rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].centroid, 
        lastWidthRectangleCentroid, 
        rtol = 1e-05), "lastWidthRectangleCentroid should be \n{}\ninstead of\n{}".format(lastWidthRectangleCentroid, rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].centroid)

    assert np.allclose(
        rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].centroid, 
        lastHeightRectangleCentroid, 
        rtol = 1e-05), "lastHeightRectangleCentroid should be \n{}\ninstead of\n{}".format(lastHeightRectangleCentroid, rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].centroid)

    assert np.allclose(
        rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].centroid, 
        lastRectangleCentroid, 
        rtol = 1e-05), "lastRectangleCentroid should be \n{}\ninstead of\n{}".format(lastRectangleCentroid, rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].centroid)

def rectangleMeshAreaTest():
    area = C3DStructure.gridStep * C3DStructure.gridStep
    allowed_error = 1e-05

    comparison = abs(rectangleMesh.C3DRM[0].area - area) <=  allowed_error
    assert comparison, "originRectangleArea should be \n{}\ninstead of\n{}".format(area, rectangleMesh.C3DRM[0].area)

    comparison = abs(rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].area - area) <=  allowed_error
    assert comparison, "lastWidthRectangleArea should be \n{}\ninstead of\n{}".format(area, rectangleMesh.C3DRM[C3DStructure.nrRectanglesInXAxis - 1].area)

    comparison = abs(rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].area - area) <=  allowed_error
    assert comparison, "lastHeightRectangleArea should be \n{}\ninstead of\n{}".format(area, rectangleMesh.C3DRM[C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis].area)

    comparison = abs(rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].area - area) <=  allowed_error
    assert comparison, "lastRectangleArea should be \n{}\ninstead of\n{}".format(area, rectangleMesh.C3DRM[C3DStructure.nrRectangles - 1].area)

def rectangleMeshBoundaryTest():
    desiredBoundarySide = 0.01
    desiredBoundarySlope = 0.0

    originIndex = 0
    originIndexNeighbours = np.array([1, C3DStructure.nrRectanglesInXAxis])

    lastWidthIndex = C3DStructure.nrRectanglesInXAxis - 1
    lastWidthIndexNeighbours = np.array([2*C3DStructure.nrRectanglesInXAxis - 1, C3DStructure.nrRectanglesInXAxis - 2])
    
    lastHeightIndex = C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis
    lastHeightIndexNeighbours = np.array([C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis + 1, C3DStructure.nrRectangles - 2 * C3DStructure.nrRectanglesInXAxis])
    
    lastIndex = C3DStructure.nrRectangles - 1
    lastIndexNeighbours = np.array([C3DStructure.nrRectangles - 2, C3DStructure.nrRectangles - C3DStructure.nrRectanglesInXAxis - 1])
    
    secondRowFirstColumnIndex = C3DStructure.nrRectanglesInXAxis
    secondRowFirstColumnIndexNeighbours = np.array([C3DStructure.nrRectanglesInXAxis + 1, 2*C3DStructure.nrRectanglesInXAxis, 0])

    secondRowSecondColumnIndex = C3DStructure.nrRectanglesInXAxis + 1
    secondRowSecondColumnIndexNeighbours = np.array([C3DStructure.nrRectanglesInXAxis + 2, 2*C3DStructure.nrRectanglesInXAxis + 1, C3DStructure.nrRectanglesInXAxis, 1])

    indexes = [originIndex, lastWidthIndex, lastHeightIndex, lastIndex, secondRowFirstColumnIndex, secondRowSecondColumnIndex]
    neighbours = [originIndexNeighbours, lastWidthIndexNeighbours, lastHeightIndexNeighbours, lastIndexNeighbours, secondRowFirstColumnIndexNeighbours, secondRowSecondColumnIndexNeighbours]
    for i in range(len(indexes)):
        comparison = rectangleMesh.C3DRM[indexes[i]].neighbours == neighbours[i]
        assert comparison.all(), "Neighbours of the {}th element should be \n{}\ninstead of\n{}".format(indexes[i], neighbours[i], rectangleMesh.C3DRM[indexes[i]].neighbours)

        if indexes[i] != secondRowSecondColumnIndex:
            assert np.allclose(
                rectangleMesh.C3DRM[indexes[i]].boundarySide, 
                desiredBoundarySide, 
                rtol = 1e-05), "Boundary Slide of the {}th element should be \n{}\ninstead of\n{}".format(indexes[i], desiredBoundarySide, rectangleMesh.C3DRM[indexes[i]].boundarySide)

            assert np.allclose(
                rectangleMesh.C3DRM[indexes[i]].boundarySlope, 
                desiredBoundarySlope, 
                rtol = 1e-05), "Boundary Slope of the {}th element should be \n{}\ninstead of\n{}".format(indexes[i], desiredBoundarySlope, rectangleMesh.C3DRM[indexes[i]].boundarySlope)

def totalAreaTest():
    totalArea = C3DStructure.gridWidth * C3DStructure.gridHeight
    allowed_error = 1e-05

    comparison = abs(C3DStructure.totalArea - totalArea) <=  allowed_error
    assert comparison, "totalArea should be \n{}\ninstead of\n{}".format(totalArea, rectangleMesh.C3DRM[0].area)

def rectangleMeshCreationTest():
    rectangleMeshIndexTest()
    rectangleMeshVertexTest()
    rectangleMeshCentroidTest()
    rectangleMeshAreaTest()
    rectangleMeshBoundaryTest()

    totalAreaTest()
    
if __name__ == "__main__":
    rectangleMesh.rectangleMeshCreation()
    rectangleMeshCreationTest()
    print("Everything passed")