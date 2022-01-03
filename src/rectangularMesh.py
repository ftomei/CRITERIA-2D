# tin.py

from math import fabs, sqrt
import numpy as np
from copy import copy
from dataStructures import *
from enum import Enum


class Neighbours(Enum):
    Right = 1
    Up = 2
    Left = 3
    Down = 4


class CheaderRM():
    xMin = NODATA
    xMax = NODATA
    yMin = NODATA
    yMax = NODATA
    zMin = NODATA
    zMax = NODATA
    dz = NODATA
    magnify = NODATA


class Crectangle:
    def __init__(self, index, v=np.zeros((C3DStructure.nrVerticesPerRectangle, C3DStructure.nrDimensions), float)):
        self.index = index
        self.v = copy(v)
        self.centroid = getCentroid3D(self.v)
        self.area = getArea2D(self.v)

        self.neighbours, self.boundarySides, self.boundarySlopes = boundaryProperties(self.index, self.v, self.centroid)

        if NOLINK in self.neighbours:
            self.isBoundary = True
            self.boundarySide = np.max(self.boundarySides)
            self.boundarySlope = np.max(self.boundarySlopes)
        else:
            self.isBoundary = False
            self.boundarySide = NODATA
            self.boundarySlope = NODATA


# global structures
header = CheaderRM()
C3DRM = []


def rectangularMeshCreation():
    index = 0
    dx = C3DStructure.gridWidth * 0.5
    dy = C3DStructure.gridHeight * 0.5
    for y in np.arange(-dy, dy, C3DStructure.gridStep):
        for x in np.arange(-dx, dx, C3DStructure.gridStep):
            rectangle = Crectangle(index, getRectangleVertices(x, y))
            C3DRM.append(copy(rectangle))
            C3DStructure.totalArea += rectangle.area
            index += 1


def getRectangleVertices(x, y):
    v = np.zeros((C3DStructure.nrVerticesPerRectangle, C3DStructure.nrDimensions), float)
    dz = min(1.0, fabs(x)) * C3DStructure.slope
    dz2 = min(1.0, fabs(x + C3DStructure.gridStep)) * C3DStructure.slope
    v[0] = [x, y, C3DStructure.z0 - dz]
    v[1] = [x + C3DStructure.gridStep, y, C3DStructure.z0 - dz2]
    v[2] = [x, y + C3DStructure.gridStep, C3DStructure.z0 - dz]
    v[3] = [x + C3DStructure.gridStep, y + C3DStructure.gridStep, C3DStructure.z0 - dz2]
    return v


def boundaryProperties(index, vertices, centroid):
    neighbours = []
    boundarySides = []
    boundarySlopes = []
    greater_than = lambda a, b: a >= b
    less_than = lambda a, b: a <= b
    plus = lambda a, b: a + b
    minus = lambda a, b: a - b
    neighbours_map = {
        Neighbours.Right: {
            "index": index + 1,
            "limit": C3DStructure.gridWidth*0.5,
            "comparator": less_than,
            "vertices": [vertices[1], vertices[3]],
            "boundary_coordinate": 0,
            "operator": plus
        },
        Neighbours.Up: {
            "index": index + C3DStructure.nrRectanglesInXAxis,
            "limit": C3DStructure.gridHeight*0.5,
            "comparator": less_than,
            "vertices": [vertices[3], vertices[2]],
            "boundary_coordinate": 1,
            "operator": plus
        },
        Neighbours.Left: {
            "index": index - 1,
            "limit": -C3DStructure.gridWidth*0.5,
            "comparator": greater_than,
            "vertices": [vertices[2], vertices[0]],
            "boundary_coordinate": 0,
            "operator": minus
        },
        Neighbours.Down: {
            "index": index - C3DStructure.nrRectanglesInXAxis,
            "limit": -C3DStructure.gridHeight*0.5,
            "comparator": greater_than,
            "vertices": [vertices[0], vertices[1]],
            "boundary_coordinate": 1,
            "operator": minus
        },
    }
    for key, value in neighbours_map.items():
        vertex_to_check = neighbours_map[key]["vertices"][0]
        coordinate_to_check = neighbours_map[key]["boundary_coordinate"]
        limit_to_check = neighbours_map[key]["limit"]
        move_coordinate = neighbours_map[key]["operator"]
        check_limit = neighbours_map[key]["comparator"]

        if check_limit(move_coordinate(
                vertex_to_check[coordinate_to_check], C3DStructure.gridStep),
                limit_to_check):
            neighbours.append(neighbours_map[key]["index"])
            boundarySides.append(NODATA)
            boundarySlopes.append(NODATA)
        else:
            neighbours.append(NOLINK)

            # compute boundary side
            boundarySides.append(distance3D(neighbours_map[key]["vertices"][0], neighbours_map[key]["vertices"][1]))
            boundaryPoint = (neighbours_map[key]["vertices"][0] + neighbours_map[key]["vertices"][1]) * 0.5

            # compute slope
            dz = centroid[2] - boundaryPoint[2]
            dxy = distance2D(centroid, boundaryPoint)
            boundarySlopes.append(dz / dxy)
    return neighbours, boundarySides, boundarySlopes


# TODO to check
def magnitude(v):
    return np.sqrt(v.dot(v))


# TODO to check
def getArea(v):
    return 0.5 * magnitude(np.cross(v[1] - v[0], v[2] - v[0]))


def getCentroid3D(vertices):
    return np.array(
        [(vertices[0][0] + vertices[1][0]) / 2,
         (vertices[0][1] + vertices[2][1]) / 2,
         (vertices[0][2] + vertices[1][2]) / 2])


def getArea2D(vertices):
    x_side = vertices[1][0] - vertices[0][0]
    y_side = vertices[2][1] - vertices[0][1]
    return x_side * y_side


def getHeader(rectangleList):
    myHeader = CheaderRM()
    myHeader.xMin = rectangleList[0].centroid[0]
    myHeader.yMin = rectangleList[0].centroid[1]
    myHeader.zMin = rectangleList[0].centroid[2]
    myHeader.xMax = myHeader.xMin
    myHeader.yMax = myHeader.yMin
    myHeader.zMax = myHeader.zMin

    for i in range(1, len(rectangleList)):
        x = rectangleList[i].centroid[0]
        y = rectangleList[i].centroid[1]
        z = rectangleList[i].centroid[2]
        myHeader.xMin = min(myHeader.xMin, x)
        myHeader.yMin = min(myHeader.yMin, y)
        myHeader.zMin = min(myHeader.zMin, z)
        myHeader.xMax = max(myHeader.xMax, x)
        myHeader.yMax = max(myHeader.yMax, y)
        myHeader.zMax = max(myHeader.zMax, z)

    halfStep = C3DStructure.gridStep * 0.5
    myHeader.xMin -= halfStep
    myHeader.xMax += halfStep
    myHeader.yMin -= halfStep
    myHeader.yMax += halfStep

    dx = myHeader.xMax - myHeader.xMin
    dy = myHeader.yMax - myHeader.yMin
    myHeader.magnify = max(1., 1. / sqrt(max(dx, dy)))
    return myHeader


def distance2D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    return sqrt(dx * dx + dy * dy)


def distance3D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    dz = fabs(v1[2] - v2[2])
    return sqrt(dx * dx + dy * dy + dz * dz)


def getAdjacentVertices(rectangle1, rectangle2):
    commonVertices = []
    for i in rectangle1.v:
        for j in rectangle2.v:
            if np.allclose(i, j, rtol=1e-05):
                commonVertices.append(i)
    if len(commonVertices) > 2:
        raise Exception('Too many common vertices')
    return commonVertices[0], commonVertices[1]


def getAdjacentSide(i, j):
    rectangle1 = C3DRM[i]
    rectangle2 = C3DRM[j]
    v1, v2 = getAdjacentVertices(rectangle1, rectangle2)
    return distance2D(v1, v2)
