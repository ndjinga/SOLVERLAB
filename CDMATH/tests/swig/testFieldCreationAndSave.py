#!/usr/bin/env python3
# -*-coding:utf-8 -*

from math import sqrt
import cdmath

print("Loading a triangular mesh of a 2D square")
filename = "./meshSquare"
M=cdmath.Mesh(filename+".med")

print("Checking boundary group names")
boundaryFaceGroupNames=M.getNameOfFaceGroups()
boundaryNodeGroupNames=M.getNameOfNodeGroups()
print(len(boundaryFaceGroupNames), " Boundary Face Group detected : ", boundaryFaceGroupNames)
print(len(boundaryNodeGroupNames), " Boundary Node Group detected : ", boundaryNodeGroupNames)

assert(len(boundaryFaceGroupNames)==5)
assert(len(boundaryNodeGroupNames)==5)

assert boundaryFaceGroupNames[0]=="Top"
assert boundaryFaceGroupNames[1]=="Right"
assert boundaryFaceGroupNames[2]=="Left"
assert boundaryFaceGroupNames[3]=="Bottom"
assert boundaryFaceGroupNames[4]=="Boundary"

assert boundaryNodeGroupNames[0]=="Top"
assert boundaryNodeGroupNames[1]=="Right"
assert boundaryNodeGroupNames[2]=="Left"
assert boundaryNodeGroupNames[3]=="Bottom"
assert boundaryNodeGroupNames[4]=="Boundary"

#Extract domain sizes
xmin = M.getXMin()
xmax = M.getXMax()
ymin = M.getYMin()
ymax = M.getYMax()

radius = min(xmax-xmin,ymax-ymin)/4
xcentre = (xmax+xmin)/2
ycentre = (ymax+ymin)/2

nbCells = M.getNumberOfCells()
nbNodes = M.getNumberOfNodes()

# Create scalar fields
temperature_field_cells = cdmath.Field("Temperature", cdmath.CELLS, M, 1)
temperature_field_nodes = cdmath.Field("Temperature", cdmath.NODES, M, 1)
Tin  = 300
Tout = 400

for i in range(nbCells):
    x = M.getCell(i).x()
    y = M.getCell(i).y()
    distance = sqrt( (x - xcentre) * (x - xcentre) + (y - ycentre) * (y - ycentre) )
    if distance < radius:
        temperature_field_cells[i] = Tin
    else:
        temperature_field_cells[i] = Tout

temperature_field_cells.writeMED(filename, False)

for i in range(nbNodes):
    x = M.getNode(i).x()
    y = M.getNode(i).y()
    distance = sqrt( (x - xcentre) * (x - xcentre) + (y - ycentre) * (y - ycentre) )
    if distance < radius:
        temperature_field_nodes[i] = Tin
    else:
        temperature_field_nodes[i] = Tout

temperature_field_nodes.writeMED(filename, False)
