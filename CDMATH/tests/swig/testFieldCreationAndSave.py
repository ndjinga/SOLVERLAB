#!/usr/bin/env python3
# -*-coding:utf-8 -*

from math import sqrt
import cdmath
import medcoupling as mc

print("Loading a triangular mesh of a 2D square")
filename = "./meshSquare"
M=cdmath.Mesh(filename+".med")

#Extract groups in the mesh
print("Checking boundary group names")
boundaryFaceGroupNames=M.getNameOfFaceGroups()
boundaryNodeGroupNames=M.getNameOfNodeGroups()
print(len(boundaryFaceGroupNames), " Boundary Face Group detected : ", boundaryFaceGroupNames)
print(len(boundaryNodeGroupNames), " Boundary Node Group detected : ", boundaryNodeGroupNames)

assert(len(boundaryFaceGroupNames)==5)
assert(len(boundaryNodeGroupNames)==5)

assert boundaryFaceGroupNames[4]=="Top"
assert boundaryFaceGroupNames[3]=="Right"
assert boundaryFaceGroupNames[2]=="Left"
assert boundaryFaceGroupNames[1]=="Bottom"
assert boundaryFaceGroupNames[0]=="Boundary"

assert boundaryNodeGroupNames[4]=="Top"
assert boundaryNodeGroupNames[3]=="Right"
assert boundaryNodeGroupNames[2]=="Left"
assert boundaryNodeGroupNames[1]=="Bottom"
assert boundaryNodeGroupNames[0]=="Boundary"

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

#### Delete mesh and still save field in a med file
m = mc.MEDCouplingCMesh()
x = mc.DataArrayDouble([0.,1.,2.])
y = mc.DataArrayDouble([0.,1.,2.])
m.setCoords(x, y)
m = m.buildUnstructured()
m.setName("mesh")

f = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
f.setMesh(m)
f.setName("F")

da = mc.DataArrayDouble([1,2,3,4])
f.setArray(da)
f.setTime(0.0,0,0)

# Maillage d'abord:
fName = "./michael.med"
mc.WriteUMesh(fName, m, True)

# Maintenant juste les champs:
ff = mc.MEDFileField1TS()
ff.setFieldNoProfileSBT(f)
ff.write(fName, 0)  

# Tue le maillage
m = 0
del m
import gc
gc.collect()  # Make sure Python interp has called mesh destructor ...
 
# Ecrit encore du champ:
da2 = da.deepCopy()
f.setTime(1.0,1,0)
da2 += 1
print(da2.getValues())
f.setArray(da2)

ff = mc.MEDFileField1TS()
ff.setFieldNoProfileSBT(f)  # le maillage n'existe plus, tant pis :-)
ff.write(fName, 0)  
ff.write(fName, 0)  # 1 append
