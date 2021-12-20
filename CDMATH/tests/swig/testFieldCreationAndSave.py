#!/usr/bin/env python3
# -*-coding:utf-8 -*

from math import sqrt
import cdmath
import medcoupling as mc

print("Loading a triangular mesh of a 2D square")
filename = "./meshSquare"
M=cdmath.Mesh(filename+".med", "Mesh_1", 0)

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

# Create solid temperature fields
temperature_field_cells = cdmath.Field("Solid temperature", cdmath.CELLS, M, 1)
temperature_field_nodes = cdmath.Field("Solid temperature", cdmath.NODES, M, 1)
Tin  = 330.
Tout = 300.

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

# Create boundary fields
## bottom
Mbottom=M.getBoundaryGroupMesh ( "Bottom" )
Mbottom.writeMED(filename, False)
Mbottom.writeVTK(filename)
temperature_bottom_cells=cdmath.Field("Bottom temperature",cdmath.CELLS,Mbottom)
for i in range(temperature_bottom_cells.getNumberOfElements()):
	temperature_bottom_cells[i]=Tout
temperature_bottom_cells.writeMED(filename,False)
temperature_bottom_cells.writeVTK(filename,True)
temperature_bottom_cells.writeCSV(filename)
temperature_bottom_nodes=cdmath.Field("Bottom temperature",cdmath.NODES,Mbottom)
for i in range(temperature_bottom_nodes.getNumberOfElements()):
	temperature_bottom_nodes[i]=Tout
temperature_bottom_nodes.writeMED(filename,False)
temperature_bottom_nodes.writeVTK(filename,True)
temperature_bottom_nodes.writeCSV(filename)
##top
Mtop=M.getBoundaryGroupMesh ( "Top" )
Mtop.writeMED(filename, False)
Mtop.writeVTK(filename)
temperature_top_cells=cdmath.Field("Top temperature",cdmath.CELLS,Mtop)
for i in range(temperature_top_cells.getNumberOfElements()):
	temperature_top_cells[i]=Tout
temperature_top_cells.writeMED(filename,False)
temperature_top_cells.writeVTK(filename,True)
temperature_top_cells.writeCSV(filename)
temperature_top_nodes=cdmath.Field("Top temperature",cdmath.NODES,Mtop)
for i in range(temperature_top_nodes.getNumberOfElements()):
	temperature_top_nodes[i]=Tout
temperature_top_nodes.writeMED(filename,False)
temperature_top_nodes.writeVTK(filename,True)
temperature_top_nodes.writeCSV(filename)
##Left
Mleft=M.getBoundaryGroupMesh ( "Left" )
Mleft.writeMED(filename, False)
Mleft.writeVTK(filename)
temperature_left_cells=cdmath.Field("Left temperature",cdmath.CELLS,Mleft)
for i in range(temperature_left_cells.getNumberOfElements()):
	temperature_left_cells[i]=Tout
temperature_left_cells.writeMED(filename,False)
temperature_left_cells.writeVTK(filename,True)
temperature_left_cells.writeCSV(filename)
temperature_left_nodes=cdmath.Field("Left temperature",cdmath.NODES,Mleft)
for i in range(temperature_left_nodes.getNumberOfElements()):
	temperature_left_nodes[i]=Tout
temperature_left_nodes.writeMED(filename,False)
temperature_left_nodes.writeVTK(filename,True)
temperature_left_nodes.writeCSV(filename)
##Right
Mright=M.getBoundaryGroupMesh ( "Right" )
Mright.writeMED(filename, False)
Mright.writeVTK(filename)
temperature_right_cells=cdmath.Field("Right temperature",cdmath.CELLS,Mright)
for i in range(temperature_right_cells.getNumberOfElements()):
	temperature_right_cells[i]=Tout
temperature_right_cells.writeMED(filename,False)
temperature_right_cells.writeVTK(filename,True)
temperature_right_cells.writeCSV(filename)
temperature_right_nodes=cdmath.Field("Right temperature",cdmath.NODES,Mright)
for i in range(temperature_right_nodes.getNumberOfElements()):
	temperature_right_nodes[i]=Tout
temperature_right_nodes.writeMED(filename,False)
temperature_right_nodes.writeVTK(filename,True)
temperature_right_nodes.writeCSV(filename)

# Create fluid temperature fields
temperature_field_cells = cdmath.Field("Fluid temperature", cdmath.CELLS, M, 1)
temperature_field_nodes = cdmath.Field("Fluid temperature", cdmath.NODES, M, 1)
Tfluid  = 300.

for i in range(nbCells):
    x = M.getCell(i).x()
    y = M.getCell(i).y()
    temperature_field_cells[i] = Tfluid

temperature_field_cells.writeMED(filename, False)

for i in range(nbNodes):
    x = M.getNode(i).x()
    y = M.getNode(i).y()
    temperature_field_nodes[i] = Tfluid

temperature_field_nodes.writeMED(filename, False)

# Create fluid enthalpy fields
enthalpy_field_cells = cdmath.Field("Fluid enthalpy", cdmath.CELLS, M, 1)
enthalpy_field_nodes = cdmath.Field("Fluid enthalpy", cdmath.NODES, M, 1)
Hin  = 1.6e6
Hout = 2.e6

for i in range(nbCells):
    x = M.getCell(i).x()
    y = M.getCell(i).y()
    distance = sqrt( (x - xcentre) * (x - xcentre) + (y - ycentre) * (y - ycentre) )
    if distance < radius:
        enthalpy_field_cells[i] = Hin
    else:
        enthalpy_field_cells[i] = Hout

enthalpy_field_cells.writeMED(filename, False)

for i in range(nbNodes):
    x = M.getNode(i).x()
    y = M.getNode(i).y()
    distance = sqrt( (x - xcentre) * (x - xcentre) + (y - ycentre) * (y - ycentre) )
    if distance < radius:
        enthalpy_field_nodes[i] = Hin
    else:
        enthalpy_field_nodes[i] = Hout

enthalpy_field_nodes.writeMED(filename, False)

# Create heat power fields
heat_field_cells = cdmath.Field("Heat power", cdmath.CELLS, M, 1)
heat_field_nodes = cdmath.Field("Heat power", cdmath.NODES, M, 1)
phi = 1.e7

for i in range(nbCells):
    x = M.getCell(i).x()
    y = M.getCell(i).y()
    heat_field_cells[i] = phi

heat_field_cells.writeMED(filename, False)

for i in range(nbNodes):
    x = M.getNode(i).x()
    y = M.getNode(i).y()
    heat_field_nodes[i] = phi

heat_field_nodes.writeMED(filename, False)

# Create pressure and velocity fields for the wave equation on CELLS
p0=155e5#reference pressure in a pressurised nuclear vessel

pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, M, 1)
velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, M, 3)

for i in range(nbCells):
    x = M.getCell(i).x()
    y = M.getCell(i).y()
    distance = sqrt( (x - xcentre) * (x - xcentre) + (y - ycentre) * (y - ycentre) )

    velocity_field[i,0] = 0
    velocity_field[i,1] = 0
    velocity_field[i,2] = 0

    if distance < radius:
        pressure_field[i] = p0
        pass
    else:
        pressure_field[i] = p0/2
        pass
    pass

pressure_field.writeMED(filename, False)
velocity_field.writeMED(filename, False)

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
fName = "/tmp/michael.med"
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
da2 *= -1
print(da2.getValues())
f.setArray(da2)

ff = mc.MEDFileField1TS()
ff.setFieldNoProfileSBT(f)  # le maillage n'existe plus, tant pis :-)
ff.write(fName, 0)    # yes 0
