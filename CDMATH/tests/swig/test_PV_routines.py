#!/usr/bin/env python
# -*-coding:utf-8 -*-

from PV_routines import *
import cdmath

#Meshes and fields initialisation
#================================

#cell field on 2D structured mesh
M1 = cdmath.Mesh(0.0, 1.0, 10, 0., 1., 5)
fieldName1="testfield1"
field1 = cdmath.Field(fieldName1, cdmath.CELLS, M1, 1)
for j in range(field1.getNumberOfComponents()):
    for i in range(field1.getNumberOfElements()):
        field1[i, j] = i + j

fileNameVTK1 = "2D_structured_cell_field"
field1.writeVTK(fileNameVTK1)

#node field on 2D unstructured mesh
M2 = cdmath.Mesh("meshSquare.med")
fieldName2="testfield2"
field2 = cdmath.Field(fieldName2, cdmath.NODES, M2, 1)
for j in range(field2.getNumberOfComponents()):
    for i in range(field2.getNumberOfElements()):
        field2[i, j] = i + j

fileNameVTK2 = "2D_unstructured_node_field"
field2.writeVTK(fileNameVTK2)

#node field on 3D unstructured mesh
M3 = cdmath.Mesh("meshCube.med")
fieldName3="testfield3"
field3 = cdmath.Field(fieldName3, cdmath.NODES, M3, 1)
for j in range(field3.getNumberOfComponents()):
    for i in range(field3.getNumberOfElements()):
        field3[i, j] = i + j

fileNameVTK3 = "3D_unstructured_node_field"
field3.writeVTK(fileNameVTK3)

#node field on sphere with unstructured mesh
M4 = cdmath.Mesh("meshSphere.med")
fieldName4="testfield4"
field4 = cdmath.Field(fieldName4, cdmath.NODES, M4, 1)
for j in range(field4.getNumberOfComponents()):
    for i in range(field4.getNumberOfElements()):
        field4[i, j] = i + j

fileNameVTK4 = "Sphere_unstructured_node_field"
field4.writeVTK(fileNameVTK4)

#cell field on 3D structured mesh
M5 = cdmath.Mesh(0.0, 1.0, 4, 0.0, 1.0, 4, 0.0, 1.0, 4)
fieldName5="testfield5"
field5 = cdmath.Field(fieldName5, cdmath.CELLS, M5, 1)
for j in range(field5.getNumberOfComponents()):
    for i in range(field5.getNumberOfElements()):
        field5[i, j] = i + j

fileNameVTK5 = "3D_structured_cell_field"
field5.writeVTK(fileNameVTK5)

#2D tests
#===========================================
point1=[1.,0.,0.]
point2=[0.,1.,0.]
resolution=100

outputFileName="Extract_PV_over_line_"+fileNameVTK1+".csv"
Extract_PV_data_over_line_to_txt_file('2D_structured_cell_field_0.vtu', outputFileName, point1, point2, resolution)
print( "Extract_VTK_over_line ok")

point=[0.5,0.5,0]
normal=[1,1,0]
outputFileName="Slice_PV_data_to_txt_file_"+fileNameVTK3+".csv"
Slice_PV_data_to_txt_file(fileNameVTK3+'_0.vtu', outputFileName, point, normal,resolution )
print( "Slice_PV_data_to_txt_file ok")

#outputFileName="Slice_field_data_to_txt_file"+fileNameVTK4+".csv"
#Slice_PV_field_data_to_txt_file(field4, outputFileName, point, normal,resolution)
#print "Slice_field_data_to_txt_file ok"

outputFileName="Clip_PV_data_to_VTK_"+fileNameVTK5
inputFileName="Clip_VTK_data_to_VTK_"+fileNameVTK5
Save_PV_data_to_picture_file(inputFileName+'_0.vtu',fieldName5,'CELLS',outputFileName)
print( "Save_PV_Clip_data_to_picture_file ok")

outputFileName="Slice_PV_data_to_VTK_"+fileNameVTK5
inputFileName="Slice_VTK_data_to_VTK_"+fileNameVTK5
Save_PV_data_to_picture_file(inputFileName+'_0.vtu',fieldName5,'CELLS',outputFileName)
print( "Save_PV_Slice_data_to_picture_file ok")

outputFileName="Save_PV_data_to_picture_file_"+fileNameVTK2
Save_PV_data_to_picture_file(fileNameVTK2+'_0.vtu',fieldName2,'NODES',outputFileName)
print( "Save_PV_data_to_picture_file " + fileNameVTK2+ " ok")

outputFileName="Save_PV_data_to_picture_file_"+fileNameVTK3
Save_PV_data_to_picture_file(fileNameVTK3+'_0.vtu',fieldName3,'NODES',outputFileName)
print( "Save_PV_data_to_picture_file " + fileNameVTK3+ " ok")

outputFileName="Save_PV_data_to_picture_file_"+fileNameVTK4
Save_PV_data_to_picture_file(fileNameVTK4+'_0.vtu',fieldName4,'NODES',outputFileName)
print( "Save_PV_data_to_picture_file " + fileNameVTK4+ " ok")

outputFileName="Save_PV_data_to_picture_file_"+fileNameVTK5
Save_PV_data_to_picture_file(fileNameVTK5+'_0.vtu',fieldName5,'CELLS',outputFileName)
print( "Save_PV_data_to_picture_file " + fileNameVTK5+ " ok")
