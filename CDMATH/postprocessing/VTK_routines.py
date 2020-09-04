#!/usr/bin/env python
# -*-coding:utf-8 -*-

import os
import numpy as np
import vtk
from vtk.util import numpy_support as npvtk 
# do I need to kill the pipeline?

def Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    probeLine = vtk.vtkLineSource()
    probeLine.SetPoint1(point1)
    probeLine.SetPoint2(point2)
    probeLine.SetResolution(resolution)
    
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(probeLine.GetOutputPort())
    probe.SetSourceData(reader.GetOutput())
    probe.Update()

    vtkarray = probe.GetOutput().GetPointData().GetArray(0) # or Slice1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array
    
def Extract_VTK_data_over_line_to_txt_file(inputFileName, outputFileName, point1, point2, resolution):

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")
   
def Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution):

    inputFileName = field.getName()#os.getcwd()+field.get_name()
    field.writeVTK(inputFileName)

    numpy_array = Extract_VTK_data_over_line_to_numpyArray(inputFileName+"_1.vtu", point1, point2, resolution)

    os.remove(inputFileName+"_1.vtu")
    return numpy_array

def Extract_field_data_over_line_to_txt_file(field, point1, point2, resolution, outputFileName):

    numpy_array = Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")

def Slice_VTK_data_to_numpyArray(inputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkFiltersCorePython.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.Update()

    vtkarray = cutter.GetOutput().GetPointData().GetArray(0)
    numpy_array = npvtk.vtk_to_numpy(vtkarray)
    
    return numpy_array

    
def Slice_VTK_data_to_txt_file(inputFileName, outputFileName,
                                           point, normal,
                                           resolution
                                           ):
    numpy_array =   Slice_VTK_data_to_numpyArray(inputFileName, point, normal, resolution )  
    
    np.savetxt(outputFileName, numpy_array, delimiter=" ")
    
     
def Slice_field_data_to_numpyArray(field,
                                   point, normal,
                                   resolution
                                   ):
    inputFileName = field.getName()
    field.writeVTK(inputFileName)
 
    numpy_array = Slice_VTK_data_to_numpyArray(inputFileName+"_1.vtu", point, normal, resolution)

    os.remove(inputFileName+"_1.vtu")
    return numpy_array

def Slice_field_data_to_txt_file(field, outputFileName,
                                        point, normal,
                                        resolution):
    numpy_array = Slice_field_data_to_numpyArray(field, point, normal, resolution)

    np.savetxt(outputFileName, numpy_array, delimiter=" ")

def Slice_VTK_data_to_VTK(inputFileName,
                             outputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkFiltersCorePython.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.Update()

    #Convert tht polydata structure générated by cutter into unstructured grid by triangulation
    triFilter = vtk.vtkDataSetTriangleFilter()
    triFilter.SetInputConnection(cutter.GetOutputPort())
    triFilter.Update()
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(triFilter.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Write()

def Clip_VTK_data_to_VTK(inputFileName,
                             outputFileName,
                                 point, normal,
                                 resolution
                                           ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    clipper = vtk.vtkClipDataSet()
    clipper.SetClipFunction(plane)
    clipper.SetInputConnection(reader.GetOutputPort())
    clipper.Update()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(clipper.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Write()

def Save_VTK_data_to_picture_file(inputFileName, field_name,
                             node_or_cell, outputFileName
                             ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(inputFileName)
    reader.Update()
    
    if node_or_cell== 'CELLS':
        reader.CellArrayStatus = [field_name]
        reader.GetOutput().GetCellData().SetActiveScalars(field_name)
    elif node_or_cell== 'NODES':
        reader.PointArrayStatus = [field_name]
        reader.GetOutput().GetPointData().SetActiveScalars(field_name)
    else:
        raise ValueError("unknown type : should be CELLS or NODES")

#-------------------------------------------------------------------------------    
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    scalarBar=vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(mapper.GetLookupTable())
    scalarBar.SetTitle(field_name)
    
    mapper.SetScalarRange(reader.GetOutput().GetScalarRange())

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    
    ren.AddViewProp(actor)
    ren.AddActor2D(scalarBar);
        
    renWin.Render()
    
    image = vtk.vtkWindowToImageFilter()
    image.SetInput(renWin)
    image.ReadFrontBufferOff()
    image.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(outputFileName+".png")
    writer.SetInputConnection(image.GetOutputPort())
    writer.Write()
    
