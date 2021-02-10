#!/usr/bin/env python
# -*- coding: utf-8 -*-
   
import vtk
import VTKReader#To generate a med file from a vtk file

# ------------------------------------------------------------
# Create Boy's surface
# ------------------------------------------------------------
boy = vtk.vtkParametricBoy()
boySource = vtk.vtkParametricFunctionSource()
boySource.SetParametricFunction(boy)
boySource.SetScalarModeToModulus()

boyMapper = vtk.vtkPolyDataMapper()
boyMapper.SetInputConnection(boySource.GetOutputPort())
boyMapper.SetScalarRange(0, 2)
boyActor = vtk.vtkActor()
boyActor.SetMapper(boyMapper)
boyActor.SetPosition(8, -4, 0)
boyActor.SetScale(1.5, 1.5, 1.5)

boyTextMapper = vtk.vtkTextMapper()
boyTextMapper.SetInput("Boy")
boyTextMapper.GetTextProperty().SetJustificationToCentered()
boyTextMapper.GetTextProperty().SetVerticalJustificationToCentered()
boyTextMapper.GetTextProperty().SetColor(1, 0, 0)
boyTextMapper.GetTextProperty().SetFontSize(14)
boyTextActor = vtk.vtkActor2D()
boyTextActor.SetMapper(boyTextMapper)
boyTextActor.GetPositionCoordinate().SetCoordinateSystemToWorld()
boyTextActor.GetPositionCoordinate().SetValue(8, -6.5, 0)

# ------------------------------------------------------------
# Create the RenderWindow, Renderer and both Actors
# ------------------------------------------------------------
ren = vtk.vtkRenderer()
ren.AddViewProp(boyActor)
ren.AddViewProp(boyTextActor)

ren.SetBackground(0.7, 0.8, 1)
ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.3)

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(500, 500)

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin) 

iren.Initialize()
renWin.Render()

image = vtk.vtkWindowToImageFilter()
image.ReadFrontBufferOff()
image.SetInput(renWin)
image.Update()

# ------------------------------------------------------------
# Save picture
# ------------------------------------------------------------

writer = vtk.vtkPNGWriter()
writer.SetFileName("BoySurface"+".png")
writer.SetInputConnection(image.GetOutputPort())
writer.Write()

# ------------------------------------------------------------
# Save polydata shape
# ------------------------------------------------------------

writer = vtk.vtkXMLPolyDataWriter()
writer.SetInputData(boySource.GetOutput())
writer.SetFileName("BoySurface"+".vtp")
writer.Write()

# ------------------------------------------------------------
# Save unstructured grid
# ------------------------------------------------------------

appendFilter=vtk.vtkAppendFilter()
appendFilter.AddInputData(boySource.GetOutput())
appendFilter.Update()

unstructuredGrid=vtk.vtkUnstructuredGrid()#appendFilter.GetOutput()#
unstructuredGrid.DeepCopy(appendFilter.GetOutput())

writer=vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("BoySurface"+".vtu")
writer.SetInputData(unstructuredGrid)
#writer.SetDataModeToBinary()
#writer.SetDataModeToAppended()
#writer.EncodeAppendedDataOn()
writer.Write()

#Generate binary vtu file in append mode
#boySurfaceBisvtu = pvs.XMLUnstructuredGridReader(FileName=['./BoySurface.vtu'])
#pvs.SaveData('./BoySurface2.vtu', proxy=boySurfaceBisvtu, DataMode='Binary',EncodeAppendedData=1)

#Generate med file
vtu = VTKReader.VTURawReader('./BoySurface.vtu')
med = vtu.loadInMEDFileDS()
med.write("./BoySurface.med", 2)
