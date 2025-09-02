#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.12.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/volatile/catB/esteban/Solverlab/INSTALL_release/share/examples')

###
### PARAVIS component
###

import pvsimple
pvsimple.ShowParaviewView()
# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from pvsimple import *


#### disable automatic camera reset on 'Show'
pvsimple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu = XMLUnstructuredGridReader(registrationName='WaveStaggered_2DCylinderDeflectionExactVelocity_0.vtu', FileName=['/volatile/catB/esteban/Solverlab/INSTALL_release/share/examples/WaveStaggered_2DCylinderDeflectionExactVelocity_0.vtu'])

# Properties modified on waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# get the material library
materialLibrary1 = GetMaterialLibrary()
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [900, 900]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.0] 
renderView1.CameraParallelScale = 7.2
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]

renderView1.OSPRayMaterialLibrary = materialLibrary1
 # init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = r'x'
renderView1.AxesGrid.YTitle = r'y'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.XTitleFontSize = 30
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.YTitleFontSize = 30
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.XLabelFontSize = 20
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.YLabelFontSize = 20
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontFile = ''
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [-4.0, -2.0, 2.0, 4.0]
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [-4.0, -2.0, 2.0, 4.0]
SetActiveView(renderView1)

# show data in view
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay = Show(waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.Representation = 'Surface'
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.Representation = 'Wireframe'
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.AmbientColor = [0.0, 0.0, 0.0]
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.Opacity = 0.4
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.LineWidth	 = 1.25
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.SetScalarBarVisibility(renderView1, False)


# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()



# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu)
calculator1.ResultArrayName = 'normeU'
calculator1.Function = 'sqrt(ExactVelocity_X*ExactVelocity_X + ExactVelocity_Y*ExactVelocity_Y)'
#calculator1.Function = 'sqrt("Velocity at cells results_Velocity at cells x_(m/s)"*"Velocity at cells results_Velocity at cells x_(m/s)" +"Velocity at cells results_Velocity at cells y_(m/s)"*"Velocity at cells results_Velocity at cells y_(m/s)")'

calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

Hide(waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu, renderView1)

calculator1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()
# Get min and maximal values
data = servermanager.Fetch(calculator1)
numCells = data.GetNumberOfCells()
#print "Number Of Cells: ", numCells
min,max = data.GetCellData().GetArray("normeU").GetRange()
#print("min : %.1e"%min)
#print("max : %.1e"%max)
min='{:0.1e}'.format(min)
max='{:0.1e}'.format(max)
# create a new 'Text'
text1 = Text()
text1.Text = 'Min: '+str(min)+'\n'+'Max:'+str(max)

# show data from text1
text1Display = Show(text1, renderView1)

# trace defaults for the display properties.
text1Display.Color = [0.0, 0.0, 0.0]
text1Display.FontFile = ''


# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(text1)
# ----------------------------------------------------------------


# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=calculator1)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')
# hide data in view
Hide(calculator1, renderView1)
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=cellDatatoPointData1)
contour1.Isosurfaces = [1.029388474872682, 0.26567167057115193, 0.4353865159714919, 0.6051013613718319, 0.7748162067721719, 0.9445310521725119, 1.114245897572852, 1.2839607429731918, 1.453675588373532, 1.6233904337738718, 1.7931052791742117]
# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
# trace defaults for the display properties.
contour1Display.Representation = 'Wireframe'
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]
contour1Display.LineWidth = 2.8
contour1Display.OSPRayScaleArray = 'MOMENT'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'Velocity'
contour1Display.ScaleFactor = 0.5158011198043824
contour1Display.SelectScaleArray = 'MOMENT'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'MOMENT'
contour1Display.GaussianRadius = 0.025790055990219118
contour1Display.SetScaleArray = ['POINTS', 'MOMENT']
contour1Display.SetScalarBarVisibility(renderView1, False)

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu)

# show data in view
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay = Show(waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtu, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
waveStaggered_2DCylinderDeflection_VelocityAtCells_35934vtuDisplay.SetScalarBarVisibility(renderView1, False)


WriteImage("ESSAI.png")
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
