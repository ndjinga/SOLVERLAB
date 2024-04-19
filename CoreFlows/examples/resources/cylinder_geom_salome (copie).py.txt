#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.12.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources')

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

model.end()

###
### SHAPERSTUDY component
###

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Disk_1 = geompy.MakeDiskR(5.5, 1)
Disk_2 = geompy.MakeDiskR(1.2, 1)
Cut_1 = geompy.MakeCutList(Disk_1, [Disk_2], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Disk_2, 'Disk_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

QuadFromMedialAxis_1D2D = smesh.CreateHypothesis('QuadFromMedialAxis_1D2D')
RadialQuadrangle_1D2D = smesh.CreateHypothesis('RadialQuadrangle_1D2D')
Number_of_Segments_1 = smesh.CreateHypothesis('NumberOfSegments')
Number_of_Segments_1.SetNumberOfSegments( 160 )
Regular_1D = smesh.CreateHypothesis('Regular_1D')
Number_of_Layers_1 = smesh.CreateHypothesis('NumberOfLayers2D')
Number_of_Layers_1.SetNumberOfLayers( 50 )
Mesh_1 = smesh.Mesh(Cut_1,'Mesh_1')
status = Mesh_1.AddHypothesis(Number_of_Layers_1)
status = Mesh_1.AddHypothesis(QuadFromMedialAxis_1D2D)
isDone = Mesh_1.Compute()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/cylinder_geom_salome.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(RadialQuadrangle_1D2D, 'RadialQuadrangle_1D2D')
smesh.SetName(QuadFromMedialAxis_1D2D, 'QuadFromMedialAxis_1D2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Number_of_Layers_1, 'Number of Layers_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
