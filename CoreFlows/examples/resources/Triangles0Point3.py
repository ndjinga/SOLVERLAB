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
Disk_1 = geompy.MakeDiskR(6, 1)
Disk_2 = geompy.MakeDiskR(0.8, 1)
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

Gmsh_Parameters = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 0 )
Gmsh_Parameters.SetMaxSize( 1e+22 )
Gmsh_Parameters.SetSizeFactor( 300 )
Gmsh_Parameters.SetIs2d( 1 )
GMSH_2D = smesh.CreateHypothesis('GMSH_2D', 'GMSHEngine')
Gmsh_Parameters.SetSizeFactor( 1 )
Gmsh_Parameters.SetMaxSize( 0.2 )
Gmsh_Parameters.SetIs2d( 1 )
Gmsh_Parameters_1 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_1.Set2DAlgo( 0 )
Gmsh_Parameters_1.SetMinSize( 0 )
Gmsh_Parameters_1.SetMaxSize( 0.01 )
Gmsh_Parameters_1.SetIs2d( 1 )
Gmsh_Parameters_2 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_2.Set2DAlgo( 0 )
Gmsh_Parameters_2.SetMinSize( 0 )
Gmsh_Parameters_2.SetMaxSize( 0.1 )
Gmsh_Parameters_2.SetIs2d( 1 )
Gmsh_Parameters_3 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_3.Set2DAlgo( 0 )
Gmsh_Parameters_3.SetMinSize( 0 )
Gmsh_Parameters_3.SetMaxSize( 1 )
Gmsh_Parameters_3.SetIs2d( 1 )
Gmsh_Parameters_4 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_4.Set2DAlgo( 0 )
Gmsh_Parameters_4.SetMinSize( 0 )
Gmsh_Parameters_4.SetMaxSize( 0.5 )
Gmsh_Parameters_4.SetIs2d( 1 )
GmshMaxSize0Point3 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
GmshMaxSize0Point3.Set2DAlgo( 0 )
GmshMaxSize0Point3.SetMinSize( 0 )
GmshMaxSize0Point3.SetMaxSize( 0.3 )
GmshMaxSize0Point3.SetIs2d( 1 )
Mesh_1 = smesh.Mesh(Cut_1,'Mesh_1')
status = Mesh_1.AddHypothesis(GmshMaxSize0Point3)
status = Mesh_1.AddHypothesis(GMSH_2D)
isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(GMSH_2D, 'GMSH_2D')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Gmsh_Parameters_1, 'Gmsh Parameters')
smesh.SetName(Gmsh_Parameters_4, 'Gmsh Parameters')
smesh.SetName(GmshMaxSize0Point3, 'GmshMaxSize0Point3')
smesh.SetName(Gmsh_Parameters_2, 'Gmsh Parameters')
smesh.SetName(Gmsh_Parameters_3, 'Gmsh Parameters')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
