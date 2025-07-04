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
""" Vertex_1 = geompy.MakeVertex(-2, -2, 0)
Vertex_2 = geompy.MakeVertex(-2, 3, 0)
Vertex_3 = geompy.MakeVertex(3, 3, 0)
Vertex_4 = geompy.MakeVertex(3, -2, 0) """
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(0, 1, 0)
Vertex_3 = geompy.MakeVertex(1, 1, 0)
Vertex_4 = geompy.MakeVertex(1, 0, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_1)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Face_1, 'Face_1' )

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
Gmsh_Parameters.SetMinSize( 0.001 )
Gmsh_Parameters.SetMaxSize( 0.05 )
Gmsh_Parameters.SetIs2d( 1 )
GMSH_2D = smesh.CreateHypothesis('GMSH_2D', 'GMSHEngine')
Quadrangle_Parameters_1 = smesh.CreateHypothesis('QuadrangleParams')
Quadrangle_2D = smesh.CreateHypothesis('Quadrangle_2D')
Quadrangle_Parameters_1.SetQuadType( smeshBuilder.QUAD_STANDARD )
Quadrangle_Parameters_1.SetTriaVertex( 4 )
Quadrangle_Parameters_1.SetCorners( [ 4, 5, 7, 9 ] )
Quadrangle_Parameters_1.SetEnforcedNodes( [], [] )
Quadrangle_Parameters_2 = smesh.CreateHypothesis('QuadrangleParams')
Quadrangle_Parameters_2.SetQuadType( smeshBuilder.QUAD_QUADRANGLE_PREF_REVERSED )
Quadrangle_Parameters_2.SetTriaVertex( -1 )
Quadrangle_Parameters_2.SetEnforcedNodes( [], [] )
Mesh_1 = smesh.Mesh(Face_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(200)
status = Mesh_1.AddHypothesis(Quadrangle_Parameters_2)
status = Mesh_1.AddHypothesis(Quadrangle_2D)
Number_of_Segments_1.SetNumberOfSegments( 55 )
isDone = Mesh_1.Compute()
isDone = Mesh_1.SplitQuadObject( Mesh_1, 1 )


## Set names of Mesh objects
smesh.SetName(GMSH_2D, 'GMSH_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D, 'Quadrangle_2D')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Quadrangle_Parameters_2, 'Quadrangle Parameters_2')
Mesh_1.ExportMED( r'/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/SquareWithTriangles'+str(55)+'.med', auto_groups=1, minor=40)


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
