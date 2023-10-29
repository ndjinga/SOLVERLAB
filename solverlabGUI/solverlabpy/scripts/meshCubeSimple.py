#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
### 
"""
dir = "/home/catB/ym268439/Documents"
import sys
sys.path.insert(0,dir)
import Mesh_cube_simple as MM
MM.myMeshCube("toto",250)

import importlib
importlib.reload(Mesh_cube_simple)
"""
import sys
import os
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/catB/ym268439/Documents')

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()
model.end()

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

def myMeshCube(name,size=200):
    geompy = geomBuilder.New()
    Box_1 = geompy.MakeBoxDXDYDZ(size, size, size)
    geompy.addToStudy( Box_1, name )

    ###
    ### SMESH component
    ###

    import  SMESH, SALOMEDS
    from salome.smesh import smeshBuilder

    smesh = smeshBuilder.New()
    #smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                     # multiples meshes built in parallel, complex and numerous mesh edition (performance)

    Mesh_1 = smesh.Mesh(Box_1)
    NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
    NETGEN_3D_Parameters_1.SetMaxSize( 34.641 )
    NETGEN_3D_Parameters_1.SetMinSize( 0.34641 )
    NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
    NETGEN_3D_Parameters_1.SetOptimize( 1 )
    NETGEN_3D_Parameters_1.SetFineness( 4 )
    NETGEN_3D_Parameters_1.SetChordalError( -1 )
    NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
    NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
    NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
    NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
    isDone = Mesh_1.Compute()
    # print("isDone %s" % isDone)
    smesh.SetName(Mesh_1, name)
    ## Set names of Mesh objects
    smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
    smesh.SetName(Mesh_1.GetMesh(), name)
    smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
    try:
        
      tmp=os.path.expandvars(r'${HOME}/Documents/%s.med' % name)
      Mesh_1.ExportMED(tmp,auto_groups=0,version=41,overwrite=1,meshPart=None,autoDimension=1)
      print("ExportMED() ok file name: '%s'" % tmp)

    except:
      print("ExportMED() failed. Invalid file name: '%s'" % name)





    if salome.sg.hasDesktop():
      salome.sg.updateObjBrowser()
