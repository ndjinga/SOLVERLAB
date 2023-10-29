#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math
import salomepy.xsalomesession as XSS

def myRun():
  geompy = XSS.getGeompy()

  O = geompy.MakeVertex(0, 0, 0)
  OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
  OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
  OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
  Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
  O_1 = geompy.MakeVertex(0, 0, 0)
  OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
  OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
  OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
  geompy.addToStudy( O, 'O' )
  geompy.addToStudy( OX, 'OX' )
  geompy.addToStudy( OY, 'OY' )
  geompy.addToStudy( OZ, 'OZ' )
  geompy.addToStudy( Box_1, 'Box_1' )
  geompy.addToStudy( O_1, 'O' )
  geompy.addToStudy( OX_1, 'OX' )
  geompy.addToStudy( OY_1, 'OY' )
  geompy.addToStudy( OZ_1, 'OZ' )

  smeshBuilder = XSS.getSmeshBuilder()
  smesh = XSS.getSmesh()
  
  Mesh_1 = smesh.Mesh(Box_1)
  NETGEN_2D3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
  NETGEN_3D_Simple_Parameters_1 = smesh.CreateHypothesis('NETGEN_SimpleParameters_3D', 'NETGENEngine')
  NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( 4 )
  NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
  NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
  NETGEN_3D_Simple_Parameters_2 = NETGEN_2D3D.Parameters(smeshBuilder.SIMPLE)
  NETGEN_3D_Simple_Parameters_2.SetNumberOfSegments( 3 )
  NETGEN_3D_Simple_Parameters_2.LengthFromEdges()
  NETGEN_3D_Simple_Parameters_2.LengthFromFaces()
  isDone = Mesh_1.Compute()


  ## Set names of Mesh objects
  smesh.SetName(NETGEN_2D3D.GetAlgorithm(), 'NETGEN_2D3D')
  smesh.SetName(NETGEN_3D_Simple_Parameters_2, 'NETGEN 3D Simple Parameters_2')
  smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
  smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

  XSS.updateObjBrowser()

if __name__ == '__main__':
  myRun()
