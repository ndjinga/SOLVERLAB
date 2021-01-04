from __future__ import division

from salome.geom import geomBuilder
from salome.smesh import smeshBuilder

import sys
import os

from math import pi, cos, sin, sqrt
import SMESH

geompy = geomBuilder.New()
smesh = smeshBuilder.New()


def create_group_from(name, mother_shape, list_elem, type="EDGE"):
    new = geompy.CreateGroup(mother_shape, geompy.ShapeType[type])
    geompy.UnionList(new, list_elem)
    new.SetName(name)
    geompy.addToStudyInFather(mother_shape, new, name)
    return new

# We build a regular tetrahedron with center of mass at (0,0,0) and 
r = 1.#distance between any two vertices
h = sqrt(6.) / 3. * r # height of the tetrahedron (distance between a vertex and the opposite plane)

points = [geompy.MakeVertex(r * cos(2 * i * pi / 3), r * sin(2 * i * pi / 3), -1./4 * h) for i in range(3)]
points.append(geompy.MakeVertex(0., 0., 3./4 * h))

faces_connectivity = [[0, 1, 2], [0, 1, 3], [1, 2, 3], [2, 3, 0]]

faces = []
for fc in faces_connectivity:
    fc2 = fc + [fc[0]]
    edges = [geompy.MakeEdge(points[fc2[i]], points[fc2[i + 1]]) for i in range(len(fc2) - 1)]
    faces.append(geompy.MakeFace(geompy.MakeWire(edges), True))

tetra_skin = geompy.MakeShell(faces)
tetra = geompy.MakeSolid(tetra_skin)
geompy.addToStudy(tetra, "Tetrahedron")
g = create_group_from("TetrahedronBoundary", tetra, [geompy.GetInPlace(tetra, tetra_skin, 1)], type="FACE")

### Mesh ###
number_of_segments = 50

TetrahedronMesh = smesh.Mesh(tetra, "Tetrahedron"+str(number_of_segments))
NETGEN_1D_2D_3D = TetrahedronMesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)

NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments(number_of_segments)
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()

TetrahedronBoundary_1 = TetrahedronMesh.GroupOnGeom(g, 'BoundaryFaces', SMESH.FACE)
TetrahedronBoundary_2 = TetrahedronMesh.GroupOnGeom(g, 'BoundaryNodes', SMESH.NODE)
isDone = TetrahedronMesh.Compute()

TetrahedronMesh.ExportMED("meshTetrahedron" + str(number_of_segments) + ".med")
