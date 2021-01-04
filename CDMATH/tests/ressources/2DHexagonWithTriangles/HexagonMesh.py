from __future__ import division

from salome.geom import geomBuilder
from salome.smesh import smeshBuilder

import sys
import os

from math import pi, cos, sin
import SMESH

geompy = geomBuilder.New()
smesh = smeshBuilder.New()


def create_group_from(name, mother_shape, list_elem, type="EDGE"):
    new = geompy.CreateGroup(mother_shape, geompy.ShapeType[type])
    geompy.UnionList(new, list_elem)
    new.SetName(name)
    geompy.addToStudyInFather(mother_shape, new, name)
    return new


r = 1.
NumberOfSegments = 10

points = [geompy.MakeVertex(r * cos(i * pi / 3), r * sin(i * pi / 3), 0) for i in range(7)]
edges = [geompy.MakeEdge(points[i], points[i + 1]) for i in range(6)]
wire = geompy.MakeWire(edges)
hexa = geompy.MakeFace(wire, True)
geompy.addToStudy(hexa, "Hexagon")
g = create_group_from("HexagonBoundary", hexa, [geompy.GetInPlace(hexa, wire, 1)])

mesh = smesh.Mesh(hexa, "HexagonWith"+str(NumberOfSegments)+"Triangles")
msurf = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters_1 = msurf.Parameters(smeshBuilder.SIMPLE)
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments( NumberOfSegments )

HexagonBoundary_1 = mesh.GroupOnGeom(g, 'BoundaryFaces', SMESH.EDGE)
HexagonBoundary_2 = mesh.GroupOnGeom(g, 'BoundaryNodes', SMESH.NODE)

mesh.Compute()
mesh.ExportMED("meshHexagonWithTriangles"+str(NumberOfSegments)+".med")
