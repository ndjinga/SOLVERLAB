from __future__ import division

from salome.geom import geomBuilder
from salome.smesh import smeshBuilder

import sys
import os

from math import pi, cos, sin

geompy = geomBuilder.New()
smesh = smeshBuilder.New()


def create_group_from(name, mother_shape, list_elem, type="EDGE"):
    new = geompy.CreateGroup(mother_shape, geompy.ShapeType[type])
    geompy.UnionList(new, list_elem)
    new.SetName(name)
    geompy.addToStudyInFather(mother_shape, new, name)
    return new


r = 1.
NumberOfSegments = 200

points = [geompy.MakeVertex(r * cos(i * pi / 3), r * sin(i * pi / 3), 0) for i in range(7)]
edges = [geompy.MakeEdge(points[i], points[i + 1]) for i in range(6)]
wire = geompy.MakeWire(edges)
hexa = geompy.MakeFace(wire, True)
geompy.addToStudy(hexa, "hexa")
g = create_group_from("boundaries", hexa, [geompy.GetInPlace(hexa, wire, 1)])

mesh = smesh.Mesh(hexa, "mesh")
msurf = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters_1 = msurf.Parameters(smeshBuilder.SIMPLE)
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments( NumberOfSegments )

mesh.Group(g)

mesh.Compute()
mesh.ExportMED("meshHexagonWithTriangles"+str(NumberOfSegments)+".med")
