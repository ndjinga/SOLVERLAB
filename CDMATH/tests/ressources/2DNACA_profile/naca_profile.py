# -*- coding: utf-8 -*-
from __future__ import print_function, division
from math import sqrt
from salome.geom  import geomBuilder
from salome.smesh import smeshBuilder

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

def create_group_near(name, mother_shape, pos, list_grp, type, add=True):
    shapes = [geompy.GetShapesNearPoint(mother_shape, P, geompy.ShapeType[type]) for P in pos]
    new = create_group_from(name, mother_shape, shapes, type)
    if add:
        list_grp.append(new)
    return new, list_grp

def create_group_from(name, mother_shape, list_elem, type="EDGE"):
    new = geompy.CreateGroup(mother_shape, geompy.ShapeType[type])
    geompy.UnionList(new, list_elem)
    new.SetName(name)
    geompy.addToStudyInFather(mother_shape, new, name)
    return new

def yt(t, x, c, finite):
    a4 = 0.1015 if finite else 0.1036
    return t / 0.2 * (0.2969 * sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c)**2 + 0.2843 * (x / c)**3 - a4 * (x / c) ** 4)

# ------------------------------------------- #
#  parametres geometriques

c = 1.   # corde
t = 0.12  # episseur max (fraction de corde)
nx = 1000  # discretisation de la forme

# taille du domaine
hx, hy = 10.0 * c, 5.0 * c

#  parametres du maillage
NumberOfSegments = 10
dx_mesh = c / NumberOfSegments
angle_mesh = 8.
# ------------------------------------------- #

dx = c / nx
xC = [i * dx for i in range(nx + 1)]
yU = [ yt(t, x, c, False) for x in xC]
yD = [-yt(t, x, c, False) for x in xC]

profile  = [geompy.MakeVertex(x, y, 0) for x, y in zip(reversed(xC), reversed(yU))]
profile += [geompy.MakeVertex(x, y, 0) for x, y in zip(xC[1:], yD[1:])]

profile_boundary = geompy.MakeWire([geompy.MakeInterpol(profile)])
profile = geompy.MakeFace(profile_boundary, 1)
geompy.addToStudy(profile, "Profile")

domain = geompy.MakeFaceHW(hx, hy, 1)
domain = geompy.MakeCut(domain, profile)

geompy.addToStudy(domain, "Domain")
group = create_group_from("Airfoil_boundary", domain, [geompy.GetInPlace(domain, profile_boundary, 1)])
groups = [group]
new, groups = create_group_near("Inlet" , domain, [geompy.MakeVertex(-hx / 2., 0., 0.)], groups, "EDGE")
new, groups = create_group_near("Outlet", domain, [geompy.MakeVertex( hx / 2., 0., 0.)], groups, "EDGE")
new, groups = create_group_near("Top"   , domain, [geompy.MakeVertex(0.,  hy / 2., 0.)], groups, "EDGE")
new, groups = create_group_near("Bottom", domain, [geompy.MakeVertex(0., -hy / 2., 0.)], groups, "EDGE")

# maillage
mesh = smesh.Mesh(domain, "MeshFlow")
for g in groups: mesh.Group(g)
msurf = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters_1 = msurf.Parameters(smeshBuilder.SIMPLE)
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments( NumberOfSegments )
mesh.Compute()

mesh.ExportMED("naca"+str(NumberOfSegments)+".med")
