# $Id: bump_geom_salome.py,v 1.1 2010/07/27 09:44:43 kumbaro Exp $
# $Name: FLIVAP-1-1-1 $
#
# Description: channel with sinusoidal bump geometry
#

import math
import geompy

### ------------------ Geometry
#
#     _____________________________
#    |        |           |        |
#    | Quad 1 | Quad bump | Quad 2 |
#    |        | _________ |        |
#    |________|/         \|________|
#
#
# Bump geometry:
def bump(x):
    if ((x >= 1.0) and (x <= 3.0)):
        y = 0.1 * (1.0 - math.cos((x - 1.0) * math.pi))
    else:
        y = 0
    return y

# Quad 1

g_q1_v1 = geompy.MakeVertex( 0.0, 0.0, 0.0 )
g_q1_v2 = geompy.MakeVertex( 1.0, 0.0, 0.0 )
g_q1_v3 = geompy.MakeVertex( 1.0, 1.0, 0.0 )
g_q1_v4 = geompy.MakeVertex( 0.0, 1.0, 0.0 )

g_q1_e1 = geompy.MakeEdge( g_q1_v1, g_q1_v2)
g_q1_e2 = geompy.MakeEdge( g_q1_v2, g_q1_v3)
g_q1_e3 = geompy.MakeEdge( g_q1_v3, g_q1_v4)
g_q1_e4 = geompy.MakeEdge( g_q1_v4, g_q1_v1)

g_q1 = geompy.MakeQuad( g_q1_e1, g_q1_e2, g_q1_e3, g_q1_e4 )

geompy.addToStudy(g_q1, "Quad_1")

# Quad 2

g_q2_v1 = geompy.MakeVertex( 3.0, 0.0, 0.0 )
g_q2_v2 = geompy.MakeVertex( 4.0, 0.0, 0.0 )
g_q2_v3 = geompy.MakeVertex( 4.0, 1.0, 0.0 )
g_q2_v4 = geompy.MakeVertex( 3.0, 1.0, 0.0 )

g_q2_e1 = geompy.MakeEdge( g_q2_v1, g_q2_v2)
g_q2_e2 = geompy.MakeEdge( g_q2_v2, g_q2_v3)
g_q2_e3 = geompy.MakeEdge( g_q2_v3, g_q2_v4)
g_q2_e4 = geompy.MakeEdge( g_q2_v4, g_q2_v1)

g_q2 = geompy.MakeQuad( g_q2_e1, g_q2_e2, g_q2_e3, g_q2_e4 )

geompy.addToStudy(g_q2, "Quad_2")

# Quad bump

g_vertices_bump = []
nx = 20
dx = 2.0 / nx
for ix in range(nx+1):
    x = 1.0 + (ix * dx)
    y = bump(x)
    print x
    print y
    g_vertices_bump.append( geompy.MakeVertex(x,y,0.0) )

g_bump = geompy.MakeInterpol(g_vertices_bump)
geompy.addToStudy(g_bump, "g_bump")

g_qb_e1 = g_bump
g_qb_e3 = geompy.MakeEdge( g_q2_v4, g_q1_v3 )

g_qb = geompy.MakeQuad( g_qb_e1, g_q2_e4, g_qb_e3, g_q1_e2 )

geompy.addToStudy(g_qb, "Quad_bump")

# Channel geometry

###g_for_mesh = geompy.MakeCompound( [g_q1, g_q2, g_qb] )
###g_for_mesh = geompy.MakeGlueFaces(g_for_mesh, 1.e-3)
g_for_mesh = geompy.MakeSewing( [g_q1, g_q2, g_qb], 1.e-3 )
id_g = geompy.addToStudy(g_for_mesh, "Channel_bump")

# Groups

def add_edge_to_group(g, e):
    v = geompy.MakeCDG(e)
    p = geompy.GetMainShape(g)
    s = geompy.GetEdgeNearPoint(p, v)
    i = geompy.GetSubShapeID(p, s)
    geompy.AddObject(g, i)

group_inlet = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_inlet, g_q1_e4)
#id_inlet = geompy.GetSubShapeID(g_for_mesh, g_q1_e4)
#geompy.AddObject(group_inlet, id_inlet)
geompy.addToStudy(group_inlet, "Inlet")

group_outlet = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_outlet, g_q2_e2)
geompy.addToStudy(group_outlet, "Outlet")

group_wall = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_wall, g_q1_e1)
add_edge_to_group(group_wall, g_qb_e1)
add_edge_to_group(group_wall, g_q2_e1)
add_edge_to_group(group_wall, g_q1_e3)
add_edge_to_group(group_wall, g_qb_e3)
add_edge_to_group(group_wall, g_q2_e3)
geompy.addToStudy(group_wall, "Wall")
