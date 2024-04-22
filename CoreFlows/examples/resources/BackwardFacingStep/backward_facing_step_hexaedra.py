# Description: channel with backward facing step geometry, catesian mesh
#
# ---------------------------------------------------------------------------- #
#									       #
# Backward-facing step 1                                                       #
#                                                                              #
# References                                                                   #
# [1] D.C. Wan, B. S. V. Patnaik, and G. W. Wei.                               #
#     Discrete Singular Convolution- Finite Subdomain Method for the Solution  #
#     of Incompressible Viscous Flows.                                         #
#     J. Comp. Phys. 180, 229-255 (2002)                                       #
# [2] B. F. Armaly, F. Durst, J. C. F. Pereira, and B. Schoung.                #
#     Experimental and theoretical investigation of backward-facing step flow. #
#     J. Fluid Mech. 127, 473 (1983)                                           #
#                                                                              #
# Re = 25								       #
# Experimental reattachment length = 1.7550                                    #
#                                                                              #
# mesh: dx = 0.2, dy = 0.2						       #
# ---------------------------------------------------------------------------- #

import salome
from salome.geom import geomBuilder
geompy = geomBuilder.New()

### ------------------ Geometry
#
#     __________________________
#    |        |                 |
#    | Quad 1 |      Quad 3     |
#    |________|_________________|
#             |                 |
#             |      Quad 2     |
#             |_________________|
#
#
# Quad 1

g_q1_v1 = geompy.MakeVertex( -4.0, 1.0, 0.0 )
g_q1_v2 = geompy.MakeVertex(  0.0, 1.0, 0.0 )
g_q1_v3 = geompy.MakeVertex(  0.0, 2.0, 0.0 )
g_q1_v4 = geompy.MakeVertex( -4.0, 2.0, 0.0 )

g_q1_e1 = geompy.MakeEdge( g_q1_v1, g_q1_v2)
g_q1_e2 = geompy.MakeEdge( g_q1_v2, g_q1_v3)
g_q1_e3 = geompy.MakeEdge( g_q1_v3, g_q1_v4)
g_q1_e4 = geompy.MakeEdge( g_q1_v4, g_q1_v1)

g_q1 = geompy.MakeQuad( g_q1_e1, g_q1_e2, g_q1_e3, g_q1_e4 )


# Quad 2

g_q2_v1 = geompy.MakeVertex(  0.0, 0.0, 0.0 )
g_q2_v2 = geompy.MakeVertex( 18.0, 0.0, 0.0 )
g_q2_v3 = geompy.MakeVertex( 18.0, 1.0, 0.0 )
g_q2_v4 = geompy.MakeVertex(  0.0, 1.0, 0.0 )

g_q2_e1 = geompy.MakeEdge( g_q2_v1, g_q2_v2)
g_q2_e2 = geompy.MakeEdge( g_q2_v2, g_q2_v3)
g_q2_e3 = geompy.MakeEdge( g_q2_v3, g_q2_v4)
g_q2_e4 = geompy.MakeEdge( g_q2_v4, g_q2_v1)

g_q2 = geompy.MakeQuad( g_q2_e1, g_q2_e2, g_q2_e3, g_q2_e4 )

# Quad 3

g_q3_v1 = geompy.MakeVertex(  0.0, 1.0, 0.0 )
g_q3_v2 = geompy.MakeVertex( 18.0, 1.0, 0.0 )
g_q3_v3 = geompy.MakeVertex( 18.0, 2.0, 0.0 )
g_q3_v4 = geompy.MakeVertex(  0.0, 2.0, 0.0 )

g_q3_e1 = geompy.MakeEdge( g_q3_v1, g_q3_v2)
g_q3_e2 = geompy.MakeEdge( g_q3_v2, g_q3_v3)
g_q3_e3 = geompy.MakeEdge( g_q3_v3, g_q3_v4)
g_q3_e4 = geompy.MakeEdge( g_q3_v4, g_q3_v1)

g_q3 = geompy.MakeQuad( g_q3_e1, g_q3_e2, g_q3_e3, g_q3_e4 )

# Channel geometry

###g_for_mesh = geompy.MakeCompound( [g_q1, g_q2, g_qb] )
###g_for_mesh = geompy.MakeGlueFaces(g_for_mesh, 1.e-3)
g_for_mesh = geompy.MakeSewing( [g_q1, g_q2, g_q3], 1.e-3 )

# Groups

def add_edge_to_group(g, e):
    v = geompy.MakeCDG(e)
    p = geompy.GetMainShape(g)
    s = geompy.GetEdgeNearPoint(p, v)
    i = geompy.GetSubShapeID(p, s)
    geompy.AddObject(g, i)

group_inlet = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_inlet, g_q1_e4)

group_outlet = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_outlet, g_q2_e2)
add_edge_to_group(group_outlet, g_q3_e2)

group_wall = geompy.CreateGroup(g_for_mesh, geompy.ShapeType["EDGE"])
add_edge_to_group(group_wall, g_q1_e1)
add_edge_to_group(group_wall, g_q3_e3)
add_edge_to_group(group_wall, g_q2_e1)
add_edge_to_group(group_wall, g_q1_e3)
add_edge_to_group(group_wall, g_q2_e4)

##################### Meshing
import  SMESH
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

mesh = smesh.Mesh(g_for_mesh, "Mesh_backward_step")

nb_cells_nx = 112#Number of segment on the top horizontal wall
nb_cells_ny = 12 #Number of segments at outlet

# Hypotheses
hyp_inlet = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( -4.0, 1.5, 0.0 ) ) )
hyp_inlet.NumberOfSegments( int(nb_cells_ny / 2) )
hyp_inlet.Propagation()

hyp_q3_outlet = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 18.0, 1.5, 0.0 ) ) )
hyp_q3_outlet.NumberOfSegments( int(nb_cells_ny / 2) )
hyp_q3_outlet.Propagation()

hyp_q2_outlet = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 18.0, 0.5, 0.0 ) ) )
hyp_q2_outlet.NumberOfSegments( int(nb_cells_ny / 2) )
hyp_q2_outlet.Propagation()

hyp_q1_bottom_wall = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( -2.0, 1.0, 0.0 ) ) )
hyp_q1_bottom_wall.NumberOfSegments( int(nb_cells_nx / 5) )
hyp_q1_bottom_wall.Propagation()

hyp_q1_top_wall = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( -2.0, 2.0, 0.0 ) ) )
hyp_q1_top_wall.NumberOfSegments( int(nb_cells_nx / 5) )
hyp_q1_top_wall.Propagation()

hyp_q2_bottom_wall = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 9.0, 0.0, 0.0 ) ) )
hyp_q2_bottom_wall.NumberOfSegments( int(nb_cells_nx * 4/ 5) )
hyp_q2_bottom_wall.Propagation()

hyp_q2_left_wall = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 0.0, 0.5, 0.0 ) ) )
hyp_q2_left_wall.NumberOfSegments( int(nb_cells_ny / 2) )
hyp_q2_left_wall.Propagation()

hyp_q3_wall = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 9.0, 2.0, 0.0 ) ) )
hyp_q3_wall.NumberOfSegments( int(nb_cells_nx * 4/ 5) )
hyp_q3_wall.Propagation()


## Set algorithm and groups
mesh.Segment().AutomaticLength( 2./nb_cells_ny )
mesh.Quadrangle()
Inlet_1 = mesh.GroupOnGeom(group_inlet, 'Inlet',SMESH.EDGE)
Outlet_1 = mesh.GroupOnGeom(group_outlet,'Outlet',SMESH.EDGE)
Wall_1 = mesh.GroupOnGeom(group_wall,'Wall',SMESH.EDGE)

# Compute mesh
mesh.Compute()

# Export MED
mesh.ExportMED( "./backward_facing_step_hexaedra.med")
