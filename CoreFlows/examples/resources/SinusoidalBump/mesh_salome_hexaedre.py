# $Id: mesh_salome_hexaedre.py,v 1.1 2010/07/27 09:44:33 kumbaro Exp $
# $Name: FLIVAP-1-1-1 $
#
# Description: channel with sinusoidal bump. SALOME meshing
# Mesh: 40x10 hexaedres

import os
import smesh

### ------------------ Geometry

dipha_root = os.getenv("DIPHA_ROOT")
dir_test   = dipha_root + "/Tests/3eqs/channel_sinus_bump/"
execfile(dir_test + "Meshing/bump_geom_salome.py")

### ------------------ Meshing

mesh = smesh.Mesh(g_for_mesh, "Mesh_channel_bump")

nb_cells_nx = 40
nb_cells_ny = 10

# Hypotheses
hyp_y = mesh.Segment()
hyp_y.NumberOfSegments(nb_cells_ny)

hyp_q1 = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 0.5, 1.0, 0.0 ) ) )
hyp_q1.NumberOfSegments(nb_cells_nx / 4)
hyp_q1.Propagation()

hyp_q2 = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 3.5, 1.0, 0.0 ) ) )
hyp_q2.NumberOfSegments(nb_cells_nx / 4)
hyp_q2.Propagation()

hyp_qb = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 2.0, 1.0, 0.0 ) ) )
hyp_qb.NumberOfSegments(nb_cells_nx / 2)
hyp_qb.Propagation()

# Groups
mesh.Group(group_inlet,  "Inlet")
mesh.Group(group_outlet, "Outlet")
mesh.Group(group_wall,   "Wall")

# Quadrangle cells
mesh.Quadrangle()
mesh.Compute()

# Export MED

mesh.ExportMED(dir_test + "Meshing/Mesh_40x10/bump_40x10_hexaedre.med")
