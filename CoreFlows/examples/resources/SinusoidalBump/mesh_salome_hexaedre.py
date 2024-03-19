# Description: channel with sinusoidal bump. SALOME meshing
# Mesh: 40x10 hexaedres


from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()


### ------------------ Geometry

exec(open( "./bump_geom_salome.py").read() )

### ------------------ Meshing

mesh = smesh.Mesh(g_for_mesh, "Mesh_channel_bump")

nb_cells_nx = 40
nb_cells_ny = 10

# Hypotheses
hyp_y = mesh.Segment()
hyp_y.NumberOfSegments(nb_cells_ny)

hyp_q1 = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 0.5, 1.0, 0.0 ) ) )
hyp_q1.NumberOfSegments( int(nb_cells_nx / 4) )
hyp_q1.Propagation()

hyp_q2 = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 3.5, 1.0, 0.0 ) ) )
hyp_q2.NumberOfSegments( int(nb_cells_nx / 4) )
hyp_q2.Propagation()

hyp_qb = mesh.Segment( geompy.GetEdgeNearPoint(g_for_mesh, geompy.MakeVertex( 2.0, 1.0, 0.0 ) ) )
hyp_qb.NumberOfSegments( int(nb_cells_nx / 2) )
hyp_qb.Propagation()

# Groups
mesh.Group(group_inlet,  "Inlet")
mesh.Group(group_outlet, "Outlet")
mesh.Group(group_wall,   "Wall")

# Quadrangle cells
mesh.Quadrangle()
mesh.Compute()

# Export MED
mesh.ExportMED("./bump_40x10_hexaedre.med")
