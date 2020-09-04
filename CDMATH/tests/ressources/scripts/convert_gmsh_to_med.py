#!/usr/bin/env python
# -*- coding: UTF8 -*-

import sys

import MEDCoupling as MC
import MEDLoader as ML

#filename = "locrafgrid_1_new.msh"
#filename = "checkerboard_2x2x2_new.msh"

if len(sys.argv) != 2:
  print "USAGE: convert_gmsh_to_med.py file.msh"
  sys.exit(-1)

filename = sys.argv[1]

print "Converting ", filename

# type de maille en fonction du nombre de noeuds.
# cf INTERP_KERNEL/CellModel.cxx
d_cell_types = {8: MC.NORM_HEXA8,
                6: MC.NORM_PENTA6}

mesh_dim = 3

read_vertices = False
coords = []

read_cells = False
read_new_cell_connectivity = False

nb_vertices = 0
nb_cells = 0
nb_nodes_in_cell = 0
cell_connectivity = []

mesh = MC.MEDCouplingUMesh.New()
mesh.setMeshDimension(mesh_dim)
mesh.setName("mesh_from_gmsh")

with open(filename, 'rb') as f:
  for line in f:
    # remove end of line character
    line = line[:-1]
    infos = line.split()
    #print infos
    if infos and infos[0] == "Vertices":
      nb_vertices = int(infos[1])
      read_vertices = True
    elif infos and infos[0] == "Volumes->faces":
      # stop reading node coords
      read_vertices = False
      meshCoords = MC.DataArrayDouble.New()
      #pdb.set_trace()
      meshCoords.setValues(coords, nb_vertices, mesh_dim)
      mesh.setCoords(meshCoords)
    elif read_vertices:
      # read node coords
      coords_i = [float(v) for v in infos]
      coords += coords_i
    elif infos and (infos[0] == "Volumes->Vertices" or infos[0] == "Volumes->Verticess"):
      # start reading cells connectivity
      nb_cells = int(infos[1])
      mesh.allocateCells(nb_cells)
      read_cells = True
      read_new_cell_connectivity = True
    elif infos and (infos[0] == "Faces->Edges" or infos[0] == "Faces->Edgess"):
      # stop reading cells connectivity
      read_cells = False
      read_faces = True
      nb_faces = int(infos[1])
    elif read_cells:
      values = [int(v) for v in infos]
      if read_new_cell_connectivity:
        #print ""
        #print "start new cell connectivity"
        nb_nodes_in_cell = int(values[0])
        cell_connectivity = values[1:]
        #print "nb_nodes_in_cell: ", nb_nodes_in_cell
        #print "cell_connectivity: ", cell_connectivity
        if len(cell_connectivity) < nb_nodes_in_cell:
          read_new_cell_connectivity = False
        else:
          read_new_cell_connectivity = True
      else:
          #print "complete the cell connectivity"
          cell_connectivity += values
          #print "cell_connectivity: ", cell_connectivity
          if len(cell_connectivity) == nb_nodes_in_cell:
            read_new_cell_connectivity = True
          else:
            read_new_cell_connectivity = False
      if read_new_cell_connectivity and cell_connectivity:
        #print "finish cell connectivity"
        #print nb_nodes_in_cell, cell_connectivity
        # start numbering at 0
        cell_connectivity = [v-1 for v in cell_connectivity]
        mesh.insertNextCell(d_cell_types[nb_nodes_in_cell], nb_nodes_in_cell, cell_connectivity)
        nb_nodes_in_cell = 0
        cell_connectivity = []

# Merge les noeuds confondus (à faire avant le conformize2D)
arr, areNodesMerged, newNbOfNodes = mesh.mergeNodes(1e-10)

# Crée des polyèdres pour rendre conforme les mailles
mesh.convertAllToPoly()
mesh.conformize3D(1e-10)
mesh.unPolyze()

# Crée les éléments 1D pour pouvoir imposer les conditions aux limites
mesh_2d = mesh.computeSkin()

# Identifie les faces de chaque côté pour créer les groupes
tol = 1e-10

barycenters = mesh_2d.computeIsoBarycenterOfNodesPerCell()
ids_left = []
ids_right = []
ids_bottom = []
ids_top = []
ids_front = []
ids_back = []
for i, coord in enumerate(barycenters):
  x, y, z = coord
  if abs(x) < tol:
    ids_left.append(i)
  elif abs(x-1) < tol:
    ids_right.append(i)
  elif abs(y) < tol:
    ids_bottom.append(i)
  elif abs(y-1) < tol:
    ids_top.append(i)
  elif abs(z) < tol:
    ids_back.append(i)
  elif abs(z-1) < tol:
    ids_front.append(i)

arr_left = MC.DataArrayInt(ids_left)
arr_right = MC.DataArrayInt(ids_right)
arr_bottom = MC.DataArrayInt(ids_bottom)
arr_top = MC.DataArrayInt(ids_top)
arr_back = MC.DataArrayInt(ids_back)
arr_front = MC.DataArrayInt(ids_front)

arr_left.setName("Left")
arr_right.setName("Right")
arr_bottom.setName("Bottom")
arr_top.setName("Top")
arr_back.setName("Back")
arr_front.setName("Front")

# Trie les cellules par type conformément à la convention MED fichier
o2n = mesh.sortCellsInMEDFileFrmt()
o2n = mesh_2d.sortCellsInMEDFileFrmt()
meshMEDFile = ML.MEDFileUMesh.New()
# Ecrit le maillage 3D
meshMEDFile.setMeshAtLevel(0,mesh)
# Ecrit le maillage 2D
meshMEDFile.setMeshAtLevel(-1,mesh_2d)
# Ecrit les groupes
meshMEDFile.addGroup(-1, arr_left)
meshMEDFile.addGroup(-1, arr_right)
meshMEDFile.addGroup(-1, arr_bottom)
meshMEDFile.addGroup(-1, arr_top)
meshMEDFile.addGroup(-1, arr_back)
meshMEDFile.addGroup(-1, arr_front)
med_filename = filename.replace(".msh", ".med")
meshMEDFile.write(med_filename,2) # 2 stands for write from scratch

print "...done"
