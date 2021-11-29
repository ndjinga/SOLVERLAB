#!/usr/bin/env python
# -*- coding: UTF8 -*-

import sys

import medcoupling as MC
import MEDLoader as ML

if len(sys.argv) != 2:
  print( "USAGE: convert_2Dgmsh_to_med.py file.typ2")
  sys.exit(-1)

filename = sys.argv[1]

print( "Converting ", filename)

# type de maille en fonction du nombre de noeuds.
# cf INTERP_KERNEL/CellModel.cxx
d_cell_types = {3: MC.NORM_TRI3,
                4: MC.NORM_QUAD4,
                5: MC.NORM_POLYGON,
                6: MC.NORM_POLYGON}

mesh_dim = 2

read_vertices = False
coords = []

read_cells = False
read_new_cell_connectivity = False


vertex_line_number=False
cell_line_number=False

nb_vertices = 0
nb_cells = 0
nb_nodes_in_cell = 0
cell_connectivity = []

mesh = MC.MEDCouplingUMesh.New()
mesh.setMeshDimension(mesh_dim)
mesh.setName("mesh_from_typ")

with open(filename, 'r', encoding="iso-8859-1") as f:
  for line in f:
    # remove end of line character
    line = line[:-1]
    infos = line.split()
    #print(infos)
    if vertex_line_number:
      print("Reading number of vertices")
      nb_vertices = int(infos[0])
      vertex_line_number = False
      read_vertices = True
      print("Start reading vertex coordinates")
      continue
    elif cell_line_number:
      print("Reading number of cells")
      nb_cells = int(infos[0])
      mesh.allocateCells(nb_cells)
      cell_line_number = False
      read_cells = True
      print("Start reading cell connectivity")
      continue
    elif infos and infos[0] == "Vertices":
      vertex_line_number = True
    elif infos and ( infos[0] == "cells" or infos[0] == "Control"):
      cell_line_number = True
      read_vertices = False
    elif read_vertices:
      # read node coords
      coords_i = [float(v) for v in infos]
      coords += coords_i
    elif read_cells:
      values = [int(v) for v in infos]
      nb_nodes_in_cell = int(values[0])
      cell_connectivity = values[1:]
      #print( "nb_nodes_in_cell: ", nb_nodes_in_cell)
      #print "cell_connectivity: ", cell_connectivity
      # start numbering at 0
      cell_connectivity = [v-1 for v in cell_connectivity]
      mesh.insertNextCell(d_cell_types[nb_nodes_in_cell], nb_nodes_in_cell, cell_connectivity)
      nb_nodes_in_cell = 0
      cell_connectivity = []
      
  meshCoords = MC.DataArrayDouble.New()
  #pdb.set_trace()
  meshCoords.setValues(coords, nb_vertices, mesh_dim)
  mesh.setCoords(meshCoords)


# Merge les noeuds confondus (à faire avant le conformize2D)
arr, areNodesMerged, newNbOfNodes = mesh.mergeNodes(1e-10)

# Crée des polyèdres pour rendre conforme les mailles
mesh.convertAllToPoly()
mesh.conformize2D(1e-10)
mesh.unPolyze()

# Crée les éléments 1D pour pouvoir imposer les conditions aux limites
mesh_1d = mesh.computeSkin()

# Identifie les faces de chaque côté pour créer les groupes
tol = 1e-10

barycenters = mesh_1d.computeIsoBarycenterOfNodesPerCell()
ids_left = []
ids_right = []
ids_bottom = []
ids_top = []
for i, coord in enumerate(barycenters):
  x, y = coord
  if abs(x) < tol:
    ids_left.append(i)
  elif abs(x-1) < tol:
    ids_right.append(i)
  elif abs(y) < tol:
    ids_bottom.append(i)
  elif abs(y-1) < tol:
    ids_top.append(i)

arr_left = MC.DataArrayInt64(ids_left)
arr_right = MC.DataArrayInt64(ids_right)
arr_bottom = MC.DataArrayInt64(ids_bottom)
arr_top = MC.DataArrayInt64(ids_top)

arr_left.setName("Left")
arr_right.setName("Right")
arr_bottom.setName("Bottom")
arr_top.setName("Top")

# Trie les cellules par type conformément à la convention MED fichier
o2n = mesh.sortCellsInMEDFileFrmt()
o2n = mesh_1d.sortCellsInMEDFileFrmt()
meshMEDFile = ML.MEDFileUMesh.New()
# Ecrit le maillage 2D
meshMEDFile.setMeshAtLevel(0,mesh)
# Ecrit le maillage 1D
meshMEDFile.setMeshAtLevel(-1,mesh_1d)
# Ecrit les groupes
meshMEDFile.addGroup(-1, arr_left)
meshMEDFile.addGroup(-1, arr_right)
meshMEDFile.addGroup(-1, arr_bottom)
meshMEDFile.addGroup(-1, arr_top)
med_filename = filename.replace(".typ2", ".med")
meshMEDFile.write(med_filename,2) # 2 stands for write from scratch

print( "...done converting ", filename, " to ", med_filename)
