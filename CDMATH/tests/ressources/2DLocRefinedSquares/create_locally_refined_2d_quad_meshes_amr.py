# -*- coding: UTF8 -*-

import medcoupling as MC
import MEDLoader as ML

## Create a 2D uniform grid mesh of a square with length x_end-x_start
# @param nb_segs_x the number of segments on the x-axis of the grid
# @param x_start the x coordinate of the first point on the x-axis of the grid
# @param x_end the x coordinate of the last point on the x-axis of the grid
# @param y_start the y coordinate of the first point on the y-axis of the grid
# @param mesh_name the name of the mesh
def createMesh(nb_segs_x, x_start=0, x_end=1., y_start=0., mesh_name="Mesh"):
  mesh_dim = 2
  nb_nodes = nb_segs_x+1
  dx = (x_end-x_start)/nb_segs_x
  mesh = MC.MEDCouplingIMesh(mesh_name, mesh_dim, [nb_nodes, nb_nodes], [x_start, y_start], [dx]*mesh_dim)
  return mesh

def createLocallyRefinedMesh(nb_segs_x, mesh_name):
  
  # First mesh
  mesh_1 = createMesh(nb_segs_x, 0., 1., 0., mesh_name)
  
  amr = MC.MEDCouplingCartesianAMRMesh(mesh_1)
  
  # 1er raffinement
  amr.addPatch([(nb_segs_x/2,nb_segs_x),(0,nb_segs_x/2)],[2,2])
  # 2eme raffinement
  amr[0].addPatch([(nb_segs_x/2,nb_segs_x),(0,nb_segs_x/2)],[2,2])
  
  # Crée un seul maillage avec tous les rafinements
  mesh = amr.buildUnstructured()
  mesh.setName(mesh_name)
  # Merge les noeuds confondus (à faire avant le conformize2D)
  arr, areNodesMerged, newNbOfNodes = mesh.mergeNodes(1e-10)
  # Crée des polygones pour rendre conforme les mailles
  mesh.conformize2D(1e-10)

  # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
  mesh_1d = mesh.computeSkin()

  # Identifie les segments de chaque côté pour créer les groupes
  tol = 1e-10
  #tol2 = 0
  #arr_left = mesh_1d.getCellsInBoundingBox([0-tol, tol, -tol, 1+tol], tol2)
  #arr_right = mesh_1d.getCellsInBoundingBox([1-tol, 1+tol, -tol, 1+tol], tol2)
  #arr_bottom = mesh_1d.getCellsInBoundingBox([0-tol, 1+tol, -tol, tol], tol2)
  #arr_top = mesh_1d.getCellsInBoundingBox([0-tol, 1+tol, 1-tol, 1+tol], tol2)
  
  # PB: getCellsInBoundingBox renvoie aussi les segments qui touchent la bounding box
  # => On boucle sur les coordonnées des barycentres

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

  arr_left = MC.DataArrayInt(ids_left)
  arr_right = MC.DataArrayInt(ids_right)
  arr_bottom = MC.DataArrayInt(ids_bottom)
  arr_top = MC.DataArrayInt(ids_top)
  
  arr_left.setName("Left")
  arr_right.setName("Right")
  arr_bottom.setName("Bottom")
  arr_top.setName("Top")

  # Trie les cellules par type conformément à la convention MED fichier
  o2n = mesh.sortCellsInMEDFileFrmt()
  meshMEDFile=ML.MEDFileUMesh.New()
  # Ecrit le maillage 2D
  meshMEDFile.setMeshAtLevel(0,mesh)
  # Ecrit le maillage 1D
  meshMEDFile.setMeshAtLevel(-1,mesh_1d)
  # Ecrit les groupes
  meshMEDFile.addGroup(-1, arr_left)
  meshMEDFile.addGroup(-1, arr_right)
  meshMEDFile.addGroup(-1, arr_bottom)
  meshMEDFile.addGroup(-1, arr_top)
  filename = mesh_name+".med"
  meshMEDFile.write(filename,2) # 2 stands for write from scratch
  
  return meshMEDFile

if __name__ == '__main__':
  createLocallyRefinedMesh(4, "squareWithLocRefSquares_1")
  createLocallyRefinedMesh(8, "squareWithLocRefSquares_2")
  createLocallyRefinedMesh(16, "squareWithLocRefSquares_3")
  createLocallyRefinedMesh(32, "squareWithLocRefSquares_4")
  createLocallyRefinedMesh(64, "squareWithLocRefSquares_5")
  createLocallyRefinedMesh(128, "squareWithLocRefSquares_6")
  createLocallyRefinedMesh(256, "squareWithLocRefSquares_7")
