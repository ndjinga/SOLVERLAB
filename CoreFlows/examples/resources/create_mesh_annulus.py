# -*-coding:utf-8 -*

#### Décomposition radiale d'un disque (r0=0) ou anneau (0<r0<r1)
### input : xcenter, ycenter, radius_min, radius_max, nr, ntheta
### output : annulusSpiderWeb.vtu, annulusSpiderWeb.med

import medcoupling as mc
import MEDLoader as ML
from math import pi


def drawPolarGrid(center_x, center_y, r0, r1, angle0, angle1, n_r, n_theta, name):
  """ Build a polar grid, centered at (center_x, center_y), with n_r subdivisions in the radial direction
  and n_theta subdivisions in the angular direction.
  The radial coordinates start at r0 and end at r1, and the angular coordinates start at angle0 and end at
  angle1.
  Angles should be given in degrees.
  You can use this script directly in the custom shape catalog.
  """
  if n_r <= 0 or n_theta <= 0:
    raise ValueError("Invalid parameter! Number of grid steps n_r and n_theta must be positive")
  if r0 >= r1:
    raise ValueError("Invalid parameter r0<r1 ! Start radius must be smaller than end radius")
  if angle0 >= angle1:
    raise ValueError("Invalid parameter angle0 < angle1 ! Start angle must be smaller than end angle")
  m = mc.MEDCouplingCMesh("spider_web")
  arr_r = mc.DataArrayDouble(n_r+1);     arr_r.iota()
  arr_t = mc.DataArrayDouble(n_theta+1); arr_t.iota()
  m.setCoordsAt(0, arr_r)
  m.setCoordsAt(1, arr_t)
  m = m.buildUnstructured()
  # Now the real job:
  coo = m.getCoords()
  dr, dtheta = (r1-r0) / float(n_r), (angle1-angle0) / float(n_theta)
  coo[:,0] = r0 + dr*coo[:,0]
  coo[:,1] = (angle0 + dtheta*coo[:,1]) * pi/180.0
  coo = coo.fromPolarToCart()
  m.setCoords(coo)
  oldNbOfNodes=m.getNumberOfNodes()
  arr, areNodesMerged, newNbOfNodes=m.mergeNodes(1e-10)
  
  print( "oldNbOfNodes=",oldNbOfNodes,"newNbOfNodes=",newNbOfNodes )
  print( "m.getNumberOfCells()=", m.getNumberOfCells() )
  
  m.checkConsistencyLight()
  
  # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
  mesh_1d = m.computeSkin()

  # Trie les cellules par type conformément à la convention MED fichier
  o2n = m.sortCellsInMEDFileFrmt()
  meshMEDFile=ML.MEDFileUMesh.New()
  # Ecrit le maillage 2D
  meshMEDFile.setMeshAtLevel(0,m)
  # Ecrit le maillage 1D
  meshMEDFile.setMeshAtLevel(-1,mesh_1d)
  # Ecrit les groupes
  arr_circle = mc.DataArrayInt(list(range(mesh_1d.getNumberOfCells())))
  arr_circle.setName("Circle")
  meshMEDFile.addGroup(-1, arr_circle)

  filename = name+"SpiderWeb"+str(n_r)+"x"+str(n_theta)+".med"
  """ # Write the result into a VTU file that can be read with ParaView
  m.writeVTK(name+"SpiderWeb"+str(n_r)+"x"+str(n_theta)+".vtu") """
  # Write the result into a MED file that can be read with Salomé
  meshMEDFile.write(filename,2) # 2 stands for write from scratch
  
  return m

if __name__ == "__main__":
  Rmax = 6.
  Rmin = 0.8
  drawPolarGrid(0., 0.,Rmin, Rmax, 0., 360., 1, 4,"Annulus")
  drawPolarGrid(0., 0.,Rmin, Rmax, 0., 360., 5, 8,"Annulus")
  drawPolarGrid(0., 0.,Rmin, Rmax, 0., 360., 3, 4,"Annulus")
  for i in range(5):
    n = 2**i
    nr = n*5
    ntheta = n*16
    drawPolarGrid(0., 0.,Rmin, Rmax, 0., 360., nr, ntheta,"Annulus")
 
 
