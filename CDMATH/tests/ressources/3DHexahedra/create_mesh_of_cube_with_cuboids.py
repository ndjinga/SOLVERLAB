# -*-coding:utf-8 -*

#### Maillage de paraléllépipèdes rectangles
### input : xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz
### output : cubeWithCuboids.vtu, cubeWithCuboids.med

import medcoupling as mc
import MEDLoader as ML

## Create a 3D grid mesh of a cuboid with lengths x_end-x_start, y_end-y_start and z_end-z_start
# @param nx the number of segments on the x-axis of the grid
# @param ny the number of segments on the y-axis of the grid
# @param nz the number of segments on the z-axis of the grid
# @param xmin the x coordinate of the first point on the x-axis of the grid
# @param xmax the x coordinate of the last point on the x-axis of the grid
# @param ymin the y coordinate of the first point on the y-axis of the grid
# @param ymax the y coordinate of the last point on the y-axis of the grid
# @param zmin the z coordinate of the first point on the z-axis of the grid
# @param zmax the z coordinate of the last point on the z-axis of the grid
# @param mesh_name the name of the mesh
def create3DGrid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, mesh_name="Mesh_cube_with_cuboids"):
  mesh_dim = 3
  dx = (xmax-xmin)/nx
  dy = (ymax-ymin)/ny
  dz = (zmax-zmin)/nz
  mesh = mc.MEDCouplingIMesh(mesh_name, mesh_dim, [nx+1,ny+1,nz+1], [xmin, ymin, zmin], [dx,dy,dz])
  return mesh

def mesh_cube_with_cuboids(xmin,xmax,nx,ymin,ymax,ny,zmin, zmax, nz, mesh_name="cubeWithCuboids"):
    
    print( "Meshing a cube with cuboids nx=",nx,"ny=",ny,"nz=",nz )
    mesh=create3DGrid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, mesh_name="Mesh_cube_with_cuboids")
    myMesh = mesh.buildUnstructured()

    #--------------- Boundary groups -----------------#
    # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
    mesh_2d = myMesh.computeSkin()
    
    # Identifie les segments de chaque côté pour créer les groupes
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
        if abs(y-ymin) < tol:
          ids_left.append(i)
        elif abs(y-ymax) < tol:
          ids_right.append(i)
        elif abs(z-zmin) < tol :
          ids_bottom.append(i)
        elif abs(z-zmax) < tol:
          ids_top.append(i)
        elif abs(x-xmin) < tol :
          ids_back.append(i)
        elif abs(x-xmax) < tol:
          ids_front.append(i)
        else:
            raise ValueError("Pb with boundary construction : barycenter does not belong to any border group")
        
    arr_left = mc.DataArrayInt(ids_left)
    arr_right = mc.DataArrayInt(ids_right)
    arr_bottom = mc.DataArrayInt(ids_bottom)
    arr_top = mc.DataArrayInt(ids_top)
    arr_front = mc.DataArrayInt(ids_front)
    arr_back = mc.DataArrayInt(ids_back)
    
    arr_left.setName("Left")
    arr_right.setName("Right")
    arr_bottom.setName("Bottom")
    arr_top.setName("Top")
    arr_front.setName("Front")
    arr_back.setName("Back")
    
    # Trie les cellules par type conformément à la convention MED fichier
    o2n = myMesh.sortCellsInMEDFileFrmt()
    meshMEDFile=ML.MEDFileUMesh.New()
    # Ecrit le maillage 2D
    meshMEDFile.setMeshAtLevel(0,myMesh)
    # Ecrit le maillage 1D
    meshMEDFile.setMeshAtLevel(-1,mesh_2d)
    # Ecrit les groupes
    meshMEDFile.addGroup(-1, arr_left)
    meshMEDFile.addGroup(-1, arr_right)
    meshMEDFile.addGroup(-1, arr_bottom)
    meshMEDFile.addGroup(-1, arr_top)
    meshMEDFile.addGroup(-1, arr_back)
    meshMEDFile.addGroup(-1, arr_front)
    
    # Check that everything is coherent (will throw if not)
    myMesh.checkConsistencyLight()
    
    filename = mesh_name+".med"
    # Write the result into a VTU file that can be read with ParaView
    #mesh.writeVTK(mesh_name+".vtu")
    # Write the result into a MED file that can be read with Salomé
    meshMEDFile.write(filename,2) # 2 stands for write from scratch

if __name__ == '__main__':
    nx=5; ny=5; nz=5
    mesh_cube_with_cuboids(0.,1.,nx,0.,1.,ny,0.,1.,nz,"squareWithCrossTriangles")
    
