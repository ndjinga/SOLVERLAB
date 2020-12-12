# -*-coding:utf-8 -*

import medcoupling as mc
import MEDLoader as ML

## Create a 2D grid mesh of a rectangle with lengths x_end-x_start and y_end-y_start
# @param nx the number of segments on the x-axis of the grid
# @param ny the number of segments on the y-axis of the grid
# @param xmin the x coordinate of the first point on the x-axis of the grid
# @param xmax the x coordinate of the last point on the x-axis of the grid
# @param ymin the y coordinate of the first point on the y-axis of the grid
# @param ymax the y coordinate of the last point on the y-axis of the grid
# @param mesh_name the name of the mesh
def mesh_rectangle_with_rectangles(xmin, xmax, nx, ymin, ymax, ny, mesh_name="Mesh_rectangle_with_rectangles"):
    mesh_dim = 2
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny
    mesh = mc.MEDCouplingIMesh(mesh_name, mesh_dim, [nx+1, ny+1], [xmin, ymin], [dx,dy]).buildUnstructured()

    #--------------- Boundary groups -----------------#
    # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
    mesh_1d = mesh.computeSkin()
    
    # Identifie les segments de chaque côté pour créer les groupes
    tol = 1e-10
    
    barycenters = mesh_1d.computeIsoBarycenterOfNodesPerCell()
    ids_left = []
    ids_right = []
    ids_bottom = []
    ids_top = []
    for i, coord in enumerate(barycenters):
        x, y = coord
        if abs(x-xmin) < tol:
          ids_left.append(i)
        elif abs(x-xmax) < tol:
          ids_right.append(i)
        elif abs(y-ymin) < tol :
          ids_bottom.append(i)
        elif abs(y-ymax) < tol:
          ids_top.append(i)
        else:
            raise ValueError("Pb with boundary construction : barycenter does not belong to any border group")
        
    arr_left = mc.DataArrayIdType(ids_left)
    arr_right = mc.DataArrayIdType(ids_right)
    arr_bottom = mc.DataArrayIdType(ids_bottom)
    arr_top = mc.DataArrayIdType(ids_top)
    
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
    
    # Check that everything is coherent (will throw if not)
    mesh.checkConsistencyLight()
    
    filename = mesh_name+".med"
    # Write the result into a VTU file that can be read with ParaView
    #mesh.writeVTK(mesh_name+".vtu")
    # Write the result into a MED file that can be read with Salomé
    meshMEDFile.write(filename,2) # 2 stands for write from scratch

if __name__ == '__main__':
    nx=5
    mesh_rectangle_with_rectangles(0.,1.,nx,0.,1.,nx,"squareWithSquares")
    
