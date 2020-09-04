# -*-coding:utf-8 -*

#### Pavage d'un rectangle avec de triangles rectangles
### input : xmin, xmax, ymin, ymax, ny
### output : squareWithRightTriangles.med

import medcoupling as mc
import math
import MEDLoader as ML

def mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=10,ny=10,mesh_name="squareWithRightTriangles"):
    dx = (xmax-xmin)/nx 
    dy = (ymax-ymin)/ny 

    print( "Meshing a rectangle with Right Triangles nx=",nx,"ny=",ny, "ncells=",2*nx*ny)

    # Building the coordinates of the first initial triangle
    d_up = mc.DataArrayDouble(3,2)
    d_up[0,0] = 0
    d_up[0,1] = 0
    d_up[1,0] = 0
    d_up[1,1] = dy
    d_up[2,0] = dx
    d_up[2,1] = 0
    d_up.setInfoOnComponents(["X [m]","Y [m]"])
    
    # Building the coordinates of the second initial triangle
    d_down = mc.DataArrayDouble(3,2)
    d_down[0,0] = dx
    d_down[0,1] = dy
    d_down[1,0] = 0
    d_down[1,1] = dy
    d_down[2,0] = dx
    d_down[2,1] = 0
    d_down.setInfoOnComponents(["X [m]","Y [m]"])

    # translations of the cells : 
    translationToPerform   = [[xmin + i*dx, ymin + j*dy] for i in range(nx) for j in range(ny)]
        
    ds = (2*len(translationToPerform))*[None]
    for pos,t in enumerate(translationToPerform):
                     ds[2*pos]    = d_up[:]        # Perform a deep copy of d_up and place it at position 'pos' in ds
                     ds[2*pos]   += t           # Adding a vector to a set of coordinates does a translation
                     ds[2*pos+1]  = d_down[:]        # Perform a deep copy of d_down and place it at position 'pos' in ds
                     ds[2*pos+1] += t           # Adding a vector to a set of coordinates does a translation
                     pass
    
    d2 = mc.DataArrayDouble.Aggregate(ds)
    # Build an unstructured mesh representing the final pattern
    mesh = mc.MEDCouplingUMesh(mesh_name,2)
    mesh.setCoords(d2)
    print( "Mesh dimension is", mesh.getMeshDimension() )
    print( "Spatial dimension is", mesh.getCoords().getNumberOfComponents() )
    mesh.allocateCells(2*nx*ny)
    for i in range(2*nx*ny):
            cell_connec = [3*i,3*i+1,3*i+2]
            mesh.insertNextCell(mc.NORM_TRI3, cell_connec)
            pass
    
    # Identifying duplicate nodes
    oldNbOfNodes=mesh.getNumberOfNodes()        
    arr, areNodesMerged, newNbOfNodes=mesh.mergeNodes(1e-10)
    print( "oldNbOfNodes=",oldNbOfNodes,"newNbOfNodes",newNbOfNodes )
    
    # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
    mesh_1d = mesh.computeSkin()
    
    # Identifie les segments de chaque côté pour créer les groupes
    tol = 1e-10
    
    # PB: getCellsInBoundingBox renvoie aussi les segments qui touchent la bounding box
    # => On boucle sur les coordonnées des barycentres
    
    barycenters = mesh_1d.computeIsoBarycenterOfNodesPerCell()
    ids_left = []
    ids_right = []
    ids_bottom = []
    ids_top = []
    #print(barycenters)
    for i, coord in enumerate(barycenters):
        x, y = coord
        if abs(x-xmin) < tol:#
          ids_left.append(i)
        elif abs(x-xmax) < tol:
          ids_right.append(i)
        elif abs(y-ymin) < tol :
          ids_bottom.append(i)
        elif abs(y-ymax) < tol :
          ids_top.append(i)
        else:
            raise ValueError("Pb with boundary construction : barycenter does not belong to any border group")
        
    arr_left = mc.DataArrayInt(ids_left)
    arr_right = mc.DataArrayInt(ids_right)
    arr_bottom = mc.DataArrayInt(ids_bottom)
    arr_top = mc.DataArrayInt(ids_top)
    
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


if __name__ == """__main__""":
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=5,ny=5,mesh_name="squareWithRightTriangles_0")
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=10,ny=10,mesh_name="squareWithRightTriangles_1")
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=20,ny=20,mesh_name="squareWithRightTriangles_2")
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=50,ny=50,mesh_name="squareWithRightTriangles_3")
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=100,ny=100,mesh_name="squareWithRightTriangles_4")
    mesh_square_with_RightTriangles(xmin=0,xmax=1.,ymin=0,ymax=1.,nx=200,ny=200,mesh_name="squareWithRightTriangles_5")
