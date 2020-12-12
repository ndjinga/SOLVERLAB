# -*-coding:utf-8 -*

#### Pavage triangulaire d'un rectangle avec des triangles équilatéraux
### input : xmin, xmax, ymin, ymax, ny
### output : squareWithEquilateralTriangles.vtu, squareWithEquilateralTriangles.med

import medcoupling as mc
import math
import MEDLoader as ML

def mesh_square_with_EquilateralTriangles(xmin=0,xmax=1,ymin=0,ymax=1,ny=10,mesh_name="squareWithEquilateralTriangles"):
    
    h = (ymax-ymin)/ny #hauteur d'un seul triangle
    lenght = 2*h/math.sqrt(3.) # longueur du coté d'un triangle
    radius=2*h/3. # rayon du cercle circonscrit
    nx = int( (xmax-xmin)/lenght )
    
    print( "Meshing a square with Equilateral Triangles nx=",nx,"ny=",ny, "ncells=",nx*ny)
    
    # Building the coordinates of the initial triangle, centered at 0,0 and pointing upward
    d_up = mc.DataArrayDouble(3,2)
    d_up[0,0] = radius*math.sqrt(3.)/2
    d_up[0,1] =-radius/2
    d_up[1,0] = 0
    d_up[1,1] = radius
    d_up[2,0] =-radius*math.sqrt(3.)/2
    d_up[2,1] =-radius/2
    d_up.setInfoOnComponents(["X [m]","Y [m]"])
    
    # Building the coordinates of the initial triangle, centered at 0,0 and pointing downward (symetry around x axis : change the sign of the y coordinate)
    d_down = mc.DataArrayDouble(3,2)
    d_down[0,0] = radius*math.sqrt(3.)/2
    d_down[0,1] = radius/2
    d_down[1,0] = 0
    d_down[1,1] =-radius
    d_down[2,0] =-radius*math.sqrt(3.)/2
    d_down[2,1] = radius/2
    d_down.setInfoOnComponents(["X [m]","Y [m]"])
    
    print( "Uniform arrays ?", d_up.magnitude().isUniform(radius,1e-12), d_down.magnitude().isUniform(radius,1e-12) )
    
    # translations of the cells : 
    translationToPerform_up   = [[xmin+radius*math.sqrt(3.)/2 +     (i%2)*lenght/2 + lenght*j,ymin+  radius/2+h*i] for i in range(ny) for j in range(nx)]
    translationToPerform_down = [[xmin+radius*math.sqrt(3.)/2 + ((i+1)%2)*lenght/2 + lenght*j,ymin+h-radius/2+h*i] for i in range(ny) for j in range(nx)]
        
    ds = (len(translationToPerform_up)+len(translationToPerform_down))*[None]
    for pos,t_up in enumerate(translationToPerform_up):
                     ds[2*pos]  = d_up[:]        # Perform a deep copy of d_up and place it at position 'pos' in ds
                     ds[2*pos] += t_up           # Adding a vector to a set of coordinates does a translation
                     pass
    for pos,t_down in enumerate(translationToPerform_down):
                     ds[2*pos+1]  = d_down[:]        # Perform a deep copy of d_down and place it at position 'pos' in ds
                     ds[2*pos+1] += t_down           # Adding a vector to a set of coordinates does a translation
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
        if abs(x-xmin-lenght/4) < tol:#
          ids_left.append(i)
        elif abs(x-xmin-lenght*nx-lenght/4) < tol:
          ids_right.append(i)
        elif abs(y-ymin) < tol :
          ids_bottom.append(i)
        elif abs(y-ymax) < tol :
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
  mesh_square_with_EquilateralTriangles(0.,1.,0.,1.,5,"squareWithEquilateralTriangles5")
  mesh_square_with_EquilateralTriangles(0.,1.,0.,1.,20,"squareWithEquilateralTriangles20")
  mesh_square_with_EquilateralTriangles(0.,1.,0.,1.,50,"squareWithEquilateralTriangles50")
  mesh_square_with_EquilateralTriangles(0.,1.,0.,1.,100,"squareWithEquilateralTriangles100")
  mesh_square_with_EquilateralTriangles(0.,1.,0.,1.,200,"squareWithEquilateralTriangles200")
    
