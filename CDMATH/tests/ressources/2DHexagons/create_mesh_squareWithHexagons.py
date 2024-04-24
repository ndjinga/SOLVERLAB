# -*-coding:utf-8 -*

#### Pavage hexagonal d'un rectangle 
### input : xmin, xmax, ymin, ymax, ny
### output : squareWithHexagons.vtu, squareWithHexagons.med

import medcoupling as mc
import math
import MEDLoader as ML

def mesh_square_with_hexagons(xmin=0,xmax=1,ymin=0,ymax=1,ny=14,mesh_name="squareWithHexagons"):
    
    radius=(ymax-ymin)/(ny*math.sqrt(3.))
    r = math.sqrt(3.)/2*radius
    nx = int( 2*(xmax-xmin)/(3.*radius) )
    
    print( "Meshing a square with hexagons nx=",nx,"ny=",ny, "ncells=",nx*ny )
    
    # Building the coordinates of the initial hexagon, centered at 0,0
    d = mc.DataArrayDouble(6,2)
    d[:,0] = radius
    d[:,1] = list(range(6))
    d[:,1] *= math.pi/3.
    d = d.fromPolarToCart()
    d.setInfoOnComponents(["X [m]","Y [m]"])
    
    print( "Uniform array ?", d.magnitude().isUniform(radius,1e-12) )
    
    # translations of the first cell
    translationToPerform = [[xmin+(1.5*j+1)*radius,ymin+(2*i+(j%2)+1)*r] for i in range(ny) for j in range(nx)]
        
    ds = len(translationToPerform)*[None]
    for pos,t in enumerate(translationToPerform):
                     ds[pos] = d[:]         # Perform a deep copy of d and place it at position 'pos' in ds
                     ds[pos] += t             # Adding a vector to a set of coordinates does a translation
                     pass
    
    d2 = mc.DataArrayDouble.Aggregate(ds)
    # Build an unstructured mesh representing the final pattern
    mesh = mc.MEDCouplingUMesh(mesh_name,2)
    mesh.setCoords(d2)
    print( "Mesh dimension is", mesh.getMeshDimension() )
    print( "Spatial dimension is", mesh.getCoords().getNumberOfComponents() )
    mesh.allocateCells(nx*ny)
    for i in range(nx*ny):
            cell_connec = [6*i,6*i+1,6*i+2,6*i+3,6*i+4,6*i+5]
            mesh.insertNextCell(mc.NORM_POLYGON, cell_connec)
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
    for i, coord in enumerate(barycenters):
        x, y = coord
        if abs(x-xmin-radius/4) < tol:#
          ids_left.append(i)
        elif abs(x-xmin-(1.5*nx+0.25)*radius) < tol:
          ids_right.append(i)
        elif abs(y-ymin) < tol or abs(y-ymin-r) < tol or abs(y-ymin-r/2) < tol:
          ids_bottom.append(i)
        elif abs(y-ymin-2*r*ny) < tol or abs(y-ymin-2*r*ny-r) < tol or abs(y-ymin-2*r*ny-r/2) < tol:
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

if __name__ == '__main__':
  mesh_square_with_hexagons(0.,1.,0.,1.,14,"squareWithHexagons")
    
