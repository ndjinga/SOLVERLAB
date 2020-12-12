# -*-coding:utf-8 -*

#### Maillage d'un cercle par une grille de carrés
### input : xcenter, ycenter, Radius, n
### output : CircleWithSquares.vtu

import medcoupling as mc
import math
import MEDLoader as ML

def mesh_disk_with_squares(xcenter=0.,ycenter=0., Radius=1.,n=17,mesh_name="diskWithSquares"):

    xmin=-Radius
    xmax=Radius
    ymin=-Radius
    ymax=Radius
    
    dx = (xmax-xmin)/n
    dy=(ymax-ymin)/n
    
    # Building the initial rectangular cell, centered at 0,0
    d = mc.DataArrayDouble(4,2)
    d[0,0] = -dx/2
    d[0,1] =  dy/2
    d[1,0] =  dx/2
    d[1,1] =  dy/2
    d[2,0] =  dx/2
    d[2,1] = -dy/2
    d[3,0] = -dx/2
    d[3,1] = -dy/2
    d.setInfoOnComponents(["X [m]","Y [m]"])
    
    print( "Uniform array ?", d.magnitude().isUniform(0.5*math.sqrt(dx*dx+dy*dy),1e-10) )
    
    # translation of the first cell
    translationToPerform = []
    for i in range(n) :
        for j in range(n):
            if (xcenter-xmin-(0.5+i)*dx)**2+(ycenter-ymin-(0.5+j)*dy)**2<Radius*Radius :
                translationToPerform.append([xmin+(0.5+i)*dx,ymin+(0.5+j)*dy] )
    
    ncells= len(translationToPerform) 
    print( "Meshing a disk with squares ",n," nb of cells=",ncells )
      
    ds = ncells*[None]
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
    mesh.allocateCells(ncells)
    for i in range(ncells):
            cell_connec = [4*i,4*i+1,4*i+2,4*i+3)
            mesh.insertNextCell(mc.NORM_QUAD4, cell_connec)
            pass
    
    # Identifying duplicate nodes
    oldNbOfNodes=mesh.getNumberOfNodes()        
    arr, areNodesMerged, newNbOfNodes=mesh.mergeNodes(1e-10)
    print( "oldNbOfNodes=",oldNbOfNodes,"newNbOfNodes",newNbOfNodes )
    
    # Check that everything is coherent (will throw if not)
    mesh.checkConsistencyLight()
    
    # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
    mesh_1d = mesh.computeSkin()
    
    # Trie les cellules par type conformément à la convention MED fichier
    o2n = mesh.sortCellsInMEDFileFrmt()
    meshMEDFile=ML.MEDFileUMesh.New()
    # Ecrit le maillage 2D
    meshMEDFile.setMeshAtLevel(0,mesh)
    # Ecrit le maillage 1D
    meshMEDFile.setMeshAtLevel(-1,mesh_1d)
    # Ecrit les groupes
    arr_circle = mc.DataArrayIdType(range(mesh_1d.getNumberOfCells()))
    arr_circle.setName("Circle")
    meshMEDFile.addGroup(-1, arr_circle)
    
    filename = mesh_name+".med"
    # Write the result into a VTU file that can be read with ParaView
    #mesh.writeVTK(mesh_name+".vtu")
    # Write the result into a MED file that can be read with Salomé
    meshMEDFile.write(filename,2) # 2 stands for write from scratch
    
if __name__ == '__main__':
  mesh_disk_with_squares(0.,0.,1.,14,"diskWithSquares")
