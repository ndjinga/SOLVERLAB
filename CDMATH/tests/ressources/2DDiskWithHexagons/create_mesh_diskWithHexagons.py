# -*-coding:utf-8 -*

#### Pavage hexagonal d'un cercle 
### input : xcenter, ycenter, Radius, ny
### output : HexagonsMesh.vtu

import medcoupling as mc
import MEDLoader as ML
import math

def mesh_disk_with_hexagons(xcenter=0., ycenter=0.,Radius=1.,ny=16,mesh_name="diskWithHexagons"):
    
    xmin=-Radius
    xmax=Radius
    ymin=-Radius
    ymax=Radius
    
    
    hradius=(ymax-ymin)/(ny*math.sqrt(3.))
    r = math.sqrt(3.)/2*hradius
    nx = int( 2*(xmax-xmin)/(3.*hradius) )
    
    # Building the coordinates of the initial hexagon, centered at 0,0
    d = mc.DataArrayDouble(6,2)
    d[:,0] = hradius
    d[:,1] = range(6)
    d[:,1] *= math.pi/3.
    d = d.fromPolarToCart()
    d.setInfoOnComponents(["X [m]","Y [m]"])
    
    print( "Uniform array ?", d.magnitude().isUniform(hradius,1e-12) )
    
    # translations of the first cell that are inside the circle
    translationToPerform = []
    for i in range(ny) :
        for j in range(nx):
            if (xcenter-xmin-(1.5*j+1)*hradius)**2+(ycenter-ymin-(2*i+(j%2)+1)*r)**2<Radius*Radius :
                translationToPerform.append([xmin+(1.5*j+1)*hradius,ymin+(2*i+(j%2)+1)*r] )
    
    ncells= len(translationToPerform) 
    print( "Meshing a disk with hexagons nx=",nx,"ny=",ny,"nb of cells=",ncells )
      
    ds = ncells*[None]
    for pos,t in enumerate(translationToPerform):
                     ds[pos] = d[:]         # Perform a deep copy of d and place it at position 'pos' in ds
                     ds[pos] += t             # Adding a vector to a set of coordinates does a translation
                     pass
    # Identifying duplicate tuples
    d2 = mc.DataArrayDouble.Aggregate(ds)
    oldNbOfTuples = d2.getNumberOfTuples()
    c,cI = d2.findCommonTuples(1e-12)
    tmp = c[cI[0]:cI[0+1]]
    print tmp
    a = cI.deltaShiftIndex()
    b = a - 1
    myNewNbOfTuples = oldNbOfTuples - sum(b.getValues())
    o2n, newNbOfTuples = mc.DataArrayInt.ConvertIndexArrayToO2N(oldNbOfTuples,c,cI)
    print( "Have I got the right number of tuples ?" )
    print( "myNewNbOfTuples = %d, newNbOfTuples = %d" % (myNewNbOfTuples, newNbOfTuples) )
    assert(myNewNbOfTuples == newNbOfTuples)
    print( "Old number of tuple was ", oldNbOfTuples )
    
    # Extracting the unique set of tuples
    d3 = d2.renumberAndReduce(o2n, newNbOfTuples)
    n2o = o2n.invertArrayO2N2N2O(newNbOfTuples)
    d3_bis = d2[n2o]
    print( "Are d3 and d3_bis equal ? %s" % (str(d3.isEqual(d3_bis, 1e-12))) )
    # Now build an unstructured mesh representing the final pattern
    mesh = mc.MEDCouplingUMesh(mesh_name,2)
    mesh.setCoords(d3)
    print( "Mesh dimension is", mesh.getMeshDimension() )
    print( "Spatial dimension is", mesh.getCoords().getNumberOfComponents() )
    mesh.allocateCells(ncells)
    for i in range(ncells):
            cell_connec = o2n[6*i,6*i+1,6*i+2,6*i+3,6*i+4,6*i+5]
            mesh.insertNextCell(mc.NORM_POLYGON, cell_connec.getValues())
            pass
    
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
    arr_circle = mc.DataArrayInt(range(mesh_1d.getNumberOfCells()))
    arr_circle.setName("Circle")
    meshMEDFile.addGroup(-1, arr_circle)
    
    filename = mesh_name+".med"
    # Write the result into a VTU file that can be read with ParaView
    #mesh.writeVTK(mesh_name".vtu")
    # Write the result into a MED file that can be read with Salomé
    meshMEDFile.write(filename,2) # 2 stands for write from scratch

if __name__ == '__main__':
  mesh_disk_with_hexagons(0.,0.,1.,14,"diskWithHexagons")
    
