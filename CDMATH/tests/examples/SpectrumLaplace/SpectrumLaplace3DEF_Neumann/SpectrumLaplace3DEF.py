# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Calcul EF du spectre de l'opérateur de Laplace 3D -\triangle avec conditions aux limites de Neumann \nabla u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Utilisation de la méthode des éléménts finis P1 avec champs discrétisés aux noeuds d'un maillage de tetraèdres
#                Création et sauvegarde des champs résultant en utilisant la librairie CDMATH
#                Résultats à consolider
#================================================================================================================================

import cdmath
import sys

if len(sys.argv) >2 :#load a mesh file
    my_mesh = cdmath.Mesh(sys.argv[1])
    mesh_name=sys.argv[2]
else :   #cuboid mesh split into triangles
    xmin=0
    xmax=1
    ymin=0
    ymax=1
    zmin=0
    zmax=1
    
    nx=15
    ny=15
    nz=15
    
    my_mesh = cdmath.Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,0)
    mesh_name="RightTetrahedra"

    eps=1e-6
    my_mesh.setGroupAtPlan(0.,0,eps,"Border")#Bord ARRIERE
    my_mesh.setGroupAtPlan(1.,0,eps,"Border")#Bord AVANT
    my_mesh.setGroupAtPlan(0.,1,eps,"Border")#Bord BAS
    my_mesh.setGroupAtPlan(1.,1,eps,"Border")#Bord HAUT
    my_mesh.setGroupAtPlan(0.,2,eps,"Border")#Bord BAS
    my_mesh.setGroupAtPlan(1.,2,eps,"Border")#Bord HAUT

if( my_mesh.getSpaceDimension()!=3 or my_mesh.getMeshDimension()!=3) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")
if(not my_mesh.isTetrahedral()) :
    raise ValueError("Wrong cell types : mesh is not made of tetrahedra")

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("Number of nodes=", nbNodes)
print("Number of cells=", nbCells)

#Détermination des noeuds intérieurs
#===================================
maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix

#parcours des noeuds pour déterminer nb max voisins d'un noeud
for i in range(nbNodes):
    Ni=my_mesh.getNode(i)
    if not my_mesh.isBorderNode(i): # Détection des noeuds frontière
        maxNbNeighbours= max(1+Ni.getNumberOfEdges(),maxNbNeighbours) 

print("Max nb of neighbours=", maxNbNeighbours)

# Construction de la matrice de rigidité
#========================================
Rigidite=cdmath.SparseMatrixPetsc(nbNodes,nbNodes,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix

nodal_volumes=cdmath.Vector(nbNodes)#Volumes of the dual cells associated to each node

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un tétraèdre
GradShapeFunc0=cdmath.Vector(3)
GradShapeFunc1=cdmath.Vector(3)
GradShapeFunc2=cdmath.Vector(3)
GradShapeFunc3=cdmath.Vector(3)

#On parcourt les tétraèdres du domaine
for i in range(nbCells):

	Ci=my_mesh.getCell(i)

	#Extraction des noeuds de la cellule
	nodeId0=Ci.getNodeId(0)
	nodeId1=Ci.getNodeId(1)
	nodeId2=Ci.getNodeId(2)
	nodeId3=Ci.getNodeId(3)
	N0=my_mesh.getNode(nodeId0)
	N1=my_mesh.getNode(nodeId1)
	N2=my_mesh.getNode(nodeId2)
	N3=my_mesh.getNode(nodeId3)

	#Formule des gradients voir EF P1 -> calcul déterminants
	GradShapeFunc0[0]= (N2.y()*N3.z()-N2.z()*N3.y()-N1.y()*N3.z()+N3.y()*N1.z()+N1.y()*N2.z()-N2.y()*N1.z())/6
	GradShapeFunc0[1]=-(N2.x()*N3.z()-N2.z()*N3.x()-N1.x()*N3.z()+N3.x()*N1.z()+N1.x()*N2.z()-N2.x()*N1.z())/6
	GradShapeFunc0[2]=(N2.x()*N3.y()-N2.y()*N3.x()-N1.x()*N3.y()+N3.x()*N1.y()+N1.x()*N2.y()-N2.x()*N1.y())/6
	GradShapeFunc1[0]=- (N2.y()*N3.z()-N2.z()*N3.y()-N0.y()*N3.z()+N3.y()*N0.z()+N0.y()*N2.z()-N2.y()*N0.z())/6
	GradShapeFunc1[1]=(N2.x()*N3.z()-N2.z()*N3.x()-N0.x()*N3.z()+N3.x()*N0.z()+N0.x()*N2.z()-N2.x()*N0.z())/6
	GradShapeFunc1[2]=-(N2.x()*N3.y()-N2.y()*N3.x()-N0.x()*N3.y()+N3.x()*N0.y()+N0.x()*N2.y()-N2.x()*N0.y())/6
	GradShapeFunc2[0]= -(N0.y()*N3.z()-N0.z()*N3.y()-N1.y()*N3.z()+N3.y()*N1.z()+N1.y()*N0.z()-N0.y()*N1.z())/6
	GradShapeFunc2[1]=(N0.x()*N3.z()-N0.z()*N3.x()-N1.x()*N3.z()+N3.x()*N1.z()+N1.x()*N0.z()-N0.x()*N1.z())/6
	GradShapeFunc2[2]= -(N0.x()*N3.y()-N0.y()*N3.x()-N1.x()*N3.y()+N3.x()*N1.y()+N1.x()*N0.y()-N0.x()*N1.y())/6
	GradShapeFunc3[0]=-(N2.y()*N0.z()-N2.z()*N0.y()-N1.y()*N0.z()+N0.y()*N1.z()+N1.y()*N2.z()-N2.y()*N1.z())/6
	GradShapeFunc3[1]=(N2.x()*N0.z()-N2.z()*N0.x()-N1.x()*N0.z()+N0.x()*N1.z()+N1.x()*N2.z()-N2.x()*N1.z())/6
	GradShapeFunc3[2]=-(N2.x()*N0.y()-N2.y()*N0.x()-N1.x()*N0.y()+N0.x()*N1.y()+N1.x()*N2.y()-N2.x()*N1.y())/6
	
	#Création d'un tableau (numéro du noeud, gradient de la fonction de forme)
	GradShapeFuncs={nodeId0 : GradShapeFunc0}
	GradShapeFuncs[nodeId1]=GradShapeFunc1
	GradShapeFuncs[nodeId2]=GradShapeFunc2
	GradShapeFuncs[nodeId3]=GradShapeFunc3

	# Remplissage de  la matrice de rigidité et du second membre
	for j in [nodeId0,nodeId1,nodeId2,nodeId3] :
		#Calcul du volume de la cellule duale pour normalisation des vecteurs propres
		nodal_volumes[j]+=Ci.getMeasure()/4
		#Contribution de la cellule tétraédrique i à la ligne j_int du système linéaire
		for k in [nodeId0,nodeId1,nodeId2,nodeId3] :
			Rigidite.addValue(j,k,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())


print("Stiffness matrix construction done")
Rigidite.viewNonZeroStructure( 0, "RigidityMatrix_FiniteElementsOn"+mesh_name+"Laplace")
Rigidite.plotEigenvalues("FiniteElementsOn"+mesh_name+"Laplace")

# Conditionnement de la matrice de rigidité
#=================================
cond = Rigidite.getConditionNumber()
print("Condition number is ",cond)

# Spectre de la matrice de rigidité
#==================================
#Homogénéisation de la matrice de rigidité (sinon elle tend vers zero de meme que ses valeurs propres)
for i in range(nbNodes):
    nodal_volumes[i]=1/nodal_volumes[i]
Rigidite.leftDiagonalScale(nodal_volumes)

nev=10
d=Rigidite.getEigenvectorsDataArrayDouble(nev)
my_eigenfield = cdmath.Field("Eigenvectors field", cdmath.NODES, my_mesh, nev)
for j in range(nbNodes):
    for k in range(nev):
      my_eigenfield[j,k]=d[j,k];#remplissage des valeurs pour les noeuds intérieurs
for k in range(nev):
    my_eigenfield.setInfoOnComponent(k,d.getInfoOnComponent(k))
# Free memory
d.decrRef()
    
# Sauvegarde du champ résultat
#===========================
my_eigenfield.writeVTK("spectrumFiniteElementsOn"+mesh_name+"Laplace")
my_eigenfield.writeMED("spectrumFiniteElementsOn"+mesh_name+"Laplace")
