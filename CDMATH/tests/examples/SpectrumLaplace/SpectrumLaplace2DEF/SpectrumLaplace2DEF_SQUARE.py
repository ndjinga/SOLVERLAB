# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Calcul EF du spectre de l'opérateur de Laplace 2D -\triangle avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2020
# Description : Utilisation de la méthode des éléménts finis P1 avec champs discrétisés aux noeuds d'un maillage triangulaire
#		        Création et sauvegarde des champs résultant en utilisant la librairie CDMATH
#================================================================================================================================

import cdmath
import sys

if len(sys.argv) >2 :#load a mesh file
    my_mesh = cdmath.Mesh(sys.argv[1])
    mesh_name=sys.argv[2]
else :   #rectangular mesh split into triangles
    xmin=0
    xmax=1
    ymin=0
    ymax=1
    
    nx=15
    ny=15
    
    my_mesh = cdmath.Mesh(xmin,xmax,nx,ymin,ymax,ny,0)
    mesh_name="RightTriangles"

if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")
if(not my_mesh.isTriangular()) :
	raise ValueError("Wrong cell types : mesh is not made of triangles")
eps=1e-6
my_mesh.setGroupAtPlan(0.,0,eps,"DirichletBorder")#Bord GAUCHE
my_mesh.setGroupAtPlan(1.,0,eps,"DirichletBorder")#Bord DROIT
my_mesh.setGroupAtPlan(0.,1,eps,"DirichletBorder")#Bord BAS
my_mesh.setGroupAtPlan(1.,1,eps,"DirichletBorder")#Bord HAUT

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("Number of nodes=", nbNodes)
print("Number of cells=", nbCells)

#Détermination des noeuds intérieurs
#===================================
nbInteriorNodes = 0
nbBoundaryNodes = 0
maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
interiorNodes=[]
boundaryNodes=[]

#parcours des noeuds pour discrétisation du second membre et extraction 1) des noeuds intérieur 2) des noeuds frontière 3) du nb max voisins d'un noeud
for i in range(nbNodes):
	Ni=my_mesh.getNode(i)
	if my_mesh.isBorderNode(i): # Détection des noeuds frontière
		boundaryNodes.append(i)
		nbBoundaryNodes=nbBoundaryNodes+1
	else: # Détection des noeuds intérieurs
		interiorNodes.append(i)
		nbInteriorNodes=nbInteriorNodes+1
		maxNbNeighbours= max(1+Ni.getNumberOfCells(),maxNbNeighbours) #true only in 2D, otherwise use function Ni.getNumberOfEdges()

print("nb of interior nodes=", nbInteriorNodes)
print("nb of boundary nodes=", nbBoundaryNodes)
print("Max nb of neighbours=", maxNbNeighbours)

# Construction de la matrice de rigidité
#========================================
Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un triangle (hypothèse 2D)
GradShapeFunc0=cdmath.Vector(2)
GradShapeFunc1=cdmath.Vector(2)
GradShapeFunc2=cdmath.Vector(2)

nodal_volumes=cdmath.Vector(nbInteriorNodes)

#On parcourt les triangles du domaine
for i in range(nbCells):

	Ci=my_mesh.getCell(i)

	#Contribution à la matrice de rigidité
	nodeId0=Ci.getNodeId(0)
	nodeId1=Ci.getNodeId(1)
	nodeId2=Ci.getNodeId(2)

	N0=my_mesh.getNode(nodeId0)
	N1=my_mesh.getNode(nodeId1)
	N2=my_mesh.getNode(nodeId2)

	#Formule des gradients voir EF P1 -> calcul déterminants
	GradShapeFunc0[0]= (N1.y()-N2.y())*0.5
	GradShapeFunc0[1]=-(N1.x()-N2.x())*0.5
	GradShapeFunc1[0]=-(N0.y()-N2.y())*0.5
	GradShapeFunc1[1]= (N0.x()-N2.x())*0.5
	GradShapeFunc2[0]= (N0.y()-N1.y())*0.5
	GradShapeFunc2[1]=-(N0.x()-N1.x())*0.5

	#Création d'un tableau (numéro du noeud, gradient de la fonction de forme
	GradShapeFuncs={nodeId0 : GradShapeFunc0}
	GradShapeFuncs[nodeId1]=GradShapeFunc1
	GradShapeFuncs[nodeId2]=GradShapeFunc2


	# Remplissage de  la matrice de rigidité et du second membre
	for j in [nodeId0,nodeId1,nodeId2] :
		if boundaryNodes.count(j)==0 : #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
			j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
			nodal_volumes[j_int]+=Ci.getMeasure()/3
			#Contribution de la cellule triangulaire i à la ligne j_int du système linéaire
			for k in [nodeId0,nodeId1,nodeId2] : 
				if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
					k_int=interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
					Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
				#else: si condition limite non nulle au bord, ajouter la contribution du bord au second membre de la cellule j

print("Stiffness matrix construction done")
quit()
# Conditionnement de la matrice de rigidité
#=================================
cond = Rigidite.getConditionNumber()
print("Condition number is ",cond)

# Spectre de la matrice de rigidité
#==================================
#Homogénéisation de la matrice de rigidité (sinon elle tend vers zero de meme que ses valeurs propres)
for i in range(nbInteriorNodes):
	nodal_volumes[i]=1/nodal_volumes[i]
Rigidite.leftDiagonalScale(nodal_volumes)

nev=10
d=Rigidite.getEigenvectorsDataArrayDouble(nev)
my_eigenfield = cdmath.Field("Eigenvectors field", cdmath.NODES, my_mesh, nev)
for j in range(nbInteriorNodes):
    for k in range(nev):
      my_eigenfield[interiorNodes[j],k]=d[j,k];#remplissage des valeurs pour les noeuds intérieurs
for j in range(nbBoundaryNodes):
    for k in range(nev):
      my_eigenfield[boundaryNodes[j],k]=0;#remplissage des valeurs pour les noeuds frontière (condition limite)
for k in range(nev):
    my_eigenfield.setInfoOnComponent(k,d.getInfoOnComponent(k))
    
# Sauvegarde du champ résultat
#===========================
my_eigenfield.writeVTK("spectrumFiniteElementsOn"+mesh_name+"Laplace")
