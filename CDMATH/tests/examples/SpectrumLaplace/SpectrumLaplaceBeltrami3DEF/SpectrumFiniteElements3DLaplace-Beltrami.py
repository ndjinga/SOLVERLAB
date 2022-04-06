# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Calcul EF du spectre de l'opérateur de Laplace-Beltrami -\triangle sur une surface en 3D
# Author      : Michael Ndjinga
# Copyright   : CEA Saclay 2020
# Description : Utilisation de la méthode des éléménts finis P1 avec champs discrétisés aux noeuds d'un maillage triangulaire
#               Création et sauvegarde des champs résultant en utilisant la librairie CDMATH
#================================================================================================================================

import cdmath
import sys

#Chargement du maillage triangulaire de la surface
#=================================================
my_mesh = cdmath.Mesh(sys.argv[1])
mesh_name=sys.argv[2]
if(not my_mesh.isTriangular()) :
	raise ValueError("Wrong cell types : mesh is not made of triangles")
if(my_mesh.getMeshDimension()!=2) :
	raise ValueError("Wrong mesh dimension : expected a surface of dimension 2")
if(my_mesh.getSpaceDimension()!=3) :
	raise ValueError("Wrong space dimension : expected a space of dimension 3")

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh building/loading done")
print("nb of nodes=", nbNodes)
print("nb of cells=", nbCells)

maxNbNeighbours = my_mesh.getMaxNbNeighbours(cdmath.NODES)+1#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix

# Construction de la matrice de rigidité
#=======================================
Rigidite=cdmath.SparseMatrixPetsc(nbNodes,nbNodes,maxNbNeighbours)# warning : third argument is number of non zero coefficients per line

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un triangle
GradShapeFunc0=cdmath.Vector(3)
GradShapeFunc1=cdmath.Vector(3)
GradShapeFunc2=cdmath.Vector(3)

normalFace0=cdmath.Vector(3)
normalFace1=cdmath.Vector(3)

nodal_volumes=cdmath.Vector(nbNodes)

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

	#Build normal to cell Ci
	normalFace0[0]=Ci.getNormalVector(0,0)
	normalFace0[1]=Ci.getNormalVector(0,1)
	normalFace0[2]=Ci.getNormalVector(0,2)
	normalFace1[0]=Ci.getNormalVector(1,0)
	normalFace1[1]=Ci.getNormalVector(1,1)
	normalFace1[2]=Ci.getNormalVector(1,2)

	normalCell = normalFace0.crossProduct(normalFace1)
	normalCell = normalCell*(1/normalCell.norm())

	cellMat=cdmath.Matrix(4)
	cellMat[0,0]=N0.x()
	cellMat[0,1]=N0.y()
	cellMat[0,2]=N0.z()
	cellMat[1,0]=N1.x()
	cellMat[1,1]=N1.y()
	cellMat[1,2]=N1.z()
	cellMat[2,0]=N2.x()
	cellMat[2,1]=N2.y()
	cellMat[2,2]=N2.z()
	cellMat[3,0]=normalCell[0]
	cellMat[3,1]=normalCell[1]
	cellMat[3,2]=normalCell[2]
	cellMat[0,3]=1
	cellMat[1,3]=1
	cellMat[2,3]=1
	cellMat[3,3]=0

	#Formule des gradients voir EF P1 -> calcul déterminants
	GradShapeFunc0[0]= cellMat.partMatrix(0,0).determinant()*0.5
	GradShapeFunc0[1]=-cellMat.partMatrix(0,1).determinant()*0.5
	GradShapeFunc0[2]= cellMat.partMatrix(0,2).determinant()*0.5
	GradShapeFunc1[0]=-cellMat.partMatrix(1,0).determinant()*0.5
	GradShapeFunc1[1]= cellMat.partMatrix(1,1).determinant()*0.5
	GradShapeFunc1[2]=-cellMat.partMatrix(1,2).determinant()*0.5
	GradShapeFunc2[0]= cellMat.partMatrix(2,0).determinant()*0.5
	GradShapeFunc2[1]=-cellMat.partMatrix(2,1).determinant()*0.5
	GradShapeFunc2[2]= cellMat.partMatrix(2,2).determinant()*0.5

	#Création d'un tableau (numéro du noeud, gradient de la fonction de forme
	GradShapeFuncs={nodeId0 : GradShapeFunc0}
	GradShapeFuncs[nodeId1]=GradShapeFunc1
	GradShapeFuncs[nodeId2]=GradShapeFunc2

	# Remplissage de  la matrice de rigidité et du second membre
	for j in [nodeId0,nodeId1,nodeId2] : 
		nodal_volumes[j]+=Ci.getMeasure()/3
		#Contribution de la cellule triangulaire i à la ligne j du système linéaire
		for k in [nodeId0,nodeId1,nodeId2] : 
			Rigidite.addValue(j,k,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())

print("Linear system matrix building done")
Rigidite.viewNonZeroStructure( 0, "RigidityMatrix_FiniteElementsOn"+mesh_name+"LaplaceBeltrami")
Rigidite.plotEigenvalues("FiniteElementsOn"+mesh_name+"LaplaceBeltrami")

# Conditionnement de la matrice de rigidité
#==========================================
cond = Rigidite.getConditionNumber(True)
print("Condition number is ",cond)

# Spectre de la matrice de rigidité
#==================================
#Homogénéisation de la matrice de rigidité (sinon elle tend vers zero de meme que ses valeurs propres)
for i in range(nbNodes):
	nodal_volumes[i]=1/nodal_volumes[i]
Rigidite.leftDiagonalScale(nodal_volumes)

nev=9
d=Rigidite.getEigenvectorsDataArrayDouble(nev)
my_eigenfield = cdmath.Field("Eigenvectors field", cdmath.NODES, my_mesh, nev)
my_eigenfield.setFieldByDataArrayDouble(d)
# Free memory
d.decrRef()
    
# Sauvegarde du champ résultat
#===========================
my_eigenfield.writeVTK("spectrumFiniteElementsOn"+mesh_name+"LaplaceBeltrami")
