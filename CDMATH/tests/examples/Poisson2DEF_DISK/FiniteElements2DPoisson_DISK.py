# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 2D -\triangle u = f sur le disque unité avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage triangulaire
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=-(r-1)*r**2*sin(theta)
#================================================================================================================================

import cdmath
from math import sin, sqrt, atan2
import numpy as np
import matplotlib.pyplot as plt
import PV_routines
import VTK_routines

#Chargement du maillage triangulaire du disque unité
#=======================================================================================
my_mesh = cdmath.Mesh("diskWithTriangles.med")
if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")
if(not my_mesh.isTriangular()) :
	raise ValueError("Wrong cell types : mesh is not made of triangles")

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("Number of nodes=", nbNodes)
print("Number of cells=", nbCells)

#Discrétisation du second membre et détermination des noeuds intérieurs
#======================================================================
my_RHSfield = cdmath.Field("RHS_field", cdmath.NODES, my_mesh, 1)
my_ExactSol = cdmath.Field("Exact_field", cdmath.NODES, my_mesh, 1)
nbInteriorNodes = 0
nbBoundaryNodes = 0
maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
interiorNodes=[]
boundaryNodes=[]

#parcours des noeuds pour discrétisation du second membre et extraction 1) des noeuds intérieur 2) des noeuds frontière 3) du nb max voisins d'un noeud
for i in range(nbNodes):
	Ni=my_mesh.getNode(i)
	x = Ni.x()
	y = Ni.y()

	r=sqrt(x*x+y*y)
	theta=atan2(y,x)

	my_RHSfield[i]=(8*r-3)*sin(theta)#mettre la fonction definie au second membre de l'edp
	my_ExactSol[i]=-(r-1)*r*r*sin(theta)#mettre la solution exacte de l'edp
    
	if my_mesh.isBorderNode(i): # Détection des noeuds frontière
		boundaryNodes.append(i)
		nbBoundaryNodes=nbBoundaryNodes+1
	else: # Détection des noeuds intérieurs
		interiorNodes.append(i)
		nbInteriorNodes=nbInteriorNodes+1
		maxNbNeighbours= max(1+Ni.getNumberOfCells(),maxNbNeighbours) #true only in 2D, otherwise use function Ni.getNumberOfEdges()

print("Right hand side discretisation done")
print("Number of interior nodes=", nbInteriorNodes)
print("Number of boundary nodes=", nbBoundaryNodes)
print("Maximum number of neighbours=", maxNbNeighbours)

# Construction de la matrice de rigidité et du vecteur second membre du système linéaire
#=======================================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbInteriorNodes)

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un triangle (hypothèse 2D)
GradShapeFunc0=cdmath.Vector(2)
GradShapeFunc1=cdmath.Vector(2)
GradShapeFunc2=cdmath.Vector(2)

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
	GradShapeFunc0[0]= (N1.y()-N2.y())/2
	GradShapeFunc0[1]=-(N1.x()-N2.x())/2
	GradShapeFunc1[0]=-(N0.y()-N2.y())/2
	GradShapeFunc1[1]= (N0.x()-N2.x())/2
	GradShapeFunc2[0]= (N0.y()-N1.y())/2
	GradShapeFunc2[1]=-(N0.x()-N1.x())/2

	#Création d'un tableau (numéro du noeud, gradient de la fonction de forme
	GradShapeFuncs={nodeId0 : GradShapeFunc0}
	GradShapeFuncs[nodeId1]=GradShapeFunc1
	GradShapeFuncs[nodeId2]=GradShapeFunc2


	# Remplissage de  la matrice de rigidité et du second membre
	for j in [nodeId0,nodeId1,nodeId2] :
		if boundaryNodes.count(j)==0 : #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
			j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
			#Ajout de la contribution de la cellule triangulaire i au second membre du noeud j 
			RHS[j_int]=Ci.getMeasure()/3*my_RHSfield[j]+RHS[j_int] # intégrale dans le triangle du produit f x fonction de base
			#Contribution de la cellule triangulaire i à la ligne j_int du système linéaire
			for k in [nodeId0,nodeId1,nodeId2] : 
				if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
					k_int=interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
					Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
				#else: si condition limite non nulle au bord, ajouter la contribution du bord au second membre de la cellule j

print("Linear system matrix building done")

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"CG","ILU")#Remplacer CG par CHOLESKY pour solveur direct
SolSyst=LS.solve()

print( "Preconditioner used : ", LS.getNameOfPc() )
print( "Number of iterations used : ", LS.getNumberOfIter() )
print( "Final residual : ", LS.getResidu() )
print("Linear system solved")

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("ResultField", cdmath.NODES, my_mesh, 1)
for j in range(nbInteriorNodes):
    my_ResultField[interiorNodes[j]]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs
for j in range(nbBoundaryNodes):
    my_ResultField[boundaryNodes[j]]=0;#remplissage des valeurs pour les noeuds frontière (condition limite)
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteElements2DPoisson_DISK_ResultField")

# Postprocessing :
#=================
# save 2D picture
PV_routines.Save_PV_data_to_picture_file("FiniteElements2DPoisson_DISK_ResultField"+'_0.vtu',"ResultField",'NODES',"FiniteElements2DPoisson_DISK_ResultField")

# extract and plot diagonal values
resolution=100
curv_abs=np.linspace(-1,1,resolution+1)
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[-0.5,-0.5,0],[0.5,0.5,0], resolution)
plt.plot(curv_abs, diag_data, label= '2D disk mesh with '+str(nbNodes) + ' nodes')
plt.legend()
plt.xlabel('Position on diagonal line')
plt.ylabel('Value on diagonal line')
plt.title('Plot over diagonal line for finite elements \n for Laplace operator on a 2D disk triangular mesh')
plt.savefig("FiniteElements2DPoisson_DISK_ResultField_"+str(nbNodes) + '_nodes'+"_PlotOverDiagonalLine.png")

print("Numerical solution of 2D Poisson equation on a disk using finite elements done")

#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
max_abs_sol_exacte=max(my_ExactSol.max(),-my_ExactSol.min())
max_sol_num=my_ResultField.max()
min_sol_num=my_ResultField.min()
erreur_abs=0
for i in range(nbNodes) :
    if  erreur_abs < abs(my_ExactSol[i] - my_ResultField[i]) :
        erreur_abs = abs(my_ExactSol[i] - my_ResultField[i])

print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
print("Maximum exact solution = ", my_ExactSol.max(), " Minimum exact solution = ", my_ExactSol.min() )

assert erreur_abs/max_abs_sol_exacte <1.
