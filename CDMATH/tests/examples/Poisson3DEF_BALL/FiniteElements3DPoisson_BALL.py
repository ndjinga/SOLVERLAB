# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 3D -\triangle u = f sur la boule unité avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage tétraédrique
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=-(r-1)*r**2*sin(theta)**2*cos(phi)
#================================================================================================================================

import cdmath
from math import sin, cos, atan2, sqrt
import numpy as np
import matplotlib.pyplot as plt
import PV_routines
import VTK_routines

#Chargement du maillage tétraédrique du domaine cubique [0,1]x[0,1]x[0,1], définition des bords
#==============================================================================================
my_mesh = cdmath.Mesh("ballWithTetrahedra.med")
if( my_mesh.getSpaceDimension()!=3 or my_mesh.getMeshDimension()!=3) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")
if(not my_mesh.isTetrahedral()) :
	raise ValueError("Wrong cell types : mesh is not made of tetrahedra")

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("Number of nodes=", nbNodes)
print("Number of cells=", nbCells)

#Discrétisation du second membre et détermination des noeuds intérieurs
#======================================================================
my_RHSfield = cdmath.Field("RHS_field", cdmath.NODES, my_mesh, 1)
my_ExactSol = cdmath.Field("EXACT_SOL", cdmath.NODES, my_mesh, 1)
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
	z = Ni.z()

	r=sqrt(x*x+y*y+z*z)
	phi=atan2(y,x)
	theta=atan2(y,z*sin(phi))

	my_RHSfield[i]=6*r*sin(theta)**2*cos(phi)+3*(r-1)*cos(phi)#mettre la fonction definie au second membre de l'edp
	my_ExactSol[i]=-(r-1)*r**2*sin(theta)**2*cos(phi)
	if my_mesh.isBorderNode(i): # Détection des noeuds frontière
		boundaryNodes.append(i)
		nbBoundaryNodes=nbBoundaryNodes+1
	else: # Détection des noeuds intérieurs
		interiorNodes.append(i)
		nbInteriorNodes=nbInteriorNodes+1
		maxNbNeighbours= max(1+Ni.getNumberOfEdges(),maxNbNeighbours)

print("Right hand side discretisation done")
print("Number of interior nodes=", nbInteriorNodes)
print("Number of boundary nodes=", nbBoundaryNodes)
print("Maximum number of neighbours per node=", maxNbNeighbours)

# Construction de la matrice de rigidité et du vecteur second membre du système linéaire
#=======================================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours) # warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbInteriorNodes)

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un hexaèdre
GradShapeFunc0=cdmath.Vector(3)
GradShapeFunc1=cdmath.Vector(3)
GradShapeFunc2=cdmath.Vector(3)
GradShapeFunc3=cdmath.Vector(3)

#On parcourt les triangles du domaine
for i in range(nbCells):

	Ci=my_mesh.getCell(i)

	#Contribution à la matrice de rigidité
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
		if boundaryNodes.count(j)==0 : #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
			j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
			#Ajout de la contribution de la cellule triangulaire i au second membre du noeud j 
			RHS[j_int]=Ci.getMeasure()/4*my_RHSfield[j]+RHS[j_int] # intégrale dans le triangle du produit f x fonction de base
			#Contribution de la cellule triangulaire i à la ligne j_int du système linéaire
 			for k in [nodeId0,nodeId1,nodeId2,nodeId3] :
				if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
					k_int=interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
					Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
				#else: si condition limite non nulle au bord, ajouter la contribution du bord au second membre de la cellule j

print("Linear system matrix building done")

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"CG","ILU")#,"ILU" Remplacer CG par CHOLESKY pour solveur direct
SolSyst=LS.solve()
print("Preconditioner used : ", LS.getNameOfPc() )
print("Number of iterations used : ", LS.getNumberOfIter() )
print("Final residual : ", LS.getResidu() )
print("Linear system solved")

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("ResultField", cdmath.NODES, my_mesh, 1)
for j in range(nbInteriorNodes):
   	my_ResultField[interiorNodes[j]]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs

for j in range(nbBoundaryNodes):
    my_ResultField[boundaryNodes[j]]=0;#remplissage des valeurs pour les noeuds frontière (condition limite)
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteElements3DPoisson_BALL_ResultField")

#Postprocessing :
#================
# save 3D picture
resolution=100
VTK_routines.Clip_VTK_data_to_VTK("FiniteElements3DPoisson_BALL_ResultField"+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "FiniteElements3DPoisson_BALL_ResultField"+'_0.vtu',[0.5,0.5,0.5], [-0.5,-0.5,-0.5],resolution )
PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"FiniteElements3DPoisson_BALL_ResultField"+'_0.vtu',"ResultField",'NODES',"Clip_VTK_data_to_VTK_"+"FiniteElements3DPoisson_BALL_ResultField")

# extract and plot diagonal values
curv_abs=np.linspace(0,sqrt(3),resolution+1)
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[-0.577,-0.577,-0.577],[0.577,0.577,0.577], resolution)
plt.plot(curv_abs, diag_data, label= str(nbNodes) + ' nodes 3D mesh')
plt.legend()
plt.xlabel('Position on diagonal line')
plt.ylabel('Value on diagonal line')
plt.title('Plot over diagonal line for finite elements \n for Laplace operator on a 3D ball tetrahedral mesh')
plt.savefig("FiniteElements3DPoisson_BALL_ResultField_"+str(nbNodes) + '_nodes'+"_PlotOverDiagonalLine.png")

print("Numerical solution of 3D Poisson equation on a 3D ball using finite elements done")


#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
max_abs_sol_exacte=max(my_ExactSol.max(),-my_ExactSol.min())
max_sol_num=my_ResultField.max()
min_sol_num=my_ResultField.min()
erreur_abs=0
for i in range(nbNodes) :
    if erreur_abs < abs(my_ExactSol[i] - my_ResultField[i]) :
        erreur_abs = abs(my_ExactSol[i] - my_ResultField[i])

print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)

assert erreur_abs/max_abs_sol_exacte <1.
