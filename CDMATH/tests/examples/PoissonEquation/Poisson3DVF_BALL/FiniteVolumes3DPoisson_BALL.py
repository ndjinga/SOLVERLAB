#!/usr/bin/env python3
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson -\triangle u = f sur la boule unité  avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage quelconque
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=-(r-1)*r**2*sin(theta)**2*cos(phi)
#================================================================================================================================

import cdmath
from math import sin, cos, pi, sqrt, atan2
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import PV_routines
import VTK_routines

import sys

#Chargement du maillage tétraédrique du domaine cubique [0,1]x[0,1]x[0,1], définition des bords
#==============================================================================================
my_mesh = cdmath.Mesh("ballWithTetrahedra.med")
if( my_mesh.getSpaceDimension()!=3 or my_mesh.getMeshDimension()!=3) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")

eps=1e-6
my_mesh.setGroupAtPlan(0,0,eps,"DirichletBorder")#Bord GAUCHE
my_mesh.setGroupAtPlan(1,0,eps,"DirichletBorder")#Bord DROIT
my_mesh.setGroupAtPlan(0,1,eps,"DirichletBorder")#Bord BAS
my_mesh.setGroupAtPlan(1,1,eps,"DirichletBorder")#Bord HAUT
my_mesh.setGroupAtPlan(0,2,eps,"DirichletBorder")#Bord AVANT 
my_mesh.setGroupAtPlan(1,2,eps,"DirichletBorder")#Bord ARRIERE 

nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("Number of cells ", nbCells)

#Discrétisation du second membre et extraction du nb max de voisins d'une cellule
#================================================================================
my_RHSfield = cdmath.Field("RHS_field", cdmath.CELLS, my_mesh, 1)
my_ExactSol = cdmath.Field("EXACT_SOL", cdmath.CELLS, my_mesh, 1)

maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
#parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
for i in range(nbCells): 
	Ci = my_mesh.getCell(i)
	x = Ci.x()
	y = Ci.y()
	z = Ci.z() 

	r=sqrt(x*x+y*y+z*z)
	phi=atan2(y,x)
	theta=atan2(y,z*sin(phi))

	my_RHSfield[i]=6*r*sin(theta)**2*cos(phi)+3*(r-1)*cos(phi)#mettre la fonction definie au second membre de l'edp
	my_ExactSol[i]=-(r-1)*r**2*sin(theta)**2*cos(phi)
	# compute maximum number of neighbours
	maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)

print("Right hand side discretisation done")
print("Maximum number of neighbours=", maxNbNeighbours)

# Construction de la matrice et du vecteur second membre du système linéaire
#===========================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbCells)
#Parcours des cellules du domaine
for i in range(nbCells):
	RHS[i]=my_RHSfield[i] #la valeur moyenne du second membre f dans la cellule i
	Ci=my_mesh.getCell(i)
	for j in range(Ci.getNumberOfFaces()):# parcours des faces voisinnes
		Fj=my_mesh.getFace(Ci.getFaceId(j))
		if not Fj.isBorder():
			k=Fj.getCellId(0)
			if k==i :
				k=Fj.getCellId(1)
			Ck=my_mesh.getCell(k)
			distance=Ci.getBarryCenter().distance(Ck.getBarryCenter())
			coeff=Fj.getMeasure()/Ci.getMeasure()/distance
			Rigidite.addValue(i,k,-coeff) # terme extradiagonal
		else:
			coeff=Fj.getMeasure()/Ci.getMeasure()/Ci.getBarryCenter().distance(Fj.getBarryCenter())
		Rigidite.addValue(i,i,coeff) # terme diagonal

print("Linear system matrix building done")

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"GMRES","ILU")
SolSyst=LS.solve()

print( "Preconditioner used : ", LS.getNameOfPc() )
print( "Number of iterations used : ", LS.getNumberOfIter() )
print( "Final residual : ", LS.getResidu() )
print("Linear system solved")

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("ResultField", cdmath.CELLS, my_mesh, 1)
for i in range(nbCells):
    my_ResultField[i]=SolSyst[i];
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteVolumes3DPoisson_BALL_ResultField")

#Postprocessing 
#==============
# save 3D picture
resolution=100
VTK_routines.Clip_VTK_data_to_VTK("FiniteVolumes3DPoisson_BALL_ResultField"+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "FiniteVolumes3DPoisson_BALL_ResultField"+'_0.vtu',[0.5,0.5,0.5], [-0.5,-0.5,-0.5],resolution )
PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"FiniteVolumes3DPoisson_BALL_ResultField"+'_0.vtu',"ResultField",'CELLS',"Clip_VTK_data_to_VTK_"+"FiniteVolumes3DPoisson_BALL_ResultField")

# extract and plot diagonal values
resolution=100
curv_abs=np.linspace(0,sqrt(3),resolution+1)
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,1,0],[1,0,0], resolution)
plt.legend()
plt.xlabel('Position on diagonal line')
plt.ylabel('Value on diagonal line')
if len(sys.argv) >1 :
    plt.title('Plot over diagonal line for finite volumes \n for Laplace operator on a 3D BALL with  mesh '+my_mesh.getName())
    plt.plot(curv_abs, diag_data, label= str(nbCells)+ ' cells mesh')
    plt.savefig("FiniteVolumes3DPoisson_BALL_ResultField_"+str(nbCells)+ '_cells'+"_PlotOverDiagonalLine.png")
else :   
    plt.title('Plot over diagonal line for finite volumes \n for Laplace operator on a 3D BALL with a rectangular grid')
    plt.plot(curv_abs, diag_data, label= str(nbCells)+ ' cells mesh')
    plt.savefig("FiniteVolumes3DPoisson_BALL_ResultField_"+str(nbCells)+ '_cells'+"_PlotOverDiagonalLine.png")

print("Numerical solution of 3D poisson equation on a 3D ball using finite volumes done")

#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
max_abs_sol_exacte=max(my_ExactSol.max(),-my_ExactSol.min())
max_sol_num=my_ResultField.max()
min_sol_num=my_ResultField.min()
erreur_abs=0
for i in range(nbCells) :
    if erreur_abs < abs(my_ExactSol[i] - my_ResultField[i]) :
        erreur_abs = abs(my_ExactSol[i] - my_ResultField[i])

print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)

assert erreur_abs/max_abs_sol_exacte <1.
