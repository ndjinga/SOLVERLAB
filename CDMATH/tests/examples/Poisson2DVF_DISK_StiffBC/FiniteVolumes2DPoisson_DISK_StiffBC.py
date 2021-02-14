# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson -\triangle u = 0 sur un disque avec conditions aux limites de Dirichlet discontinues
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage quelconque
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=atan(2*x/(x**2+y**2-1))=atan(2 r cos(theta)/(r*2-1))
#================================================================================================================================

import cdmath
from math import atan, pi, sqrt
from numpy import linspace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import PV_routines
import VTK_routines
import sys

if len(sys.argv) >1 :#non rectangular mesh
    meshName=sys.argv[1]
    my_mesh = cdmath.Mesh(sys.argv[1])
else :
    raise ValueError("Give an input mesh of the disk")
    
nbCells = my_mesh.getNumberOfCells()

if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")

print( "Mesh loading done" )
print( "Number of cells  = ", nbCells )

#Discrétisation du second membre et extraction du nb max de voisins d'une cellule
#================================================================================
my_ExactSol = cdmath.Field("Exact_field", cdmath.CELLS, my_mesh, 1)
eps=1e-6#For coarse meshes

maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
#parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
for i in range(nbCells): 
    Ci = my_mesh.getCell(i)
    x = Ci.x()
    y = Ci.y()

    #Robust calculation of atan(2x/(x**2+y**2-1)
    if x**2+y**2-1 > eps :
        print("!!! Warning Mesh ", meshName," !!! Cell is not in the unit disk."," eps=",eps, ", x**2+y**2-1=",x**2+y**2 - 1)
        #raise ValueError("x**2+y**2 > 1 !!! Domain should be the unit disk.")
    if x**2+y**2-1 < -eps :
        my_ExactSol[i] = atan(2*x/(x**2+y**2-1))
    elif x>0 : #x**2+y**2-1>=0
        my_ExactSol[i] = -pi/2
    elif x<0 : #x**2+y**2-1>=0
        my_ExactSol[i] =  pi/2
    else : #x=0
        my_ExactSol[i] = 0
        
    # compute maximum number of neighbours
    maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)

print("Right hand side discretisation done")
print( "Maximum number of neighbours = ", maxNbNeighbours)

# Construction de la matrice et du vecteur second membre du système linéaire
#===========================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbCells)

#Parcours des cellules du domaine
for i in range(nbCells):
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
            #For the particular case where the mesh boundary does not coincide with the domain boundary
            x=Fj.getBarryCenter().x()
            y=Fj.getBarryCenter().y()
            if x**2+y**2-1 > eps :
                print("!!! Warning Mesh ", meshName," !!! Face is not in the unit disk.",", eps=",eps, ", x**2+y**2-1=",x**2+y**2 - 1)
                #raise ValueError("!!! Domain should be the unit disk.")
            if x**2+y**2-1 < -eps :
                RHS[i]+= coeff*atan(2*x/(x**2+y**2-1))
            elif x>0 : #x**2+y**2-1>=0
                RHS[i]+= coeff*(-pi/2)
            elif x<0 : #x**2+y**2-1>=0
                RHS[i]+= coeff*pi/2
            else : #x=0
                RHS[i]+=  0
        Rigidite.addValue(i,i,coeff) # terme diagonal

print("Linear system matrix building done")

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"GMRES","ILU")
SolSyst=LS.solve()

print( "Preconditioner used : ", LS.getNameOfPc() )
print( "Number of iterations used : ", LS.getNumberOfIter() )
print( "Final residual : ", LS.getResidu() )
print( "Linear system solved" )

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("ResultField", cdmath.CELLS, my_mesh, 1)
for i in range(nbCells):
    my_ResultField[i]=SolSyst[i];
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteVolumes2DPoisson_DISK_ResultField")

#Postprocessing : 
#===============
# save 2D picture
PV_routines.Save_PV_data_to_picture_file("FiniteVolumes2DPoisson_DISK_ResultField"+'_0.vtu',"ResultField",'CELLS',"FiniteVolumes2DPoisson_DISK_ResultField")

# extract and plot diagonal values
resolution=100
curv_abs=linspace(0,sqrt(2),resolution+1)
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,-1,0],[0,1,0], resolution)
plt.legend()
plt.xlabel('Position on diagonal line')
plt.ylabel('Value on diagonal line')
if len(sys.argv) >1 :
    plt.title('Plot over diagonal line for finite Volumes \n for Laplace operator on a 2D disk mesh '+my_mesh.getName())
    plt.plot(curv_abs, diag_data, label= str(nbCells)+ ' cells mesh')
    plt.savefig("FiniteVolumes2DPoisson_DISK_ResultField_"+str(nbCells)+ '_cells'+"_PlotOverDiagonalLine.png")
else :   
    plt.title('Plot over diagonal line for finite Volumes \n for Laplace operator on a 2D disk with a rectangular grid')
    plt.plot(curv_abs, diag_data, label= str(nx) +'x'+str(ny)+ ' cells mesh')
    plt.savefig("FiniteVolumes2DPoisson_DISK_ResultField_"+str(nx) +'x'+str(ny)+ '_cells'+"_PlotOverDiagonalLine.png")

print("Numerical solution of 2D Poisson equation on a disk using finite volumes done")

#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
l2_norm_sol_exacte=my_ExactSol.normL2()[0]
l2_error = (my_ExactSol - my_ResultField).normL2()[0]

print("L2 absolute error = norm( exact solution - numerical solution ) = ",l2_error )
print("L2 relative error = norm( exact solution - numerical solution )/norm( exact solution ) = ",l2_error/l2_norm_sol_exacte ) 
print("Maximum numerical solution = ", my_ResultField.max(), " Minimum numerical solution = ", my_ResultField.min() )
print("Maximum exact solution = ", my_ExactSol.max(), " Minimum exact solution = ", my_ExactSol.min() )

assert l2_error/l2_norm_sol_exacte <1.
