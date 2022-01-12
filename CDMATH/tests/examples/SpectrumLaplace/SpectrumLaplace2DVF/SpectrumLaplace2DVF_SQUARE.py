# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Calcul VF du spectre de l'opérateur de Laplace 2D -\triangle avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2020
# Description : Utilisation de la méthode des volumes finis P1 avec champs discrétisés aux cellules d'un maillage quelconque
#		        Création et sauvegarde des champs résultant en utilisant la librairie CDMATH
#================================================================================================================================

import cdmath
import sys

if len(sys.argv) >2 :#non rectangular mesh
    my_mesh = cdmath.Mesh(sys.argv[1])
    mesh_name=sys.argv[2]
else :   #rectangular mesh
# Création d'un maillage cartésien du domaine carré [0,1]x[0,1], définition des bords
#====================================================================================
    xmin=0
    xmax=1
    ymin=0
    ymax=1
    
    nx=15
    ny=15
    
    my_mesh = cdmath.Mesh(xmin,xmax,nx,ymin,ymax,ny)
    mesh_name="RegularGrid"

    eps=1e-6
    my_mesh.setGroupAtPlan(0,0,eps,"DirichletBorder")#Bord GAUCHE
    my_mesh.setGroupAtPlan(1,0,eps,"DirichletBorder")#Bord DROIT
    my_mesh.setGroupAtPlan(0,1,eps,"DirichletBorder")#Bord BAS
    my_mesh.setGroupAtPlan(1,1,eps,"DirichletBorder")#Bord HAUT

nbCells = my_mesh.getNumberOfCells()

if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")

print("Mesh loading done")
print("Number of cells  = ", nbCells)

#Détermination du nb max de voisins d'une cellule
#================================================
maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite volumes rigidity matrix
#parcours des cellules pour extraction du nb max de voisins d'une cellule
for i in range(nbCells): 
    Ci=my_mesh.getCell(i)
    # compute maximum number of neighbours
    maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)

print("Max nb of neighbours = ", maxNbNeighbours)

# Construction de la matrice et du vecteur second membre du système linéaire
#===========================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
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
        Rigidite.addValue(i,i,coeff) # terme diagonal

print("Stiffness matrix construction done")
Rigidite.viewMatrix(True, 0, "RigidityMatrix_FiniteVolumesOn"+mesh_name+"Laplace")
Rigidite.plotEigenvalues("FiniteVolumesOn"+mesh_name+"Laplace")

# Conditionnement de la matrice de rigidité
#=================================
cond = Rigidite.getConditionNumber()
print("Condition number is ",cond)

# Spectre de la matrice de rigidité
#==================================
nev=10
d=Rigidite.getEigenvectorsDataArrayDouble(nev)
my_eigenfield = cdmath.Field("Eigenvectors field", cdmath.CELLS, my_mesh, nev)
my_eigenfield.setFieldByDataArrayDouble(d)

# Sauvegarde du champ résultat
#===========================
my_eigenfield.writeVTK("spectrumFiniteVolumesOn"+mesh_name+"Laplace")
