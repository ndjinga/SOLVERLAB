#!/usr/bin/env python3
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Laplace-Beltrami -\triangle u = f sur un tore
# Author      : Michael Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage triangulaire
#               Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Référence : M. A. Olshanskii, A. Reusken, and J. Grande. A finite element method for elliptic equations
#                           on surfaces. SIAM J. Num. Anal., 47, p. 3355
#               Résolution d'un système linéaire à matrice singulière : les vecteurs constants sont dans le noyau
#================================================================================================================================

import cdmath
from math import sin, cos, atan2, sqrt
import PV_routines
import VTK_routines
import paraview.simple as pvs

#Chargement du maillage triangulaire du tore
#=======================================================================================
my_mesh = cdmath.Mesh("meshTorus.med")
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

# Torus radii (calculation  will fail if the mesh is not correct)
R=1 #Grand rayon
r=0.6 #Petit rayon

#Discrétisation du second membre, de la solution exacte et détermination des noeuds intérieurs
#======================================================================
my_RHSfield = cdmath.Field("RHS field", cdmath.NODES, my_mesh, 1)
exactSolField = cdmath.Field("Exact solution field", cdmath.NODES, my_mesh, 1)

maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix

#parcours des noeuds pour discrétisation du second membre et extraction du nb max voisins d'un noeud
for i in range(nbNodes):
	Ni=my_mesh.getNode(i)
	x = Ni.x()
	y = Ni.y()
	z = Ni.z()
      
	theta=atan2(z,sqrt(x*x+y*y)-R)
	phi=atan2(y,x)

	exactSolField[i] = sin(3*phi)*cos(3*theta+ phi) # for the exact solution we use the funtion given in the article of Olshanskii, Reusken 2009, page 19
	my_RHSfield[i] = 9*sin(3*phi)*cos(3*theta+ phi)/(r*r) + (10*sin(3*phi)*cos(3*theta+ phi) + 6*cos(3*phi)*sin(3*theta+ phi))/((R+r*cos(theta))*(R+r*cos(theta))) - 3*sin(theta)*sin(3*phi)*sin(3*theta+ phi)/(r*(R+r*cos(theta))) #for the right hand side we use the function given in the article of Olshanskii, Reusken 2009, page 19
	if my_mesh.isBorderNode(i): # Détection des noeuds frontière
		raise ValueError("Mesh should not contain borders")
	else:
		maxNbNeighbours = max(1+Ni.getNumberOfCells(),maxNbNeighbours) #true only for planar cells, otherwise use function Ni.getNumberOfEdges()

print("Right hand side discretisation done")
print("Max nb of neighbours=", maxNbNeighbours)
print("Integral of the RHS", my_RHSfield.integral(0))

# Construction de la matrice de rigidité et du vecteur second membre du système linéaire
#=======================================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbNodes,nbNodes,maxNbNeighbours)
RHS=cdmath.Vector(nbNodes)

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un triangle
GradShapeFunc0=cdmath.Vector(3)
GradShapeFunc1=cdmath.Vector(3)
GradShapeFunc2=cdmath.Vector(3)

normalFace0=cdmath.Vector(3)
normalFace1=cdmath.Vector(3)

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
		#Ajout de la contribution de la cellule triangulaire i au second membre du noeud j 
		RHS[j]=Ci.getMeasure()/3*my_RHSfield[j]+RHS[j] # intégrale dans le triangle du produit f x fonction de base
		#Contribution de la cellule triangulaire i à la ligne j du système linéaire
		for k in [nodeId0,nodeId1,nodeId2] : 
			Rigidite.addValue(j,k,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())

print("Linear system matrix building done")

# Conditionnement de la matrice de rigidité
#=================================
cond = Rigidite.getConditionNumber(True)
print("Condition number is ",cond)

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"GMRES","ILU")#Remplacer CG par CHOLESKY pour solveur direct
LS.setMatrixIsSingular()#En raison de l'absence de bord
SolSyst=LS.solve()
print( "Preconditioner used : ", LS.getNameOfPc() )
print( "Number of iterations used : ", LS.getNumberOfIter() )
print( "Final residual : ", LS.getResidu() )
print("Linear system solved")

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("Numerical result field", cdmath.NODES, my_mesh, 1)
for j in range(nbNodes):
    my_ResultField[j]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteElementsOnTorusPoisson")

print("Integral of the numerical solution", my_ResultField.integral(0))
print("Numerical solution of Poisson equation on a torus using finite elements done")

#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
max_sol_exacte=exactSolField.normMax()[0]
erreur_max=(exactSolField - my_ResultField).normMax()[0]
max_sol_num=my_ResultField.max()
min_sol_num=my_ResultField.min()

print("Relative error =  max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_max/max_sol_exacte)
print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
print("Maximum exact solution = ", exactSolField.max(), " Minimum exact solution = ", exactSolField.min())

#Postprocessing : 
#================
# Save 3D picture
PV_routines.Save_PV_data_to_picture_file("FiniteElementsOnTorusPoisson"+'_0.vtu',"Numerical result field",'NODES',"FiniteElementsOnTorusPoisson")
resolution=100
VTK_routines.Clip_VTK_data_to_VTK("FiniteElementsOnTorusPoisson"+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "FiniteElementsOnTorusPoisson"+'_0.vtu',[0.25,0.25,0.25], [-0.5,-0.5,-0.5],resolution )
PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"FiniteElementsOnTorusPoisson"+'_0.vtu',"Numerical result field",'NODES',"Clip_VTK_data_to_VTK_"+"FiniteElementsOnTorusPoisson")

# Plot  over slice circle
finiteElementsOnTorus_0vtu = pvs.XMLUnstructuredGridReader(FileName=["FiniteElementsOnTorusPoisson"+'_0.vtu'])
slice1 = pvs.Slice(Input=finiteElementsOnTorus_0vtu)
slice1.SliceType.Normal = [0.5, 0.5, 0.5]
renderView1 = pvs.GetActiveViewOrCreate('RenderView')
finiteElementsOnTorus_0vtuDisplay = pvs.Show(finiteElementsOnTorus_0vtu, renderView1)
pvs.ColorBy(finiteElementsOnTorus_0vtuDisplay, ('POINTS', 'Numerical result field'))
slice1Display = pvs.Show(slice1, renderView1)
pvs.SaveScreenshot("./FiniteElementsOnTorusPoisson"+"_Slice"+'.png', magnification=1, quality=100, view=renderView1)
plotOnSortedLines1 = pvs.PlotOnSortedLines(Input=slice1)
lineChartView2 = pvs.CreateView('XYChartView')
plotOnSortedLines1Display = pvs.Show(plotOnSortedLines1, lineChartView2)
plotOnSortedLines1Display.UseIndexForXAxis = 0
plotOnSortedLines1Display.XArrayName = 'arc_length'
plotOnSortedLines1Display.SeriesVisibility = ['Numerical result field (1)']
pvs.SaveScreenshot("./FiniteElementsOnTorusPoisson"+"_PlotOnSortedLine_"+'.png', magnification=1, quality=100, view=lineChartView2)
pvs.Delete(lineChartView2)

assert erreur_max/max_sol_exacte <1.
