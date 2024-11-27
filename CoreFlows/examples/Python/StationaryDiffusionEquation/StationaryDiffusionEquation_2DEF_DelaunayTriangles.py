#!/usr/bin/env python
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 2D -\triangle T = f sur un carré avec conditions aux limites de Dirichlet variables
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2024
# Description : Utilisation de la méthode des éléments finis avec champs T et f discrétisés aux noeuds d'un maillage de triangles de Delaunay
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution numérique avec la solution exacte T=cos(pi*x)*cos(pi*y)
#================================================================================================================================

import CoreFlows as cf
import cdmath as cm
from math import cos, pi
import os

def StationaryDiffusionEquation_2DFE_EquilateralTriangles_cosinus():
	spaceDim = 2;
	# Prepare for the mesh
	print("Loading mesh " );
	M=cm.Mesh('../resources/squareWithEquilateralTriangles20.med')#Delaunay triangular mesh
	
	print( "Loaded 2D equilateral triangle mesh with ", M.getNumberOfNodes(), " nodes")

	# set the limit field 
	boundaryNodes = M.getBoundaryNodeIds()
	boundaryValues = {}
	print("Setting Dirichlet boundary values")
	for i in boundaryNodes :
		Ni=M.getNode(i)
		x=Ni.x()
		y=Ni.y()
		boundaryValues[i] = cos(pi*x)*cos(pi*y)
		
	FEComputation=True
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(M);
	myProblem.setDirichletValues(cf.MapIntDouble(boundaryValues))
	
	#Set the right hand side function
	my_RHSfield = cm.Field("RHS_field", cm.NODES, M, 1)
	for i in range(M.getNumberOfNodes()):
		Ni= M.getNode(i)
		x = Ni.x()
		y = Ni.y()

		my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y)#mettre la fonction definie au second membre de l'edp
	
	myProblem.setHeatPowerField(my_RHSfield)
	myProblem.setLinearSolver(cf.GMRES,cf.ILU);

	# name of result file
	fileName = "StationnaryDiffusion_2DFE_DelaunayTriangles";

	# computation parameters
	myProblem.setFileName(fileName);

	# Run the computation
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.solveStationaryProblem();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );
		####################### Postprocessing #########################
		my_ResultField = myProblem.getOutputTemperatureField()
		#The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi
		max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(2*pi*pi)
		max_sol_num=my_ResultField.max()
		min_sol_num=my_ResultField.min()
		erreur_abs=0
		for i in range(M.getNumberOfNodes()) :
			if erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) :
				erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i])
		
		print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
		print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
		
		assert erreur_abs/max_abs_sol_exacte <1.
		pass

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	StationaryDiffusionEquation_2DFE_EquilateralTriangles_cosinus()
