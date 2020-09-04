#!/usr/bin/env python
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson 2D -\triangle T = f sur un carré avec conditions aux limites de Dirichlet T=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des volumes finis avec champs T et f discrétisés aux cellules d'un maillage de triangules équilatéraux
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution numérique avec la solution exacte T=sin(pi*x)*sin(pi*y)
#================================================================================================================================

import CoreFlows as cf
import cdmath as cm
from math import sin, pi
import os

def StationaryDiffusionEquation_2DFV_EquilateralTriangles():
	spaceDim = 2;
	# Prepare for the mesh
	print("Loading mesh " );
	M=cm.Mesh(os.environ['CDMATH_INSTALL']+'/share/meshes/2DEquilateralTriangles/squareWithEquilateralTriangles20.med')#Equilateral triangular mesh
	
	print( "Loaded 2D equilateral triangle mesh with ", M.getNumberOfCells(), " cells")

	# set the limit field 
	boundaryFaces = M.getBoundaryFaceIds()
	boundaryValues = {}
	print("Setting Dirichlet boundary values")
	for i in boundaryFaces :
		Fi=M.getFace(i)
		x=Fi.x()
		y=Fi.y()
		boundaryValues[i] = sin(pi*x)*sin(pi*y)
		
	FEComputation=False
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(M);
	myProblem.setDirichletValues(cf.MapIntDouble(boundaryValues))
	
	#Set the right hand side function
	my_RHSfield = cm.Field("RHS_field", cm.CELLS, M, 1)
	for i in range(M.getNumberOfCells()):
		Ci= M.getCell(i)
		x = Ci.x()
		y = Ci.y()

		my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
	
	myProblem.setHeatPowerField(my_RHSfield)
	myProblem.setLinearSolver(cf.GMRES,cf.ILU);

	# name of result file
	fileName = "StationnaryDiffusion_2DFV_StructuredTriangles";

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
		for i in range(M.getNumberOfCells()) :
			if erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) :
				erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i])
		
		print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
		print ("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
		
		assert erreur_abs/max_abs_sol_exacte <1.
        pass

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	StationaryDiffusionEquation_2DFV_EquilateralTriangles()
