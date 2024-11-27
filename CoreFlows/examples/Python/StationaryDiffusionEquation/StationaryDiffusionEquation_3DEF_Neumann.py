#!/usr/bin/env python3
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 3D -\triangle T = f avec conditions aux limites de Neumann \nabla T=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Utilisation de la méthode des éléménts finis P1 avec champs T et f discrétisés aux noeuds d'un maillage tétraédrique
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#================================================================================================================================

import CoreFlows as cf
import cdmath as cm
from math import exp

def StationaryDiffusionEquation_3DEF_BALL_Neumann():
	spaceDim = 3;
	# Prepare for the mesh
	print("Loading mesh " );
	ball_mesh = cm.Mesh("../resources/Ball_1.med")
	if( ball_mesh.getSpaceDimension()!=3 or ball_mesh.getMeshDimension()!=3) :
	    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")
	if(not ball_mesh.isTetrahedral()) :
		raise ValueError("Wrong cell types : mesh is not made of tetrahedra")
	
	nbNodes = ball_mesh.getNumberOfNodes()
	nbCells = ball_mesh.getNumberOfCells()

	print( "Loaded a 3D ball mesh with ", nbNodes, " nodes and ",nbCells, " cells")

	FEComputation=True
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(ball_mesh);

    # set the Neumann boundary condition
	myProblem.setNeumannBoundaryCondition("Boundary")

	#Set the right hand side function
	my_RHSfield = cm.Field("Heat power", cm.NODES, ball_mesh, 1)
	for i in range(ball_mesh.getNumberOfNodes()):
		Ni= ball_mesh.getNode(i)
		x = Ni.x()
		y = Ni.y()
		z = Ni.z()

		my_RHSfield[i]=1000*exp(-(x*x+y*y+z*z))#mettre la fonction definie au second membre de l'edp
	my_RHSfield.writeVTK("HeatPowerField")
	
	myProblem.setHeatPowerField(my_RHSfield)
	myProblem.setLinearSolver(cf.GMRES,cf.ILU);

    # name of result file
	fileName = "3DEF_BALL_Neumann";

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
		####################### Postprocessing #########################
		my_ResultField = myProblem.getOutputTemperatureField()
		my_ResultField.writeMED("TemperatureField_3D_BALL_Neumann")#Pour s'en servir comme donnée initiale
		pass

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    StationaryDiffusionEquation_3DEF_BALL_Neumann()
