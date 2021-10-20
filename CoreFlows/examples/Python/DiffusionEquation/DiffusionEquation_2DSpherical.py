#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab

#===============================================================================================================================
# Name        : Finite Elements simulation of the 2D heat equation -\triangle T = f with Neumann boundary condition
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : Test solving the diffusion of the temperature T in a solid (default is Uranium). 
#		        \rho cp dT/dt-\lambda\Delta T=\Phi(T) + \lambda_{sf} (T_{fluid}-T)
#		        Heat capacity cp, density \rho, and conductivity \lambda MUST be defined
#		        The solid may be extra refrigerated by a fluid with transfer coefficient \lambda_{sf} (functions setFluidTemperature and setHeatTransfertCoeff)
#		        The solid may receive some extra heat power \Phi due to nuclear fissions (function setHeatSource)
#================================================================================================================================


def DiffusionEquation_2DSpherical(FECalculation):

    # Prepare for the mesh
	inputfile="../resources/BoxWithMeshWithTriangularCells";
	fieldName="Temperature";
	spaceDim=2
	
    # Mandatory physical values
	cp_ur=300# heat capacity
	rho_ur=10000# density
	lambda_ur=5# conductivity

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,rho_ur,cp_ur,lambda_ur);

    #Optional physical values
	fluidTemp=573.;#fluid mean temperature
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient
	phi=1e5;#heat power ddensity
	myProblem.setFluidTemperature(fluidTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	myProblem.setHeatSource(phi);

    #Initial field load
	time_iteration=0
	print("Loading unstructured mesh and initial data" )
	if( FECalculation):
		myProblem.setInitialField(inputfile,fieldName,time_iteration, solverlab.NODES)
	else:
		myProblem.setInitialField(inputfile,fieldName,time_iteration, solverlab.CELLS)

    # the boundary conditions :
	if( FECalculation):
		boundaryNodeGroupNames=myProblem.getMesh().getNameOfNodeGroups()
		print(len(boundaryNodeGroupNames), " Boundary Node Group detected : ", boundaryNodeGroupNames)
	else:
		boundaryFaceGroupNames=myProblem.getMesh().getNameOfFaceGroups()
		print(len(boundaryFaceGroupNames), " Boundary Face Group detected : ", boundaryFaceGroupNames)

	myProblem.setNeumannBoundaryCondition("GAUCHE");
	myProblem.setNeumannBoundaryCondition("DROITE");
	myProblem.setNeumannBoundaryCondition("HAUT");
	myProblem.setNeumannBoundaryCondition("BAS");

    # set the numerical method
	myProblem.setTimeScheme( solverlab.Explicit);
	myProblem.setLinearSolver(solverlab.GMRES,solverlab.ILU,True);

    # name of result file
	if( FECalculation):
		fileName = "2DSpherical_FE";
	else:
		fileName = "2DSpherical_FV";

    # computation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 100000000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

    # evolution
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        FECalculation = bool(int(sys.argv[1]))
        DiffusionEquation_2DSpherical(FECalculation)
    else :
        raise ValueError("DiffusionEquation_2DSpherical : missing one argument")
