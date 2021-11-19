#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab

#===============================================================================================================================
# Name        : Finite Elements simulation of the 2D heat equation -\triangle T = f with Neumann boundary condition
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
#================================================================================================================================


def DiffusionEquation_2DSpherical(FECalculation):

	""" Description : Test solving the diffusion of the temperature T in a solid (default is Uranium). 
		Equation : Thermal diffusion equation  \rho cp dT/dt-\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T)
		        Heat capacity, density, and conductivity of the solid MUST be defined
		        The solid may be extra refrigerated by a fluid with transfer coefficient using functions setFluidTemperature and setHeatTransfertCoeff
		        The solid may receive some extra heat power due to nuclear fissions using function setHeatSource
	"""
	
    # Prepare for the mesh and initial data
	inputfile="../resources/BoxWithMeshWithTriangularCells";
	fieldName="Temperature";
	spaceDim=2
	
    # Mandatory physical values
	specific_heat=300# specific heat capacity
	density=10000# density
	conductivity=5# conductivity

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,density,specific_heat,conductivity);

    # Optional physical values (default value is zero)
	fluidTemperature=573.;#fluid mean temperature
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient
	constant_heat=1e5;#heat power ddensity
	myProblem.setFluidTemperature(fluidTemperature);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	myProblem.setHeatSource(constant_heat);

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
	myProblem.setLinearSolver(solverlab.GMRES,solverlab.ILU);

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
