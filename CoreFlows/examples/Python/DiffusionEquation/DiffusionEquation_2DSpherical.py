#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab

#===============================================================================================================================
# Name        : Finite Elements simulation of the 2D heat equation -\triangle T = f with Neumann boundary condition
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
#================================================================================================================================


def DiffusionEquation_2DSpherical(FECalculation, fileName):

	""" Description : Test solving the diffusion of the temperature T in a solid (default is Uranium). 
		Equation : Thermal diffusion equation  \rho cp dT/dt-\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T)
		        Heat capacity, density, and conductivity of the solid MUST be defined
		        The solid may be extra refrigerated by a fluid with transfer coefficient using functions setFluidTemperature and setHeatTransfertCoeff
		        The solid may receive some extra heat power due to nuclear fissions using function setHeatSource
	"""
	#Space dimension of the problem
	spaceDim=2
	
    # Mandatory physical values
	solid_specific_heat=300# specific heat capacity
	solid_density=10000# density
	solid_conductivity=5# conductivity

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,solid_density,solid_specific_heat,solid_conductivity);

    # Optional physical values (default value is zero)
	fluidTemperature=573.;#fluid mean temperature
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient
	myProblem.setFluidTemperature(fluidTemperature);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);

	# Definition of field support parameter
	if( FECalculation):
		supportOfField=solverlab.NODES
	else:
		supportOfField=solverlab.CELLS
		
	# Loading heat power field
	heat_power_inputfile="../resources/BoxWithMeshWithTriangularCells";
	heat_power_fieldName="Heat power field";
	heat_power_time_iteration=0
	heat_power_time_sub_iteration=0
	heat_power_meshLevel=0
	
	print("Loading field :", heat_power_fieldName, " in file ", heat_power_inputfile)
	heatPowerField=solverlab.Field(heat_power_inputfile, supportOfField, heat_power_fieldName, heat_power_time_iteration, heat_power_time_sub_iteration, heat_power_meshLevel)
	myProblem.setHeatPowerField(heatPowerField)
	
    # Prepare for the mesh and initial data
	initial_data_inputfile="../resources/BoxWithMeshWithTriangularCells";
	initial_data_fieldName="Temperature";

    #Initial field load
	print("Loading unstructured mesh and initial data", " in file ", initial_data_inputfile )
	initial_data_time_iteration=0
	myProblem.setInitialField(initial_data_inputfile, initial_data_fieldName, initial_data_time_iteration, supportOfField)

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
        # name of result file
        if( FECalculation):
        	fileName = "2DSpherical_FE";
        else:
        	fileName = "2DSpherical_FV";
        DiffusionEquation_2DSpherical(FECalculation, fileName)
    else :
        raise ValueError("DiffusionEquation_2DSpherical : missing one argument")
