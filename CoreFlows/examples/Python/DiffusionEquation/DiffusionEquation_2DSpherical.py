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
	solid_specific_heat=300# specific heat capacity, default value 300
	solid_density=10000# density, default value 10000
	solid_conductivity=5# conductivity, default value 5

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,solid_density,solid_specific_heat,solid_conductivity);

	# Definition of field support parameter
	if( FECalculation):
		supportOfField=solverlab.NODES
	else:
		supportOfField=solverlab.CELLS	
	
    # Set the mesh and initial data
	initial_data_inputfile="../resources/BoxWithMeshWithTriangularCells";
	initial_data_fieldName="Temperature";
	print("Loading unstructured mesh and initial data", " in file ", initial_data_inputfile )
	initial_data_time_iteration=0# default value is 0
	myProblem.setInitialField(initial_data_inputfile, initial_data_fieldName, initial_data_time_iteration, supportOfField)

    #### Optional physical values (default value is zero) : fluid temperature field, heat transfert coefficient, heat power field 
	# Loading and setting fluid temperature field
	fluid_temperature_inputfile="../resources/BoxWithMeshWithTriangularCells";
	fluid_temperature_fieldName="Fluid temperature field";
	fluid_temperature_time_iteration=0# default value is 0
	fluid_temperature_time_sub_iteration=0# default value is 0
	fluid_temperature_meshLevel=0# default value is 0
	print("Loading field :", fluid_temperature_fieldName, " in file ", fluid_temperature_inputfile)
	fluidTemperatureField=solverlab.Field(fluid_temperature_inputfile, supportOfField, fluid_temperature_fieldName, fluid_temperature_time_iteration, fluid_temperature_time_sub_iteration, fluid_temperature_meshLevel)
	myProblem.setFluidTemperatureField(fluidTemperatureField)
	# Setting heat transfert coefficient
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient, default value is 0
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	# Loading heat power field
	heat_power_inputfile="../resources/BoxWithMeshWithTriangularCells";
	heat_power_fieldName="Heat power field";
	heat_power_time_iteration=0# default value is 0
	heat_power_time_sub_iteration=0# default value is 0
	heat_power_meshLevel=0# default value is 0
	print("Loading field :", heat_power_fieldName, " in file ", heat_power_inputfile)
	heatPowerField=solverlab.Field(heat_power_inputfile, supportOfField, heat_power_fieldName, heat_power_time_iteration, heat_power_time_sub_iteration, heat_power_meshLevel)
	myProblem.setHeatPowerField(heatPowerField)

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
	MaxNbOfTimeStep = 3 ;# default value is 10
	freqSave = 1;# default value is 1
	cfl = 0.95;# default value is 1
	maxTime = 100000000;# default value is 10
	precision = 1e-6;# default value is 1e-6

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
        	fileName = "2DSpherical_FE";# default value is ""
        else:
        	fileName = "2DSpherical_FV";# default value is ""
        DiffusionEquation_2DSpherical(FECalculation, fileName)
    else :
        raise ValueError("DiffusionEquation_2DSpherical : missing one argument")
