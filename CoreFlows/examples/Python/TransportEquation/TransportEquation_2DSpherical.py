#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab

#===============================================================================================================================
# Name        : Simulation of a 2D transport equation 
# Description : dh/dt + \div h \vec v  = \Phi + \lambda_{sf} (T_{solid} - T_{fluid}) 
#               Neumann or Dirichlet boundary conditions
#               Finite volumes or finite elements
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
#================================================================================================================================


def TransportEquation_2DSpherical(  fileName):

	""" Description : Test solving the transport of the enthalpy h in a fluid (water or steam). 
		Equation : Linear transport equation  dh/dt + \div h \vec v  = \lambda_{sf} (T_{solid} - T_{fluid}) 
					Phase (water or steam), pressure of the fluid MUST be defined
		        The fluid may be heated by a solid with temperature T_{solid} and transfer coefficient \lambda_{sf} using functions setRodTemperature and setHeatTransfertCoeff
		        The solid may receive some extra heat power due to nuclear fissions using function setHeatSource
	"""
	#Space dimension of the problem
	spaceDim=2
	
    # Mandatory physical values
	transport_velocity=[5,5]# Constant velocity vector of dimension spaceDim
	fluid_phase=solverlab.Water# ou bien solverlab.Air
	pressure_regime=solverlab.around155bars600K# ou bien solverlab.around1bar300KTransport
	myProblem = solverlab.TransportEquation(fluid_phase,pressure_regime,transport_velocity);

    # Set the mesh and initial data
	initial_data_inputfile="../resources/meshSquare";
	initial_data_fieldName="Fluid temperature";
	print("Loading unstructured mesh and initial data", " in file ", initial_data_inputfile )
	initial_data_time_iteration=0# default value is 0
	initial_data_time_sub_iteration=0# default value is 0
	initial_data_time_meshLevel=0# default value is 0
	myProblem.setInitialField(initial_data_inputfile, initial_data_fieldName, initial_data_time_iteration, initial_data_time_sub_iteration, initial_data_time_meshLevel, solverlab.CELLS)

    #### Optional physical values (default value is zero) : fluid temperature field, heat transfert coefficient, heat power field 
	# Loading and setting fluid temperature field
	solid_temperature_inputfile="../resources/meshSquare";
	solid_temperature_fieldName="Fluid temperature";
	solid_temperature_time_iteration=0# default value is 0
	solid_temperature_time_sub_iteration=0# default value is 0
	solid_temperature_meshLevel=0# default value is 0
	print("Loading field :", solid_temperature_fieldName, " in file ", solid_temperature_inputfile)
	solidTemperatureField=solverlab.Field(solid_temperature_inputfile, solverlab.CELLS, solid_temperature_fieldName, solid_temperature_time_iteration, solid_temperature_time_sub_iteration, solid_temperature_meshLevel)
	myProblem.setRodTemperatureField(solidTemperatureField)
	# Setting heat transfert coefficient
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient, default value is 0
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	# Loading heat power field
	heat_power_inputfile="../resources/meshSquare";
	heat_power_fieldName="Heat power";
	heat_power_time_iteration=0# default value is 0
	heat_power_time_sub_iteration=0# default value is 0
	heat_power_meshLevel=0# default value is 0
	print("Loading field :", heat_power_fieldName, " in file ", heat_power_inputfile)
	heatPowerField=solverlab.Field(heat_power_inputfile, solverlab.CELLS, heat_power_fieldName, heat_power_time_iteration, heat_power_time_sub_iteration, heat_power_meshLevel)
	myProblem.setHeatPowerField(heatPowerField)

    # the boundary conditions :
	boundaryGroupNames=myProblem.getMesh().getNameOfFaceGroups()
	print(len(boundaryGroupNames), " Boundary Face Group detected : ", boundaryGroupNames)

	# for each boundary we load the boundary field (replace by a loop over the boundaries)
	boundary1_type=solverlab.NeumannTransport
	boundary1_inputfile="../resources/meshSquare";
	boundary1_fieldName="Left temperature";
	boundary1_time_iteration=0# default value is 0
	boundary1_time_sub_iteration=0# default value is 0
	boundary1_meshLevel=0# default value is 0
	print("Boundary ", boundaryGroupNames[3], ", loading field :", boundary1_fieldName, " in file ", boundary1_inputfile)
	boundary1Field=solverlab.Field(boundary1_inputfile, solverlab.CELLS, boundary1_fieldName, boundary1_time_iteration, boundary1_time_sub_iteration, boundary1_meshLevel)
	boundary2_type=solverlab.DirichletTransport
	boundary2_inputfile="../resources/meshSquare";
	boundary2_fieldName="Right temperature";
	boundary2_time_iteration=0# default value is 0
	boundary2_time_sub_iteration=0# default value is 0
	boundary2_meshLevel=0# default value is 0
	print("Boundary ", boundaryGroupNames[2], ", loading field :", boundary2_fieldName, " in file ", boundary2_inputfile)
	boundary2Field=solverlab.Field(boundary2_inputfile, solverlab.CELLS, boundary2_fieldName, boundary2_time_iteration, boundary2_time_sub_iteration, boundary2_meshLevel)
	boundary3_type=solverlab.NeumannTransport
	boundary3_inputfile="../resources/meshSquare";
	boundary3_fieldName="Top temperature";
	boundary3_time_iteration=0# default value is 0
	boundary3_time_sub_iteration=0# default value is 0
	boundary3_meshLevel=0# default value is 0
	print("Boundary ", boundaryGroupNames[4], ", loading field :", boundary3_fieldName, " in file ", boundary3_inputfile)
	boundary3Field=solverlab.Field(boundary3_inputfile, solverlab.CELLS, boundary3_fieldName, boundary3_time_iteration, boundary3_time_sub_iteration, boundary3_meshLevel)
	boundary4_type=solverlab.DirichletTransport
	boundary4_inputfile="../resources/meshSquare";
	boundary4_fieldName="Bottom temperature";
	boundary4_time_iteration=0# default value is 0
	boundary4_time_sub_iteration=0# default value is 0
	boundary4_meshLevel=0# default value is 0
	print("Boundary ", boundaryGroupNames[1], ", loading field :", boundary4_fieldName, " in file ", boundary4_inputfile)
	boundary4Field=solverlab.Field(boundary4_inputfile, solverlab.CELLS, boundary4_fieldName, boundary4_time_iteration, boundary4_time_sub_iteration, boundary4_meshLevel)

	# for each boundary we need to know if we want a Neumann or a Dirichlet boundary condition
	if boundary1_type==solverlab.NeumannTransport :
		myProblem.setNeumannBoundaryCondition("Left", boundary1Field)
	elif boundary1_type==solverlab.DirichletTransport :
		myProblem.setDirichletBoundaryCondition("Left", boundary1Field)
	if boundary2_type==solverlab.NeumannTransport :
		myProblem.setNeumannBoundaryCondition("Right", boundary2Field)
	elif boundary2_type==solverlab.DirichletTransport :
		myProblem.setDirichletBoundaryCondition("Right", boundary2Field)
	if boundary3_type==solverlab.NeumannTransport :
		myProblem.setNeumannBoundaryCondition("Top", boundary3Field)
	elif boundary3_type==solverlab.DirichletTransport :
		myProblem.setDirichletBoundaryCondition("Top", boundary3Field);
	if boundary4_type==solverlab.NeumannTransport :
		myProblem.setNeumannBoundaryCondition("Bottom", boundary4Field)
	elif boundary4_type==solverlab.DirichletTransport :
		myProblem.setDirichletBoundaryCondition("Bottom", boundary4Field);

    # set the numerical method
	myProblem.setTimeScheme( solverlab.Explicit)# Otherwise solverlab.Implicit
	max_nb_its_lin_solver = 50
	myProblem.setLinearSolver(solverlab.GMRES, solverlab.ILU, max_nb_its_lin_solver );

    # computation parameters
	MaxNbOfTimeStep = 3 ;# default value is 10
	freqSave = 1;# default value is 1
	cfl = 0.95;# default value is 1
	maxTime = 100000000;# default value is 10
	precision = 1e-6;# default value is 1e-6
	result_directory="."# default value = current directory

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setResultDirectory(result_directory)
	myProblem.setSaveFileFormat(solverlab.MED)#default value is solverlab.VTK

    # evolution
	myProblem.initialize();
	print("Running python test "+ fileName );

	ok = myProblem.run();
	if (ok):
		print( "Python simulation " + fileName + " is successful !" );
		pass
	else:
		print( "Python simulation " + fileName + "  failed ! " );
		pass

	print( "------------ End of simulation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	# name of result file
	fileName = "2DSpherical";# default value is ""
	TransportEquation_2DSpherical( fileName)
