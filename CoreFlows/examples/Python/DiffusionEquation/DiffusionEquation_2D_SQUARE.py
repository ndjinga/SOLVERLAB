#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab
from math import pi, sin, cos, exp

#===============================================================================================================================
# Name        : Simulation of a 2D heat equation 
# Description : Test solving the diffusion of the temperature T in a solid subject to a heat source
#               \rho cp dT/dt-\lambda\Delta T=\Phi 
#               Neumann or Dirichlet boundary conditions
#               Finite elements or finite volumes
#               The heat source function is an eigenvector of the Laplacean, so the exact solution is known using separation of variables
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
#================================================================================================================================


def DiffusionEquation_2D_SQUARE(FECalculation, fileName, DirichletBC):

	""" Description : Test solving the diffusion of the temperature T in a solid subject to a heat source. 
		Equation : Thermal diffusion equation  \rho cp dT/dt-\lambda\Delta T=\Phi 
		        Heat capacity cp, density \rho, and conductivity \lambda of the solid MUST be defined.
		        The solid may receive some extra heat power due to nuclear fissions or magnetic waves using function setHeatSource.
                The heat source function is an eigenvector of the Laplacean, so the exact solution is known using separation of variables.
	"""
	#Space dimension of the problem
	spaceDim=2
	
	# Prepare for the mesh
	print("Building mesh " )
	xmin=0
	xmax=1
	nx=20
	ymin=0
	ymax=1
	ny=20
	
	if(FECalculation):
		square_mesh = solverlab.Mesh(xmin,xmax,nx,ymin,ymax,ny,0)#right triangle mesh build from spliting the cells of a nx*ny cartesian mesh
	else:
		square_mesh = solverlab.Mesh(xmin,xmax,nx,ymin,ymax,ny)#nx*ny cartesian finite volume mesh
	
	if( square_mesh.getSpaceDimension()!=2 or square_mesh.getMeshDimension()!=2) :
		raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")
	if(FECalculation and not square_mesh.isTriangular()) :
		raise ValueError("Wrong cell types : mesh is not made of triangles")
	
	nbNodes = square_mesh.getNumberOfNodes()
	nbCells = square_mesh.getNumberOfCells()

	print( "Built a 2D square mesh with ", nbNodes, " nodes and ",nbCells, " cells")

	# Definition of field support parameter
	if( FECalculation):
		supportOfField=solverlab.NODES
	else:
		supportOfField=solverlab.CELLS	

	#Set the right hand side function
	heatPowerField = solverlab.Field("Heat power", supportOfField, square_mesh, 1)
	exactSolution = solverlab.Field("Exact temperature", supportOfField, square_mesh, 1)#eigenvector associated to -laplacian on a square
	if(DirichletBC):
		if( FECalculation):
			for i in range(square_mesh.getNumberOfNodes()):
				Ni= square_mesh.getNode(i)
				x = Ni.x()
				y = Ni.y()
				z = Ni.z()
		
				heatPowerField[i]=sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
				exactSolution[i]=sin(pi*x)*sin(pi*y)
		else:
			for i in range(square_mesh.getNumberOfCells()):
				Ci= square_mesh.getCell(i)
				x = Ci.x()
				y = Ci.y()
				z = Ci.z()

				heatPowerField[i]=sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
				exactSolution[i]=sin(pi*x)*sin(pi*y)
	else:
		if( FECalculation):
			for i in range(square_mesh.getNumberOfNodes()):
				Ni= square_mesh.getNode(i)
				x = Ni.x()
				y = Ni.y()
				z = Ni.z()
		
				heatPowerField[i]=cos(pi*x)*cos(pi*y)#mettre la fonction definie au second membre de l'edp
				exactSolution[i]=cos(pi*x)*cos(pi*y)
		else:
			for i in range(square_mesh.getNumberOfCells()):
				Ci= square_mesh.getCell(i)
				x = Ci.x()
				y = Ci.y()
				z = Ci.z()
		
				heatPowerField[i]=cos(pi*x)*cos(pi*y)#mettre la fonction definie au second membre de l'edp
				exactSolution[i]=cos(pi*x)*cos(pi*y)
				
	heatPowerField.writeVTK("HeatPowerField")
	
	initialTemperature=0
	
    # Mandatory physical values
	solid_specific_heat=300# specific heat capacity, default value 300
	solid_density=10000# density, default value 10000
	solid_conductivity=5# conductivity, default value 5
	diffusivity=solid_conductivity/(solid_density*solid_specific_heat)

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,solid_density,solid_specific_heat,solid_conductivity);
	myProblem.setInitialFieldConstant(square_mesh,[initialTemperature],supportOfField);

    # set the Neumann boundary condition
	if(DirichletBC):
		myProblem.setDirichletBoundaryCondition("Boundary")
	else:
		myProblem.setNeumannBoundaryCondition("Boundary")
		
	myProblem.setHeatPowerField(heatPowerField)
	
    # set the numerical method
	myProblem.setTimeScheme( solverlab.Implicit)
	max_nb_its_lin_solver = 50
	myProblem.setLinearSolver(solverlab.GMRES, solverlab.ILU, max_nb_its_lin_solver );

    # computation parameters
	MaxNbOfTimeStep = 3 ;# default value is 10
	freqSave = 1;# default value is 1
	cfl = 100;# default value is 1
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
	#myProblem.setSaveFileFormat(solverlab.MED)#default value is solverlab.VTK

    # evolution
	myProblem.initialize();
	print("Running python test "+ fileName );

	ok = myProblem.run();
	if (ok):
		print( "Python simulation " + fileName + " is successful !" );
		pass
		
		final_time = myProblem.getTime()
		Lambda = 2*pi*pi#eigenvalue associated to the heat power field
		time_factor = 1/(Lambda*diffusivity)*(1-exp(-Lambda*diffusivity*final_time))#Exact solution is time factor times heat source function
		if( FECalculation):
			for i in range(square_mesh.getNumberOfNodes()):
				exactSolution[i]*=time_factor
		else:
			for i in range(square_mesh.getNumberOfCells()):
				exactSolution[i]*=time_factor
		exactSolution.writeVTK("exactSolution")
	else:
		print( "Python simulation " + fileName + "  failed ! " );
		pass

	print( "------------ End of simulation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    if len(sys.argv) >2 :
        FECalculation = bool(int(sys.argv[1]))
        DirichletBC   = bool(int(sys.argv[2]))
        # name of result file
        if( FECalculation):
            if(DirichletBC):
                fileName = "DiffusionEquation_2D_SQUARE_FE_Dirichlet";
            else:
                fileName = "DiffusionEquation_2D_SQUARE_FE_Neumann";
        else:
            if(DirichletBC):
                fileName = "DiffusionEquation_2D_SQUARE_FV_Dirichlet"
            else:
                fileName = "DiffusionEquation_2D_SQUARE_FV_Neumann";
        DiffusionEquation_2D_SQUARE(FECalculation, fileName,DirichletBC)
    else :
        raise ValueError("DiffusionEquation_2D_SQUARE_Neumann : missing arguments")
