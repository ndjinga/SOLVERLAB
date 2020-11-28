#!/usr/bin/env python3
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_2DAirfoil():
	spaceDim = 2;

    # Prepare the mesh data
	print("Loading the mesh" );
	filename="./resources/naca10.med"
	my_mesh=cm.Mesh(filename)
       
	# Boundary data
	inletVelocityX   = 175;
	inletVelocityY   = 0.;
	inletTemperature = 300
	outletPressure   = 1e5
	wallVelocityX    = 0;
	wallVelocityY    = 0;
	wallTemperature  = 300;

    # physical constants
	gravite  = [0] * spaceDim
	gravite[0]= 0;
	gravite[1]= -10;
	viscosity = 15e-6
	viscosite = [viscosity];

	myProblem = cf.SinglePhase(cf.Gas,cf.around1bar300K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletVelocityY;
	VV_Constant[3] = inletTemperature ;

    #Initial field creation
	print("Setting initial data" );
	myProblem.setInitialFieldConstant(my_mesh,VV_Constant)

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("Inlet", outletPressure);
	myProblem.setInletBoundaryCondition( "Outlet", inletTemperature, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("Top",     wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("Bottom",  wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("Airfoil_boundary", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
    
	# name file save
	fileName = "2DAirfoil";

	# parameters calculation
	MaxNbOfTimeStep = 10000 ;
	freqSave = 100;
	cfl = 50;
	maxTime = 5000;
	precision = 1e-4;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setNewtonSolver(precision,50);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass

	# evolution
	myProblem.initialize();

	ok = myProblem.run();
	if (ok):
		print( "Simulation " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation " + fileName + "  failed ! " );
		pass

	print( "------------ End of Simulation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_2DAirfoil()
