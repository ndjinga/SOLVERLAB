#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def SinglePhase_2DPoiseuilleFlow():
	spaceDim = 2;

    # Prepare the mesh data
	xinf = 0 ;
	xsup=1.0;
	yinf=0.0;
	ysup=4;
	nx=10;
	ny=40; 

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=573;

    # physical constants
	viscosite=[0.025];

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	initialVelocityX=0;
	initialVelocityY=1;
	initialTemperature=573;
	initialPressure=155e5;
	VV_Constant[0] = initialPressure ;
	VV_Constant[1] = initialVelocityX;
	VV_Constant[2] = initialVelocityY;
	VV_Constant[3] = initialTemperature ;

    #Initial field creation
	print("Building mesh and initial data" );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
										  yinf,ysup,ny,"neumann","neumann", 
										  0.0,0.0,  0,  "", "")

    # the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setNeumannBoundaryCondition("neumann");

    # set physical parameters
	myProblem.setViscosity(viscosite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    
	# name file save
	fileName = "2DPoiseuilleFlowExplicit_Phase1";

	# parameters calculation
	MaxNbOfTimeStep = 10
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5000;
	precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setNewtonSolver(1e-3,20);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass

	# evolution
	myProblem.initialize();

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of Phase 1 !!! -----------" );

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);

	# name file save
	fileName = "2DPoiseuilleFlowImplicit_Phase2";

	cfl = 0.5;
	MaxNbOfTimeStep = 20 ;
	myProblem.setCFL(cfl);
	myProblem.setNewtonSolver(1e-6,20);
	myProblem.setFileName(fileName);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	
	# evolution
	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of phase 2 !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_2DPoiseuilleFlow()
