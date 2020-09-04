#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf

def FiveEqsTwoFluid_2DInclinedSedimentation():
	spaceDim = 2;

    # Prepare for the mesh
	xinf = 0 ;
	xsup=2.0;
	yinf=0.0;
	ysup=4.0;
	nx=20;
	ny=40; 

    # set the limit field for each boundary
	wallVelocityX=[0];
	wallVelocityY=[0];
	wallTemperature=573;

    # set the initial field  
	initialVoidFraction=0.5;
	initialVelocityX=[0]*2;
	initialVelocityY=[1]*2;
	initialTemperature=573;
	initialPressure=155e5;

	# physical constants
	gravite = [0] * spaceDim
    
	gravite[1]=-7;
	gravite[0]=7;

	myProblem = cf.FiveEqsTwoFluid(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = initialVoidFraction;
	VV_Constant[1] = initialPressure ;
	VV_Constant[2] = initialVelocityX[0];
	VV_Constant[3] = initialVelocityY[0];
	VV_Constant[4] = initialVelocityX[1];
	VV_Constant[5] = initialVelocityY[1];
	VV_Constant[6] = initialTemperature ;

    #Initial field creation
	print("Building initial data " );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
					  yinf,ysup,ny,"wall","wall", 
					  0.0,0.0,  0,  "", "")

    # the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setEntropicCorrection(True);
    
	# name of result file
	fileName = "2DInclinedSedimentation";

	# simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.25;
	maxTime = 5;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass

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
    FiveEqsTwoFluid_2DInclinedSedimentation()
