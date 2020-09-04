#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf

def FiveEqsTwoFluid_2DInclinedBoilingChannel():
	spaceDim = 2;

    # Prepare for the mesh
	xinf = 0 ;
	xsup=2;
	yinf=0.0;
	ysup=4;
	nx=20;
	ny=40; 

    # set the limit field for each boundary
	wallVelocityX=[0];
	wallVelocityY=[0];
	wallTemperature=573;
	inletVoidFraction=0;
	inletVelocityX=[0]*2;
	inletVelocityY=[1]*2;
	inletTemperature=573;
	outletPressure=155e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[1]=-7;
	gravite[0]=7;

	heatPower=5e7;

	myProblem = cf.FiveEqsTwoFluid(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = inletVoidFraction;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX[0];
	VV_Constant[3] = inletVelocityY[0];
	VV_Constant[4] = inletVelocityX[1];
	VV_Constant[5] = inletVelocityY[1];
	VV_Constant[6] = inletTemperature ;

    #Initial field creation
	print("Building mesh and initial data " );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
					  yinf,ysup,ny,"inlet","outlet", 
					  0.0,0.0,  0,  "", "")

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);
	myProblem.setInletBoundaryCondition("inlet", inletVoidFraction, inletTemperature, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
	myProblem.setEntropicCorrection(True);
	#myProblem.setWellBalancedCorrection(True);    
    
	# name of result file
	fileName = "2DInclinedBoilingChannel";

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
	#myProblem.saveConservativeField(True);
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
    FiveEqsTwoFluid_2DInclinedBoilingChannel()
