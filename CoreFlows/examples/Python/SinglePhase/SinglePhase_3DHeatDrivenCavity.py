#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf

def SinglePhase_3DHeatDrivenCavity():
 
	spaceDim = 3;
    #Preprocessing: mesh data
	xinf=0;
	xsup=1;
	yinf=0;
	ysup=1;
	zinf=0;
	zsup=1;
	
	nx=10;
	ny=10;
	nz=10;

    # set the limit field for each boundary

	coldWallVelocityX=0;
	coldWallVelocityY=0;
	coldWallVelocityZ=0;
	coldWallTemperature=563;

	hotWallVelocityX=0;
	hotWallVelocityY=0;
	hotWallVelocityZ=0;
	hotWallTemperature=613;
	
    # physical constants

	gravite = [0] * spaceDim
	gravite[2]=-10;
	gravite[1]=0;
	gravite[0]=0;
	viscosite=[8.85e-5];
	conductivite=[1000];#Wall heat transfert due to nucleate boiling.


	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar = myProblem.getNumberOfVariables();

	#Initial field creation
	print("Building initial data " ); 
	
    # Prepare for the initial condition
    
	VV_Constant = [0] * nVar

	# constant vector
	VV_Constant[0] = 155e5;
	VV_Constant[1] = 0 ;
	VV_Constant[2] = 0;
	VV_Constant[3] = 0;
	VV_Constant[4] = 573;

    #Initial field creation
	print("Setting mesh and initial data" );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"hotWall","hotWall",yinf,ysup,ny,"hotWall","hotWall",zinf,zsup,nz, "hotWall", "coldWall");

    # Set the boundary conditions
	myProblem.setWallBoundaryCondition("coldWall", coldWallTemperature, coldWallVelocityX, coldWallVelocityY, coldWallVelocityZ);
	myProblem.setWallBoundaryCondition("hotWall", hotWallTemperature, hotWallVelocityX, hotWallVelocityY, hotWallVelocityZ);

    # set physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
	myProblem.setEntropicCorrection(False);
	myProblem.setWellBalancedCorrection(False);   
    
    # name file save
	fileName = "3DHeatDrivenCavity";

    # simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 10;
	maxTime = 50;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,50);
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

	print( "------------ End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_3DHeatDrivenCavity()
