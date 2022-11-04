#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl

def IsothermalSinglePhase_2DLidDrivenCavity():
 
	spaceDim = 2;
    #Preprocessing: mesh data
	xinf=0;
	xsup=1;
	yinf=0;
	ysup=1;
	
	nx=50;
	ny=50;

    # set the limit field for each boundary

	fixedWallVelocityX=0;
	fixedWallVelocityY=0;

	movingWallVelocityX=1;
	movingWallVelocityY=0;
	
    # physical constants
	viscosite=[0.025];

	myProblem = svl.IsothermalSinglePhase(svl.Gas,svl.around1bar300K,spaceDim);#,False
	nVar = myProblem.getNumberOfVariables();

	#Initial field creation
	print("Building initial data " ); 
	
    # Prepare for the initial condition
    
	VV_Constant = [0] * nVar

	# constant vector
	VV_Constant[0] = 1e5;
	VV_Constant[1] = 0 ;
	VV_Constant[2] = 0;

    #Initial field creation
	print("Building mesh and initial data" );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"fixedWall","fixedWall",yinf,ysup,ny,"fixedWall","movingWall");

    # Set the boundary conditions
	myProblem.setWallBoundaryCondition("fixedWall", fixedWallVelocityX, fixedWallVelocityY);
	myProblem.setWallBoundaryCondition("movingWall", movingWallVelocityX, movingWallVelocityY);

    # set physical parameters
	myProblem.setViscosity(viscosite);

    # set the numerical method
	myProblem.setNumericalScheme(svl.staggered, svl.Implicit);
	myProblem.setLinearSolver(svl.GMRES,svl.LU);
   
    # name file save
	fileName = "2DLidDrivenCavity_Compressible_50x50_CFL10";

    # simulation parameters
	MaxNbOfTimeStep = 100001 ;
	freqSave = 100;
	cfl = 10;
	maxTime = 50000000000;
	precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(float('inf'),50);#1000000*
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
    IsothermalSinglePhase_2DLidDrivenCavity()
