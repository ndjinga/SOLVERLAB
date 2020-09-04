#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def IsothermalTwoFluid_1DSedimentation():
 
	spaceDim = 1;
    # Prepare for the mesh
	xinf = 0. ;
	xsup=1.0;
	nx=50; 

    # set the limit field for each boundary
	wallVelocityX=[0] * 2;
	wallPressure=1.e5;

	initialVoidFraction=0.5;

    # physical constants
	gravite=[-10];

	myProblem = cf.IsothermalTwoFluid(cf.around1bar300K,spaceDim);
	nVar = myProblem.getNumberOfVariables();
	print(nVar)

    # Prepare for the initial condition
	VV_Constant =[0] *nVar

	# constant vector
	VV_Constant[0] = initialVoidFraction;
	VV_Constant[1] = wallPressure;
	VV_Constant[2] = wallVelocityX[0];
	VV_Constant[3] = wallVelocityX[1];

    #Initial field creation
	print("Building mesh and initial data " );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"wall","wall")

    # the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallVelocityX);

    # set physical parameters
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
	myProblem.setEntropicCorrection(True);
    
    # name file save
	fileName = "1DSedimentation";

    # simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 5.;
	maxTime = 5.;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

    # evolution
	myProblem.initialize();
	print("Running python"+ fileName );

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
    IsothermalTwoFluid_1DSedimentation()
