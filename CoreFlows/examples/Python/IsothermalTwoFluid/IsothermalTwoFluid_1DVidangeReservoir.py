#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def IsothermalTwoFluid_1DVidangeReservoir():
 
	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0. ;
	xsup=1.0;
	nx=50; 
	M=cm.Mesh(xinf,xsup,nx)
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xinf,0,eps,"wall")
	M.setGroupAtPlan(xsup,0,eps,"inlet")

    # set the limit field for each boundary
	wallVelocityX=[0] * 2

    # set the limit field for each boundary
	initialAlphaTop=1;
	initialAlphaBottom=0.;
	initialVelocityX=[0]*2;
	initialPressure=1e5;


    # physical constants
	gravite=[-10];

	myProblem = cf.IsothermalTwoFluid(cf.around1bar300K,spaceDim);
	nVar = myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_top =cm.Vector(nVar)
	VV_bottom =cm.Vector(nVar)

	# top and bottom vectors
	VV_top[0] = initialAlphaTop ;
	VV_top[1] = initialPressure ;
	VV_top[2] = initialVelocityX[0];
	VV_top[3] = initialVelocityX[1];

	VV_bottom[0] = initialAlphaBottom ;
	VV_bottom[1] = initialPressure ;
	VV_bottom[2] = initialVelocityX[0];
	VV_bottom[3] = initialVelocityX[1];

    #Initial field creation
	print("Building initial data " );
	myProblem.setInitialFieldStepFunction( M, VV_bottom, VV_top, .8, 0);

    # the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallVelocityX);
	myProblem.setNeumannBoundaryCondition("inlet");

    # set physical parameters
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	#myProblem.setEntropicCorrection(True);
	myProblem.setNonLinearFormulation(cf.VFFC) 
    
    # name file save
	fileName = "1DVidangeReservoir";

    # simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = .1;
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
    IsothermalTwoFluid_1DVidangeReservoir()
