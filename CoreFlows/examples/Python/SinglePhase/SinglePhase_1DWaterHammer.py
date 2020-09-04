#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_1DWaterHammer():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;
	M=cm.Mesh(xinf,xsup,nx)

    # set the initial field

	initialPressure=155e5;
	initialVelocityX=-5;
	initialTemperature=573;

   # set the limit field for each boundary

	wallVelocityX=0;
	wallTemperature=573;

 	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = initialPressure ;
	VV_Constant[1] = initialVelocityX;
	VV_Constant[2] = initialTemperature ;


    #Initial field creation
	print("Building initial data" ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"wall","neumann");

    # set the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX);
	myProblem.setNeumannBoundaryCondition("neumann");


    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    
    # name of result file
	fileName = "1DSinglePhase_1DWaterHammer";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 500;
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
    SinglePhase_1DWaterHammer()
