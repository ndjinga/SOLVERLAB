#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def DriftModel_1DDepressurisation():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;

    # set the boundary field for each boundary
	initialConc=0;
	initialVelocityX=0;
	initialTemperature=565;
	initialPressure=155e5
    # set the boundary field for each boundary
	wallVelocityX=0;
	wallTemperature=565;
	outletPressure=1e5;

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = initialConc;
	VV_Constant[1] = initialPressure ;
	VV_Constant[2] = initialVelocityX;
	VV_Constant[3] = initialTemperature ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"wall","outlet");

    # set the boundary conditions
	myProblem.setWallBoundaryCondition("wall",wallTemperature,wallVelocityX)
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup]);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit); 
    	myProblem.setEntropicCorrection(True);

    # name of result file
	fileName = "1DDepressurisation";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 1;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision*1e7,20);
	myProblem.saveConservativeField(True);
	myProblem.saveAllFields(True);
 
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
    DriftModel_1DDepressurisation()
