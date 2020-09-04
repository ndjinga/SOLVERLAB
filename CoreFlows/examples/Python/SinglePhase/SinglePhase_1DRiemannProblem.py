#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;
	discontinuity=(xinf+xsup)/2
	M=cm.Mesh(xinf,xsup,nx)
	eps=1e-6
	M.setGroupAtPlan(xsup,0,eps,"RightBoundary")
	M.setGroupAtPlan(xinf,0,eps,"LeftBoundary")

    # set the limit field for each boundary

	initialVelocity_Left=1;
	initialTemperature_Left=565;
	initialPressure_Left=155e5;

	initialVelocity_Right=1;
	initialTemperature_Right=565;
	initialPressure_Right=150e5;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

        # Prepare for the initial condition
	VV_Left =cm.Vector(nVar)
	VV_Right =cm.Vector(nVar)
	
	# left and right constant vectors		
	VV_Left[0] = initialPressure_Left;
	VV_Left[1] = initialVelocity_Left;
	VV_Left[2] = initialTemperature_Left ;

	VV_Right[0] = initialPressure_Right;
	VV_Right[1] = initialVelocity_Right;
	VV_Right[2] = initialTemperature_Right ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary");
	myProblem.setNeumannBoundaryCondition("RightBoundary");

    # set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);
    
    # name of result file
	fileName = "1DRiemannProblem";

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
	myProblem.usePrimitiveVarsInNewton(True)
 
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
    SinglePhase_1DRiemannProblem()
