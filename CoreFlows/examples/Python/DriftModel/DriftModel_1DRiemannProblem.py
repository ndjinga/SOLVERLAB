#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def DriftModel_1DRiemannProblem():

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

	inletConc_Left=0;
	inletVelocity_Left=1;
	inletTemperature_Left=565;
	outletPressure_Left=155e5;

	inletConc_Right=0;
	inletVelocity_Right=1;
	inletTemperature_Right=565;
	outletPressure_Right=155.1e5;

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

        # Prepare for the initial condition
	VV_Constant_Left =cm.Vector(nVar)
	VV_Constant_Right =cm.Vector(nVar)
	
	# constant vectors		
	VV_Constant_Left[0] = inletConc_Left ;
	VV_Constant_Left[1] = outletPressure_Left;
	VV_Constant_Left[2] = inletVelocity_Left;
	VV_Constant_Left[3] = inletTemperature_Left ;

	VV_Constant_Right[0] = inletConc_Right ;
	VV_Constant_Right[1] = outletPressure_Right;
	VV_Constant_Right[2] = inletVelocity_Right;
	VV_Constant_Right[3] = inletTemperature_Right ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,VV_Constant_Left,VV_Constant_Right,discontinuity)

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary")
	myProblem.setNeumannBoundaryCondition("RightBoundary");

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    
    # name of result file
	fileName = "1DRiemannProblemUpwind";

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
	myProblem.saveConservativeField(True);
 
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
    DriftModel_1DRiemannProblem()
