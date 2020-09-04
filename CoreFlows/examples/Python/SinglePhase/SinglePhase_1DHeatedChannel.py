#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def SinglePhase_1DHeatedChannel():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;

    # set the limit field for each boundary

	inletVelocityX=5;
	inletTemperature=565;
	outletPressure=155e5;

    # physical parameters
	heatPower=1e7;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletTemperature ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"inlet","outlet");

    # set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet",inletTemperature,inletVelocityX)
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup]);

    # set physical parameters
	myProblem.setHeatSource(heatPower);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setWellBalancedCorrection(True);  
    
    # name of result file
	fileName = "1DHeatedChannelUpwindWB";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 500;
	precision = 1e-7;

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
    SinglePhase_1DHeatedChannel()
