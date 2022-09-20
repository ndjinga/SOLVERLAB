#!/usr/bin/env python3
# -*-coding:utf-8 -*

import CoreFlows as cf

def IsothermalSinglePhase_1DChannelGravity():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1;
	nx=10;


    # physical parameters
	gravite=[-10];

	myProblem = cf.IsothermalSinglePhase(cf.Liquid,cf.around155bars600K,spaceDim,False);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = 1e5 ;
	VV_Constant[1] = 0;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Top","Bottom");

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("Top")
	myProblem.setNeumannBoundaryCondition("Bottom");

    # set physical parameters
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);
    
    # name of result file
	fileName = "1DChannelGravityStaggered_Incompressible";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 1;
	maxTime = 500;
	precision = 1e-7;

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
    IsothermalSinglePhase_1DChannelGravity()
