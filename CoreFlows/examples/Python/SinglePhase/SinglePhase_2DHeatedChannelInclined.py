#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def SinglePhase_2DHeatedChannelInclined():
	spaceDim = 2;

    # Prepare the mesh data
	xinf = 0 ;
	xsup=3.0;
	yinf=0.0;
	ysup=5.0;
	nx=50;
	ny=50; 

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=573;
	inletVelocityX=0;
	inletVelocityY=0.5;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[1]=-7;
	gravite[0]=7;

	heatPower=1e8;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletVelocityY;
	VV_Constant[3] = inletTemperature ;

    #Initial field creation
	print("Building mesh and initial data" );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
					  yinf,ysup,ny,"inlet","outlet", 
					  0.0,0.0,  0,  "", "")

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);
    # set physical parameters

	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);
	myProblem.setWellBalancedCorrection(False);
    
	# name file save
	fileName = "2DInclinedHeatedChannel";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
    SinglePhase_2DHeatedChannelInclined()
