#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def DriftModel_2DInclinedBoilingChannel():
	spaceDim = 2;

    # Prepare for the mesh
	xinf = 0 ;
	xsup=1.0;
	yinf=0.0;
	ysup=4.0;
	nx=10;
	ny=40; 

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=563;
	inletConcentration=0;
	inletVelocityX=0;
	inletVelocityY=1;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[1]=-8.5;
	gravite[0]=5;

	heatPower=0e8;

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = inletConcentration;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletVelocityY;
	VV_Constant[4] = inletTemperature ;

    #Initial field creation
	print("Building mesh and initial data " );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
					  yinf,ysup,ny,"inlet","outlet", 
					  0.0,0.0,  0,  "", "")

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConcentration, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);
	myProblem.setWellBalancedCorrection(True);    
	myProblem.setNonLinearFormulation(cf.VFFC) 
    
	# name of result file
	fileName = "2DInclinedChannelGravity";

	# simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
    DriftModel_2DInclinedBoilingChannel()
