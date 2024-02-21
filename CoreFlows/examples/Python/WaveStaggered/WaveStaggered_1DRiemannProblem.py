#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math


def WaveStaggered_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1
	nx=50;
	discontinuity=(xinf+xsup)/2
	M=svl.Mesh(xinf,xsup,nx)


    # set the limit field for each boundary
	kappa = 2;
	rho = 5;

	initialVelocity_Left=1;
	initialPressure_Left=155e5;

	initialVelocity_Right=1;
	initialPressure_Right=150e5;

	myProblem = svl.WaveStaggered(spaceDim, rho, kappa);

        # Prepare for the initial condition
	Pressure_Left =svl.Vector(1);
	Pressure_Right =svl.Vector(1);
	Velocity_Left =svl.Vector(1);
	Velocity_Right =svl.Vector(1);
	
	# left and right constant vectors		
	Pressure_Left[0] = initialPressure_Left;
	Pressure_Right[0] = initialPressure_Right;
	Velocity_Left[0] = initialVelocity_Left;
	Velocity_Right[0] = initialVelocity_Right;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,Pressure_Left,Pressure_Right,discontinuity, svl.CELLS);
	myProblem.setInitialFieldStepFunction(M,Velocity_Left,Velocity_Right,discontinuity, svl.FACES);

    # set the boundary conditions
	#TODO : set boundary cond for Riemann problem ?
	# set the boundary conditions
	def boundPressure(x):
		if x < discontinuity:
			return initialPressure_Left
		else :
			return initialPressure_Right

	def boundVelocity(x):
		if x < discontinuity:
			return initialVelocity_Left
		else:
			return initialVelocity_Right

	wallPressureMap = {};
	wallVelocityMap = {}; 
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		isBoundary = Fj.isBorder;
		if (isBoundary == True):
			wallPressureMap[j] = boundPressure(Fj.x()) ;
			wallVelocityMap[j] = boundVelocity(Fj.x()) ;

	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);

    
    # name of result file
	fileName = "WaveStaggered_1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.2;
	maxTime = 500;
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
    WaveStaggered_1DRiemannProblem()
