#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as sl
import math

def WaveStaggered_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=100;
	M=svl.Mesh(xinf,xsup,nx)

    # Prepare initial data
	kappa =2;
	rho = 5;
	myProblem = svl.WaveStaggered(spaceDim, kappa, rho );

    # Prepare for the initial condition
	initialVelocity=1;
	initialPressure=155e7;

    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant(M, initialVelocity, FACES);
	myProblem.setInitialFieldConstant(M, initialPressure, CELLS);

    # set the boundary conditions
	def boundPressure(x):
		return sin(x)

	def boundVelocity(x):
		return cos(x)

	wallPressureMap = [];
	wallVelocityMap = [];
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		isBoundary = Fj.isBorder;
		if (isBoundary == True):
			wallPressureMap.append([j, boundPressure(Fj.x()) ]);
			wallVelocityMap.append([j, boundVelocity(Fj.x()) ]);

	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
    
    # name of result file
	fileName = "WaveStaggered_1DRiemannProblem";

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
	myProblem.setSaveFileFormat(cf.CSV)

 
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
