#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math

def WaveStaggered_Statio():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1.0;
	nx=4;
	M=svl.Mesh(xinf,xsup,nx)

    # Prepare initial data
	kappa = 2;
	rho = 5;
	myProblem = svl.WaveStaggered(spaceDim, kappa, rho );

    # Prepare for the initial condition
	initialVelocity=1;
	initialPressure=155e7;

    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant(M, [initialVelocity], svl.FACES);
	myProblem.setInitialFieldConstant(M, [initialPressure], svl.CELLS);

    # set the boundary conditions
	def boundPressure(x):
		return 155e7

	def boundVelocity(x):
		return 1

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
	fileName = "WaveStaggered_Statio";

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
	myProblem.setSaveFileFormat(svl.CSV)
	myProblem.setVerbose(True)

 
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
	assert(myProblem.isStationary()==True);

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    WaveStaggered_Statio()
