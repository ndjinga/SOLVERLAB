#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
def WaveStaggered_2D_StructuredSquares():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0.0;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 1.0;
	nx=30;
	ny=30; 
	M=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)#Regular square mesh
	discontinuity=(xinf+xsup)/2 + 0.75/nx #TODO à définir de manière vectorielle

	
	print( "Built a regular 2D square mesh with ", nx,"x" ,ny, " cells")
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);
	#myProblem.setMesh(M);

	# Prepare for the initial condition
	initialVelocity_Left=3;
	initialPressure_Left=-1;

	initialVelocity_Right=1;
	initialPressure_Right=3;

	def initialPressure(x,y=0):
		if x < discontinuity:
			return initialPressure_Left
		else :
			return initialPressure_Right

	def initialVelocity(x,y=0):
		if x < discontinuity:
			return initialVelocity_Left
		else:
			return initialVelocity_Right

	 #Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	PressureMap = {};
	VelocityMap = {}; 
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		if (Fj.getNumberOfCells()==1):
			wallPressureMap[j] = initialPressure(Fj.x()) ;
			wallVelocityMap[j] = initialPressure(Fj.x()) ;
		else:
			PressureMap[j] = initialPressure(Fj.x()) ;
			VelocityMap[j] = initialVelocity(Fj.x()) ;

	myProblem.setInitialFieldFunction(M, PressureMap, svl.CELLS, "pressure");
	myProblem.setInitialFieldFunction(M, VelocityMap, svl.FACES, "velocity");
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2D_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 10 ;
	freqSave = 1;
	cfl = 0.4; 
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.CSV)
	myProblem.saveVelocity();
	myProblem.savePressure();
	myProblem.setVerbose(False);

	# Run the computation
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.solveStationaryProblem();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	WaveStaggered_2D_StructuredSquares()
