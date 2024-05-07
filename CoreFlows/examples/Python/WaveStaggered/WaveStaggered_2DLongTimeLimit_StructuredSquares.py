#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2DLongTimeLimit_StructuredSquares():
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

	
	print( "Built a regular 2D square mesh with ", nx,"x" ,ny, " cells")
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the initial interior conditions
	def initialPressure(x,y):
		x1 = x -0.5
		y1 = y - 0.5
		norm = math.sqrt(x1*x1 + y1*y1)
		if norm <0.2:
			return 1/2
		else :
			return 0	
	def initialVelocity(x,y):
		vec = np.array([0, 0])
		return vec
	# set the boundary conditions
	def initialBoundPressure(x,y):
		return 0		
	def initialBoundVelocity(x,y):
		vec = np.array([0, 0])
		return vec
	#Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1);
	
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		vec_normal_sigma = np.zeros(2)
		Ctemp1 = M.getCell(idCells[0]);
		
		for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
		myProblem.setOrientation(j,vec_normal_sigma)
		if(Fj.getNumberOfCells()==2):
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y()) 
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp2.y())	
			Velocity0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()), vec_normal_sigma)  ;
		else:
			wallPressureMap[j] = initialBoundPressure(Ctemp1.x(),Ctemp1.y()) ;
			wallVelocityMap[j] = np.dot(initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma) ;


	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap)

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DLongTimeLimit_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 50000 ;
	freqSave = 80;
	cfl = 0.4; 
	maxTime = 10;
	precision = 1e-3;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.VTK)
	myProblem.saveVelocity();
	myProblem.savePressure();
	myProblem.setVerbose(False);

	# Run the computation
	myProblem.initialize();

	ok = myProblem.run();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	WaveStaggered_2DLongTimeLimit_StructuredSquares()
