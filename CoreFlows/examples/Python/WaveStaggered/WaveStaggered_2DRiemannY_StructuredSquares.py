#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2DRiemannY_StructuredSquares():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0.0;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 1.0;
	discontinuity = (yinf + ysup)/2.0
	nx=100;
	ny=100
	M=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)#Regular square mesh
	print( "Built a regular 2D square mesh with ", nx,"x" ,ny, " cells")
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the boundary conditions
	initialVelocity_Left=3;
	initialPressure_Left=-1;

	initialVelocity_Right=1;
	initialPressure_Right=3;

	def initialPressure(Z):
		if Z < discontinuity:
			return initialPressure_Left
		else :
			return initialPressure_Right

	def initialVelocity(vec_normal,Z):
		vec_y = np.array([0,1])
		if (np.dot(vec_normal, vec_y)== 0):
			return 0
		else :
			if Z < discontinuity:
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
		idCells = Fj.getCellsId();
		vec_normal = np.zeros(2)
		if(Fj.getNumberOfCells()==2):
			Ctemp1 = M.getCell(idCells[0]);
			Ctemp2 = M.getCell(idCells[1]);
			PressureMap[idCells[0]] = initialPressure(Ctemp1.y()) ;
			PressureMap[idCells[1]] = initialPressure(Ctemp2.y());
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
				
			VelocityMap[j] = initialVelocity(vec_normal,Fj.y()) ;
		else:
			Ctemp1 = M.getCell(idCells[0]);
			wallPressureMap[j] = initialPressure(Ctemp1.y()) ;
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
			wallVelocityMap[j] = initialVelocity(vec_normal,Fj.y()) ;

	myProblem.setInitialFieldFunction(M, PressureMap, svl.CELLS, "pressure");
	myProblem.setInitialFieldFunction(M, VelocityMap, svl.FACES, "velocity");
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DRiemannY_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 6000 ;
	freqSave = 80;
	cfl = 0.4; 
	maxTime = 10;
	precision = 1e-6;

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
	WaveStaggered_2DRiemannY_StructuredSquares()
