#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2DRiemannX_StructuredSquares():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0.0;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 1.0;
	discontinuity = (xinf + xsup)/2.0
	nx=70;
	ny=70; 
	M=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)#Regular square mesh

	
	print( "Built a regular 2D square mesh with ", nx,"x" ,ny, " cells")
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the boundary conditions

	def initialPressure(Z):
		if Z < discontinuity:
			return -3
		else :
			return 0

	def initialVelocity(vec_normal,Z):
		vec_y = np.array([1,0])
		if (np.dot(vec_normal, vec_y)== 0):
			return 0
		else :
			if Z < discontinuity:
				return -4
			else:
				return -1
		

	#Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1);
	
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		vec_normal = np.zeros(2)
		if(Fj.getNumberOfCells()==2):
			Ctemp1 = M.getCell(idCells[0]);
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x()) ;
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x());
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
				
			Velocity0[j] = initialVelocity(vec_normal,Fj.x()) ;
		else:
			Ctemp1 = M.getCell(idCells[0]);
			wallPressureMap[j] = initialPressure(Ctemp1.x()) ;
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
			wallVelocityMap[j] = initialVelocity(vec_normal,Fj.x()) ;

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

	# Set periodicity on top/bottom faces
	myProblem.setVerticalPeriodicFaces()

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DRiemannX_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 2000 ;
	freqSave = 20;
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
	WaveStaggered_2DRiemannX_StructuredSquares()
