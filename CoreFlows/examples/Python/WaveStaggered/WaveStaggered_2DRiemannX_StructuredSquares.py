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
	nx=30;
	ny=3; 
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
			return 2
		else :
			return 1

	def initialVelocity(vec_normal_sigma,Z):
		vec_x = np.array([1,0])
		if (np.dot(vec_normal_sigma, vec_x)== 0):
			return 0
		else :
			if Z < discontinuity:
				return 0
			else: 
				return 1
		

	#Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1);

	# Set periodicity on top/bottom faces
	myProblem.setVerticalPeriodicFaces()

	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		vec_normal_sigma = np.zeros(2)
		Ctemp1 = M.getCell(idCells[0]);
		for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);

		if(Fj.getNumberOfCells()==2):
			myProblem.setInteriorIndex(j);
			myProblem.setOrientation(j,vec_normal_sigma)
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x()) ;
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x());
			Velocity0[j] = initialVelocity(vec_normal_sigma,Fj.x()) ;
		else:
			for idim in range(spaceDim):
				if vec_normal_sigma[idim] < 0:	
					vec_normal_sigma[idim] = -vec_normal_sigma[idim]
			myProblem.setOrientation(j,vec_normal_sigma)
			if (Fj.x() <(xsup - xinf)/(3*nx) ) : # TODO on left side ?
				myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] = 0
			else :
				myProblem.setSteggerBoundIndex(j) #TODO trouver le moyen d'éviter la face de bord qui n'est pas calculée en périodique
				wallPressureMap[j] = initialPressure(Ctemp1.x()) ;
				wallVelocityMap[j] = initialVelocity(vec_normal_sigma,Fj.x()) ;

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);


    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DRiemannX_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 1500 ;
	freqSave = 30;
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