#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2D_StructuredSquares():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0.0;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 1.0;
	nx=90;
	ny=90; 
	M=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)#Regular square mesh

	
	print( "Built a regular 2D square mesh with ", nx,"x" ,ny, " cells")
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);
	#myProblem.setMesh(M);

	# Prepare for the initial condition

	def initialPressure(x,y):
		return math.sin(2*math.pi*x)

	def initialVelocity(vec_y,vec_normal,x=0,y=0):
		if (np.dot(vec_normal,vec_y) ==0):
			return 0
		else: 
			return 1

	 #Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	PressureMap = {};
	VelocityMap = {}; 
	vec_normal = np.zeros(2)
	vec_y = np.array([0,1])
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		if(Fj.getNumberOfCells()==2):
			Ctemp1 = M.getCell(idCells[0]);
			Ctemp2 = M.getCell(idCells[1]);
			PressureMap[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.x())/2.0 ;
			PressureMap[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp1.y())/2.0 ;
			for l in range(Ctemp1.getNumberOfFaces()):
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
			VelocityMap[j] = initialVelocity(vec_y, vec_normal, Fj.x(),Fj.y()) ;
		else:
			Ctemp1 = M.getCell(idCells[0]);
			wallPressureMap[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y())/2.0 ;
			wallVelocityMap[j] = initialVelocity(vec_y, vec_normal, Fj.x(),Fj.y()) ;

	myProblem.setInitialFieldFunction(M, PressureMap, svl.CELLS, "pressure");
	myProblem.setInitialFieldFunction(M, VelocityMap, svl.FACES, "velocity");
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2D_StructuredSquares";

	# computation parameters
	MaxNbOfTimeStep = 2000 ;
	freqSave = 80;
	cfl = 0.4; 
	maxTime = 500;
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
	WaveStaggered_2D_StructuredSquares()
