#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2DCylinderDeflection():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	inputfile="/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/AnnulusSpiderWeb5x16.med"

	M=svl.Mesh(inputfile);
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the boundary conditions
	def initialPressure(x,y):
		return 10
	def initialBoundPressure(x,y):
		return 4
	def initialVelocity(x,y):
		return [-math.cos(x)*y,y*y]#[ x/np.sqrt((x*x)+ (y*y)),y/np.sqrt((x*x)+ (y*y))]
	def initialBoundVelocity(x,y):
		return [-3,10]#[x/np.sqrt((x*x)+ (y*y)),y/np.sqrt((x*x)+ (y*y))]
	
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
			Velocity0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) 
		elif (Fj.getNumberOfCells()==1):
			wallPressureMap[j] = initialBoundPressure(Ctemp1.x(),Ctemp1.y()) 
			if ( np.sqrt( Ctemp1.x()**2 + Ctemp1.y()**2 ) <= 2): #in fact 1.2 is enough since radius of small circle (on which we impose wallbound conditions)is 1.2
				#myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] = 0
			else :
				wallVelocityMap[j] = np.dot(initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma)	
		
	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

	

    # set the numerical metho50
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DCylinderDeflection";

	# computation parameers
	MaxNbOfTimeStep = 200000
	freqSave = 70
	maxTime = 200
	cfl =0.4
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
	WaveStaggered_2DCylinderDeflection()
