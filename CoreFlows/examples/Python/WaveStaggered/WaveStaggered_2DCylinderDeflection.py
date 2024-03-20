#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
import numpy as np
def WaveStaggered_2DCylinderDeflection():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	inputfile="/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/cylinder_geom_salome.med"

	M=svl.Mesh(inputfile);
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the boundary conditions

	def initialPressure(x,y):
		return 0
	def initialBoundPressure(x,y):
		return 

	def initialVelocity(vec_normal,x,y):
		return 0
	def initialBoundVelocity(x,y):
		return 4
		
	

		

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
			Pressure0[idCells[0]] = initialPressure(Ctemp1.y()) ;
			Pressure0[idCells[1]] = initialPressure(Ctemp2.y());
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
				
			Velocity0[j] = initialVelocity(vec_normal,Fj.y()) ;
		else:
			Ctemp1 = M.getCell(idCells[0]);
			wallPressureMap[j] = initialBoundPressure(Ctemp1.y()) ;
			for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
			wallVelocityMap[j] = initialBoundVelocity(vec_normal,Fj.y()) ;

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DCylinderDeflection";

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
	WaveStaggered_2DCylinderDeflection()
