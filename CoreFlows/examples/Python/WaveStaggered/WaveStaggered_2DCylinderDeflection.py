#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
import sys
def WaveStaggered_2DCylinderDeflection():
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	inputfile="/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/AnnulusSpiderWeb5x16.med"
	r0 = 0.8
	r1 = 6

	M=svl.Mesh(inputfile);
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim,rho, kappa);

	# Prepare for the initial condition
	# set the boundary conditions
	def ExactVelocity(r, theta, r1, r0):
		return [r1*r1/(r1*r1 -r0*r0)*(1 - r0*r0/(r*r) * math.cos(2*theta)),  r1*r1/(r1*r1 -r0*r0)*(- r0*r0/(r*r) * math.sin(2*theta))]
	def initialPressure(x,y):
		return 0
	def initialBoundPressure(x,y):
		return 0
	def initialVelocity(x,y):
		return [0,0] 
	def initialBoundVelocity(x,y):
		return [1,0]
	
	#Initial field creation
	print("Building initial data " ); 
	wallPressureMap = {};
	wallVelocityMap = {}; 
	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1);
	ExactVelocityInftyAtCells = svl.Field("ExactVelocityInftyAtCells", svl.CELLS, M, 3); 
	ExactVelocityInftyAtFaces = svl.Field("ExactVelocityInftyAtFaces", svl.FACES, M, 1)
	ExactVelocityInftyInterpolate = svl.Field("ExactVelocityInftyAtInterpolate", svl.CELLS, M, 3);
	for l in range(M.getNumberOfCells()):
		for k in range(spaceDim):
			ExactVelocityInftyInterpolate[l, k] =0;
	
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		vec_normal_sigma = np.zeros(2)
		Ctemp1 = M.getCell(idCells[0]);
		#set orientation
		found = False
		for l in range( Ctemp1.getNumberOfFaces()) :
			if (j == Ctemp1.getFacesId()[l]):
				found = True
				for idim in range(spaceDim):
					vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
		if (not found):
			sys.exit()

		myProblem.setOrientation(j,vec_normal_sigma)

		if(Fj.getNumberOfCells()==2):
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y()) 
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp2.y())	
			Velocity0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma )

			r =  np.sqrt( Fj.x()**2 + Fj.y()**2 )
			theta = np.arctan(Fj.y()/Fj.x())
			ExactVelocityInftyAtFaces[j] = np.dot(ExactVelocity(r, theta, r1, r0),vec_normal_sigma ) 
			for k in range(spaceDim):
					ExactVelocityInftyInterpolate[idCells[0], k] += ExactVelocityInftyAtFaces[j] * vec_normal_sigma[k]/Ctemp1.getNumberOfFaces()
					ExactVelocityInftyInterpolate[idCells[1], k] += ExactVelocityInftyAtFaces[j] * vec_normal_sigma[k]/Ctemp2.getNumberOfFaces() 
		elif (Fj.getNumberOfCells()==1):
			# if face is on interior (wallbound condition) r_int = 1.2 ou 0.8 selon le maillage
			if ( np.sqrt( Fj.x()**2 + Fj.y()**2 )  ) <= (r0 +r1)/2.0:  
				myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] =  0 # np.dot(initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma)	
			# if face is on exterior (stegger condition) 
			else : 											
				wallVelocityMap[j] = np.dot(initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma)	
				wallPressureMap[j] = initialBoundPressure(Ctemp1.x(),Ctemp1.y()) 	
			# building exact solution at faces and its interpolation at cells			
			ExactVelocityInftyAtFaces[j] = wallVelocityMap[j]
			for k in range(spaceDim):
				ExactVelocityInftyInterpolate[idCells[0], k] += ExactVelocityInftyAtFaces[j] * vec_normal_sigma[k]/Ctemp1.getNumberOfFaces()

	
	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical metho50
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DCylinderDeflection";

	# computation parameers
	MaxNbOfTimeStep = 1
	freqSave = 1
	maxTime = 10000
	cfl =0.6	 #Computed CFL = d/2 = 0.12 in quad 
	precision = 1e-8;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.VTK)
	myProblem.saveVelocity(True);
	myProblem.savePressure(True);
	myProblem.setVerbose(False);

	for l in range(M.getNumberOfCells()):
		Ctemp1 = M.getCell(l)
		rayon1 =  np.sqrt( Ctemp1.x()**2 + Ctemp1.y()**2 )
		theta1 = np.arctan(Ctemp1.y()/Ctemp1.x())
		for k in range(spaceDim):
			exa = ExactVelocity(rayon1, theta1, r1, r0)
			ExactVelocityInftyAtCells[l,k] = exa[k]
	myProblem.setExactVelocityFieldAtCells(ExactVelocityInftyAtCells)

	
	# Run the computation
	myProblem.initialize();

	ok = myProblem.run();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.ErrorRelativeVelocityInfty(ExactVelocityInftyAtFaces);
	normL2 = myProblem.ErrorL2VelocityInfty(ExactVelocityInftyAtFaces)
	sizeMesh = M.getNumberOfCells()
	print("nbmailles =", sizeMesh, "norme L2 =", normL2)
	
	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	WaveStaggered_2DCylinderDeflection()
