#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np

def WaveStaggered_2DCylinderDeflection(n):
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	nr = 5*n
	ntheta = 16*n
	inputfile="/volatile/catC/esteban/SOLVERLAB/SOLVERLAB_SRC/CoreFlows/examples/resources/AnnulusSpiderWeb"+ str(nr)+"x" +str(ntheta)+".med"
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
	ExactVelocityInftyAtFaces = svl.Field("ExactVelocityInftyAtFaces", svl.FACES, M, 1)
	ExactVelocityInftyAtCells = svl.Field("ExactVelocityInftyAtCells", svl.CELLS, M, 3)
	
	
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		vec_normal_sigma = np.zeros(2)
		Ctemp1 = M.getCell(idCells[0]);
		#set orientation
		for l in range( Ctemp1.getNumberOfFaces()) :
				if (j == Ctemp1.getFacesId()[l]):
					for idim in range(spaceDim):
						vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
		myProblem.setOrientation(j,vec_normal_sigma)

		### Exact long time limit at CELLS ###
		r =  np.sqrt( Ctemp1.x()**2 + Ctemp1.y()**2 )
		theta = np.arctan(Ctemp1.y()/Ctemp1.x())
		for k in range(2):
			ExactVelocityInftyAtFaces[idCells[0],k] = ExactVelocity(r, theta, r1, r0)[k]

		if(Fj.getNumberOfCells()==2):
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y()) 
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp2.y())	
			Velocity0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) 
			### Exact long time limit at FACES ###
			r =  np.sqrt( Fj.x()**2 + Fj.y()**2 )
			theta = np.arctan(Fj.y()/Fj.x())
			ExactVelocityInftyAtFaces[j] = np.dot(ExactVelocity(r, theta, r1, r0),vec_normal_sigma ) 
			### Exact long time limit at CELLS ###
			r =  np.sqrt( Ctemp2.x()**2 + Ctemp2.y()**2 )
			theta = np.arctan(Ctemp2.y()/Ctemp2.x())
			for k in range(2):
				ExactVelocityInftyAtFaces[idCells[1], k] = ExactVelocity(r, theta, r1, r0)[k]

		elif (Fj.getNumberOfCells()==1):
			# if face is on interior (Wall boundary condition) r_int = 0.6
			if ( np.sqrt( Fj.x()**2 + Fj.y()**2 )  ) <= (r0 +r1)/2.0:  
				myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] = 0
			# if face is on exterior (Stegger-Warming condition) 
			else : 											
				wallVelocityMap[j] = np.dot(initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma)	
				wallPressureMap[j] = initialBoundPressure(Ctemp1.x(),Ctemp1.y()) 				
			ExactVelocityInftyAtFaces[j] = wallVelocityMap[j]

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical metho50
	myProblem.setTimeScheme(svl.Explicit);
	# name of result file
	fileName = "WaveStaggered_2DCylinderDeflection";

	# computation parameers
	MaxNbOfTimeStep = 100000000
	freqSave = 40000
	maxTime = 2000
	cfl =0.6
	precision = 1e-10;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.VTK)
	myProblem.saveVelocity();
	myProblem.savePressure(True);
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

	normL2 = myProblem.ErrorL2VelocityInfty(ExactVelocityInftyAtFaces, ExactVelocityInftyAtCells)
	sizeMesh = M.getNumberOfCells()
	myProblem.terminate();
	normL2faces = normL2[0]
	normL2cells = normL2[1]
	return [sizeMesh, normL2faces, normL2cells]

if __name__ == """__main__""":
	NormL2faces = []
	NormL2cells = []
	sizeMesh = []
	for i in range(5):
		N = 2**i
		result = WaveStaggered_2DCylinderDeflection(N)
		print("mesh size = ", result[0], "L2 error on velocity at infinity at FACES= ", result[1], "L2 error on velocity at infinity at CELLS= ", result[2])
		sizeMesh.append(result[0])
		NormL2faces.append(result[1])
		NormL2cells.append(result[2])
	plt.figure()
	plt.loglog(sizeMesh, NormL2cells,label = "L2 error on velocity infty at CELLS")
	plt.loglog(sizeMesh, NormL2faces,label = "L2 error on velocity infty at FACES")
	plt.loglog(sizeMesh, sizeMesh,label = "order 1")
	plt.legend()
	plt.title("error on velocity infty")
	plt.savefig("error on velocity infty")
