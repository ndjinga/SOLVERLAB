#!/usr/bin/env python
# -*-coding:utf-8 -*


import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np

RHO_b = 2
GAMMA = 2
C_b = np.sqrt(GAMMA * pow(RHO_b, GAMMA-1) ) 
M_REF = 1e-4
KAPPA = 1
def EulerBarotropicStaggered_2DCylinderDeflection(n):
	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	nr = 5*n
	ntheta = 16*n
	inputfile="/volatile/catB/esteban/Solverlab/INSTALL/share/examples/resources/AnnulusSpiderWeb"+ str(nr)+"x" +str(ntheta)+".med"
	r0 = 0.8
	r1 = 6

	M=svl.Mesh(inputfile);
	myProblem =  svl.EulerBarotropicStaggered(svl.GasStaggered, svl.around1bar300K, KAPPA, GAMMA, spaceDim )

	# Prepare for the initial condition
	# set the boundary conditions
	def ExactStationnaryVelocity(r, theta, r1, r0):
		return np.array([M_REF * C_b * r1*r1/(r1*r1 -r0*r0)*(1 - r0*r0/(r*r) * math.cos(2*theta)),  r1*r1/(r1*r1 -r0*r0)*(- r0*r0/(r*r) * math.sin(2*theta))])
	def initialDensity(x,y):
		return RHO_b
	def initialVelocity(x,y):
		return np.array([0, 0])
	def initialBoundVelocity(x,y):
		return np.array([M_REF * C_b,0])
		
	#Initial field creation
	print("Building initial data " ); 
	wallDensityMap = {};
	wallVelocityMap = {}; 
	wallVelocityVector = np.zeros(spaceDim)
	Density0 = svl.Field("Density", svl.CELLS, M, 1);
	Momentum0 = svl.Field("velocity", svl.FACES, M, 1);
	ExactVelocity = np.zeros(M.getNumberOfFaces())
	
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
		
		#Todo : orientaiton
		if (  Fj.x() >1e-10 and abs( np.arctan(Fj.y()/Fj.x()) ) <1e-10 ):
			vec_normal_sigma[0] *= -1;
			vec_normal_sigma[1] *= -1
		myProblem.setOrientation(j,vec_normal_sigma)

		ExactVelocity[j] = np.dot(ExactStationnaryVelocity( np.sqrt( Fj.x()**2 + Fj.y()**2 ), np.arctan(Fj.y()/Fj.x()), r1, r0),vec_normal_sigma ) 
		Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y()) 
		
		if(Fj.getNumberOfCells()==2):
			myProblem.setInteriorIndex(j);
			Density0[idCells[1]] = initialDensity(M.getCell(idCells[1]).x(),M.getCell(idCells[1]).y())	
			Momentum0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * ( Density0[idCells[0]] + Density0[idCells[1]]  )/2;
		elif (Fj.getNumberOfCells()==1):
			Momentum0[j] = np.dot(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * initialDensity(Fj.x(),Fj.y());
			# if face is on interior (Wall boundary condition) r_int = 0.8
			if ( np.sqrt( Fj.x()**2 + Fj.y()**2 )  ) <= (r0 +r1)/2.0:  
				myProblem.setWallBoundIndex(j);
				for idim in range(spaceDim):
					wallVelocityVector[idim] = 0;
			# if face is on exterior (Stegger-Warming condition) 
			else : 											
				myProblem.setSteggerBoundIndex(j);								
				wallDensityMap[j] = initialDensity(Fj.x(),Fj.y());
				for idim in range(spaceDim):
					wallVelocityVector[idim] = initialBoundVelocity(Fj.x(), Fj.y())[idim];
			wallVelocityMap[j] = np.dot( wallVelocityVector, vec_normal_sigma );
			myProblem.setboundaryVelocityVector(j, wallVelocityVector);

	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

     #name of result file
	fileName = "EulerBarotropicStaggered_2DCylinderDeflection";

    #parameters calculation
	MaxNbOfTimeStep = 1000000	;
	cfl = 100;
	precision = 1e-8;
	freqSave = 100;
	maxTime = 50;


	myProblem.setTimeScheme(svl.Implicit);
	myProblem.setLinearSolver(svl.GMRES, svl.LU, 60);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.VTK);
	myProblem.saveVelocity(True);
	myProblem.savePressure(True);
	myProblem.setVerbose(False);
	
	myProblem.initialize();
	myProblem.InterpolateFromFacesToCells(ExactVelocity);

	ok = myProblem.run();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );

	print( "------------ !!! End of calculation !!! -----------" );

	normL2 = myProblem.ErrorVelocity(ExactVelocity)
	sizeMesh = M.getNumberOfCells()
	myProblem.terminate();
	normL2faces = normL2[0]
	normL2cellsx = normL2[1]
	normL2cellsy = normL2[2]
	return [sizeMesh, normL2faces, normL2cellsx, normL2cellsy]

if __name__ == """__main__""":
	NormL2faces = []
	NormL2cellsx = []
	NormL2cellsy = []
	sizeMesh = []
	for i in range(5):
		N = 2**i
		result = EulerBarotropicStaggered_2DCylinderDeflection(N)
		print("mesh size = ", result[0], "L2 error on velocity at infinity at FACES= ", result[1], "L2 error on velocity at infinity at CELLS= ", result[2])
		sizeMesh.append(result[0])
		NormL2faces.append(result[1])
		NormL2cellsx.append(result[2])
		NormL2cellsy.append(result[3])
	plt.figure()
	plt.loglog(sizeMesh, NormL2cellsx,label = "L2 error on velocity X component at CELLS")
	plt.loglog(sizeMesh, NormL2cellsy,label = "L2 error on velocity Y component at CELLS")
	plt.loglog(sizeMesh, NormL2faces,label = "L2 error on velocity infty at FACES")
	plt.loglog(sizeMesh, sizeMesh,label = "order 1")
	plt.legend()
	plt.title("error on the stationnary velocity")
	plt.savefig("error on the stationnary velocity")
