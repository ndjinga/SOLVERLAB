#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
import exact_rs_euler_barotropic


def EulerBarotropicStaggered_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1
	nx=200;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 + 0.75/nx

    # set the limit field for each boundary
	a=1.0
	gamma = 2.0
	myProblem = svl.EulerBarotropicStaggered(svl.GasStaggered, svl.around1bar300K, a, gamma, spaceDim ); 

    # Prepare for the initial condition

	print("Building initial data " ); 
	initialDensity_Left = 2
	initialDensity_Right = 2

	initialVelocity_Left = 1
	initialVelocity_Right = 1

	
	def initialDensity(x):
		if x < discontinuity:
			return initialDensity_Left
		elif discontinuity < x:
			return initialDensity_Right

	def initialVelocityForPb(x): # in order to test th wall boundary cond
		if x < discontinuity:
			return initialVelocity_Left
		elif discontinuity <= x:
			return initialVelocity_Right

	Density0 = svl.Field("Density", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1); 
	wallDensityMap = {};
	wallVelocityMap = {}; 

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
			Density0[idCells[0]] = initialDensity(Ctemp1.x()) ;
			Density0[idCells[1]] = initialDensity(Ctemp2.x());
			Velocity0[j] = initialVelocityForPb(Fj.x())
		elif (Fj.getNumberOfCells()==1):
			for idim in range(spaceDim):
				if vec_normal_sigma[idim] < 0:	
					vec_normal_sigma[idim] = -vec_normal_sigma[idim]
			myProblem.setOrientation(j,vec_normal_sigma)
			if (j == M.getNumberOfFaces()-1 ):
				myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] =0
			else : 
				myProblem.setSteggerBoundIndex(j) 
				wallVelocityMap[j] =initialVelocityForPb(Fj.x())
				wallDensityMap[j] = initialDensity(Fj.x()) ;
			
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
    
    # name of result file
	fileName = "EulerBarotropicStaggered_1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 20000;
	freqSave = 400;
	cfl = 0.5
	maxTime = 0.1;
	precision = 1e-10;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.CSV)
	myProblem.saveVelocity();
	myProblem.savePressure();
	myProblem.setVerbose(False);

    # evolution
	myProblem.initialize();

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );


		
	Tmax = myProblem.getTime();
	time = myProblem.getTimeEvol();

	myProblem.terminate();
	if not os.path.exists(fileName):
		os.mkdir(fileName)
	
	#print only at final time 
	velocitydata = pd.read_csv(fileName + "_Velocity_" + str(len(time) -1)+ ".csv", sep='\s+')
	velocitydata.columns =['x','velocity', 'index']
	Densitydata = pd.read_csv(fileName + "_Pressure_" + str(len(time) -1)+ ".csv", sep='\s+')
	Densitydata.columns =['x','pressure', 'index']
	
	myEOS = myProblem.getBarotropicEOS(0)#
	initialPressure_Right = myEOS.getPressure(initialDensity_Right)
	a = myEOS.constante("a");
	gamma = myEOS.constante("gamma");
	initalVelocityWall = -initialVelocity_Left
	exactDensity, exactVelocity = exact_rs_euler_barotropic.exact_sol_Riemann_problem(xinf, xsup, Tmax, gamma, a , [initialPressure_Right , initialVelocity_Left ], [ initialPressure_Right, initalVelocityWall ], xsup, nx)
	
	plt.figure()
	plt.subplot(121)
	plt.plot(Densitydata['x'], exactDensity,  label = "exact Density")
	plt.plot(Densitydata['x'], Densitydata['pressure'],  label = "pressure results")
	plt.legend()
	plt.subplot(122)
	plt.plot(Densitydata['x'], exactVelocity,label = "exact velocity")
	plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
	plt.legend()
	plt.title("Data at time step"+str(len(time) -1)+"t ="+str(Tmax))
	plt.savefig(fileName + "/Data at time step"+str(len(time) -1))

	return ok

if __name__ == """__main__""":
    EulerBarotropicStaggered_1DRiemannProblem()