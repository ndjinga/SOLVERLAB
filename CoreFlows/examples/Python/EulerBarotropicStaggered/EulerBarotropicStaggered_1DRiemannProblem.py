#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np


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
	myProblem = svl.EulerBarotropicStaggered(Gas, around1bar300K, spaceDim )

    # Prepare for the initial condition

	print("Building initial data " ); 
		
	def initialPressure(x):
		if x < discontinuity:
			return 6
		elif discontinuity < x:
			return 6

	def initialVelocity(x): # in order to compute exacte solution at the wall bond cond
		if x < xinf:
			return -2
		elif xinf <= x:
			return 2

	def initialVelocityForPb(x): # in order to test th wall boundary cond
		if x < discontinuity:
			return 2
		elif discontinuity <= x:
			return 2

	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1); 
	wallPressureMap = {};
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
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x()) ;
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x());
			Velocity0[j] = initialVelocityForPb(Fj.x())
		elif (Fj.getNumberOfCells()==1):
			for idim in range(spaceDim):
				if vec_normal_sigma[idim] < 0:	
					vec_normal_sigma[idim] = -vec_normal_sigma[idim]
			myProblem.setOrientation(j,vec_normal_sigma)
			if ( j== 0 ): 
				myProblem.setWallBoundIndex(j) 
				wallVelocityMap[j] = 0
			else :
				myProblem.setSteggerBoundIndex(j) 
				wallVelocityMap[j] =initialVelocityForPb(Fj.x()) ;
				wallPressureMap[j] = initialPressure(Fj.x()) ;
			

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
    
    # name of result file
	fileName = "1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 200;
	freqSave = 5;
	cfl = 0.4 
	maxTime = 20;
	precision = 1e-6;

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

	dt = myProblem.getTimeStep()
	Tmax = myProblem.getTime();
	myProblem.terminate();
	time = 0
	i=0
	if not os.path.exists("EulerBarotropicStaggered_"+fileName):
		os.mkdir("EulerBarotropicStaggered_"+fileName)
	while time < Tmax:
		velocitydata = pd.read_csv("EulerBarotropicStaggered_"+fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		pressuredata = pd.read_csv("EulerBarotropicStaggered_"+fileName + "_Pressure_" + str(i)+ ".csv", sep='\s+')
		pressuredata.columns =['x','pressure', 'index']
			
		plt.figure()
		plt.subplot(121)
		plt.plot(pressuredata['x'], pressuredata['pressure'],  label = "pressure results")
		plt.legend()
		plt.subplot(122)
		plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
		plt.legend()
		plt.title("Data at time step"+str(i))
		plt.savefig("EulerBarotropicStaggered_"+fileName + "/Data at time step"+str(i))
		i+=freqSave
		time += freqSave*dt

	return ok

if __name__ == """__main__""":
    EulerBarotropicStaggered_1DRiemannProblem()
