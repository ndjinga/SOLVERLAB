#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np


def WaveStaggered_1DLongTimeLimit():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1
	nx=90;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 + 0.75/nx

    # set the limit field for each boundary
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim, rho, kappa);

    # Prepare for the initial condition

	print("Building initial data " ); 
	initialVelocity_Left=-12;
	initialVelocity_Right=-1;
	
	def initialPressure(x):
		return 10
	def initialVelocity(x):
		if xinf<x < discontinuity:
			return initialVelocity_Left
		elif x <=xinf :
			return initialVelocity_Left
		elif discontinuity < x < xsup:
			return initialVelocity_Right
		elif xsup <= x:
			return initialVelocity_Left


	Pressure0 = svl.Field("pressure", svl.CELLS, M, 1);
	Velocity0 = svl.Field("velocity", svl.FACES, M, 1);
	wallPressureMap = {};
	wallVelocityMap = {}; 
	
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId();
		if(Fj.getNumberOfCells()==2):
			Ctemp1 = M.getCell(idCells[0]);
			Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x()) ;
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x());
			Velocity0[j] = initialVelocity(Fj.x()) ;
		else:
			wallPressureMap[j] = initialPressure(Fj.x()) ;
			wallVelocityMap[j] = initialVelocity(Fj.x()) ;

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);

    
    # name of result file
	fileName = "1DLongTimeLimit";

    # simulation parameters 
	MaxNbOfTimeStep = 10000;
	freqSave = 200;
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
	if not os.path.exists("WaveStaggered_"+fileName):
		os.mkdir("WaveStaggered_"+fileName)
	while time < Tmax:
		velocitydata = pd.read_csv("WaveStaggered_"+fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		pressuredata = pd.read_csv("WaveStaggered_"+fileName + "_Pressure_" + str(i)+ ".csv", sep='\s+')
		pressuredata.columns =['x','pressure', 'index']
	
		plt.figure()
		plt.subplot(121)
		plt.plot(pressuredata['x'], pressuredata['pressure'],  label = "pressure results")
		plt.legend()
		plt.subplot(122)
		plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
		plt.legend()
		plt.title("Data at time step"+str(i))
		plt.savefig("WaveStaggered_"+fileName + "/Data at time step"+str(i))
		i+=freqSave
		time += freqSave*dt

	return ok

if __name__ == """__main__""":
    WaveStaggered_1DLongTimeLimit()