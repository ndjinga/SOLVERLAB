#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np


def WaveStaggered_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1
	nx=400;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 + 0.75/nx

    # set the limit field for each boundary
	kappa = 1;
	rho = 1;
	c = math.sqrt(kappa/rho)
	myProblem = svl.WaveStaggered(spaceDim, rho, kappa);

    # Prepare for the initial condition
	Pressure_Left =svl.Vector(1);
	Pressure_Right =svl.Vector(1);
	Velocity_Left =svl.Vector(1);
	Velocity_Right =svl.Vector(1);

	initialVelocity_Left=3;
	initialPressure_Left=-1;

	initialVelocity_Right=1;
	initialPressure_Right=3;
	
	# left and right constant vectors		
	Pressure_Left[0] = initialPressure_Left;
	Pressure_Right[0] = initialPressure_Right;
	Velocity_Left[0] = initialVelocity_Left;
	Velocity_Right[0] = initialVelocity_Right;


    #Initial field creation
	print("Building initial data " ); 
	direction = 0;
	myProblem.setInitialFieldStepFunction(M,Pressure_Left,Pressure_Right,discontinuity, direction, svl.CELLS);
	myProblem.setInitialFieldStepFunction(M,Velocity_Left,Velocity_Right,discontinuity, direction, svl.FACES);

	# set the boundary conditions
	def initialPressure(x):
		if x < discontinuity:
			return initialPressure_Left
		else :
			return initialPressure_Right

	def initialVelocity(x):
		if x < discontinuity:
			return initialVelocity_Left
		else:
			return initialVelocity_Right
	
	# Define the exact solution of the 1d Problem 
	def ExactPressure(x,t):
		return (initialPressure(x - c * t) + initialPressure(x + c * t))/2.0 + (initialVelocity(x-c*t) -initialVelocity(x+c*t))/(2*rho*c)
	def ExactVelocity(x,t):
		return (initialVelocity(x - c * t) + initialVelocity(x + c * t))/2.0 + rho*c*(initialPressure(x-c*t) -initialPressure(x+c*t))/2.0

	wallPressureMap = {};
	wallVelocityMap = {}; 
	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		if (Fj.getNumberOfCells()==1):
			wallPressureMap[j] = initialPressure(Fj.x()) ;
			wallVelocityMap[j] = initialVelocity(Fj.x()) ;

	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);

    
    # name of result file
	fileName = "1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 700 ;
	freqSave = 1;
	cfl = 0.4 
	maxTime = 500;
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
	myProblem.terminate();
	time = 0
	i=0
	if not os.path.exists("WaveStaggered_"+fileName):
		os.mkdir("WaveStaggered_"+fileName)
	for t in range(MaxNbOfTimeStep):
		velocitydata = pd.read_csv("WaveStaggered_"+fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		pressuredata = pd.read_csv("WaveStaggered_"+fileName + "_Pressure_" + str(i)+ ".csv", sep='\s+')
		pressuredata.columns =['x','pressure', 'index']
		
		pressure = np.zeros(nx)
		velocity = np.zeros(nx+1)
		for j in range(nx):
			pressure[j] = ExactPressure(xinf + j*(xsup - xinf)/nx + (xsup - xinf)/(2*nx),time)
		for j in range(nx+1):
			velocity[j] = ExactVelocity(xinf + j*(xsup - xinf)/nx,time)
			
		plt.figure()
		plt.subplot(121)
		plt.plot(pressuredata['x'], pressure, 'k-',  label = "exact pressure")
		plt.plot(pressuredata['x'], pressuredata['pressure'],  label = "pressure results")
		plt.legend()
		plt.subplot(122)
		plt.plot(velocitydata['x'], velocitydata['velocity'], 'k-',  label = "velocity results")
		plt.plot(velocitydata['x'], velocity, label = "exact velocity")	
		plt.legend()
		plt.title("Data at time step"+str(i))
		plt.savefig("WaveStaggered_"+fileName + "/Data at time step"+str(i))
		i+=1
		time += dt

	return ok

if __name__ == """__main__""":
    WaveStaggered_1DRiemannProblem()
