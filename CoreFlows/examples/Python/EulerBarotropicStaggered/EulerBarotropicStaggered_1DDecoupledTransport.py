#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
import exact_rs_stiffenedgas

# Impose u_0 = cst, gradp = 0, c=1 -> we have a transport equation on rho
# Impose rho_0 = 1.0, div(rhou) =0, gradp = 0, u_0 discontinuous, c=1 -> Burgers equation on u
def EulerBarotropicStaggered_1DDecoupledTransport():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = -1 ;
	xsup=1
	nx=200;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 #+ 0.75/nx

    # set the limit field for each boundary
	a=1.0
	gamma = 2.0
	myProblem = svl.EulerBarotropicStaggered(svl.GasStaggered, svl.around1bar300K, a, gamma, spaceDim ); 

    # Prepare for the initial condition

	print("Building initial data " ); 
	initialDensity_Left = 1
	initialDensity_Right = 1

	initialVelocity_Left = 1.5
	initialVelocity_Right = -3

	
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
			myProblem.setSteggerBoundIndex(j) 
			wallVelocityMap[j] =initialVelocityForPb(Fj.x()) ;
			wallDensityMap[j] = initialDensity(Fj.x()) ;
			

	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
    
    # name of result file
	fileName = "EulerBarotropicStaggered_1DDecoupledTransport";

    # simulation parameters 
	MaxNbOfTimeStep = 100;
	freqSave = 1;
	cfl = 0.99
	maxTime = 0.05;
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


	# Define the exact solution of the 1d Problem 
	def initialDensity(x):
		if x <0:
			return initialDensity_Left
		else :
			return initialDensity_Right

	def ExactDensity(x,t, velocity):
		return initialDensity(x - velocity * t) 

	def ExactVelocityFromRiemannProblem(u1,u2):
		def u(x,t):
			if t== 0:
					if x <0:
						return u1
					else :
						return u2
			else :
				if u1 <= u2:
					if  x/t <= u1:
						return u1
					elif u1 <= x/t <= u2:
						return x/t
					elif u2 <= x/t: 
						return u2
				if u2 < u1:
					if x/t<(u1+u2)/2.0:
						return u1
					elif (u1+u2)/2.0 < x/t:
						return u2
		return u
		
	Tmax = myProblem.getTime();
	time = myProblem.getTimeEvol();
	
	myProblem.terminate();
	if not os.path.exists(fileName):
		os.mkdir(fileName)
	i=0
	while time[i] < Tmax:
		velocitydata = pd.read_csv(fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		Densitydata = pd.read_csv(fileName + "_Pressure_" + str(i)+ ".csv", sep='\s+')
		Densitydata.columns =['x','pressure', 'index']
	
		exactDensity = np.zeros(nx)
		exactVelocity = np.zeros(nx)
		ExactVelocity = ExactVelocityFromRiemannProblem(initialVelocity_Left,initialVelocity_Right)
		for j in range(nx):
			exactVelocity[j] = ExactVelocity(xinf + j*(xsup - xinf)/nx,time[i])
		for j in range(nx):
			exactDensity[j] = ExactDensity(xinf + j*(xsup - xinf)/nx + (xsup - xinf)/(2*nx),time[i], initialVelocity_Left)
		
			
		plt.figure()
		plt.subplot(121)
		plt.plot(Densitydata['x'], exactDensity,  label = "exact Density")
		plt.plot(Densitydata['x'], Densitydata['pressure'],  label = "pressure results")
		plt.legend()
		plt.subplot(122)
		plt.plot(Densitydata['x'], exactVelocity,label = "exact velocity")
		plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
		plt.legend()
		plt.title("Data at time step"+str(i))
		plt.savefig(fileName + "/Data at time step"+str(i))
		i+= freqSave
	return ok

if __name__ == """__main__""":
    EulerBarotropicStaggered_1DDecoupledTransport()
