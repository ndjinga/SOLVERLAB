#!/usr/bin/env python
# -*-coding:utf-8 -*

import solverlab as svl
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
import exact_rs_stiffenedgas


def EulerBarotropicStaggered_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = -1 ;
	xsup=1
	nx=30;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 #+ 0.75/nx

    # set the limit field for each boundary
	myProblem = svl.EulerBarotropicStaggered(svl.GasStaggered, svl.around1bar300K, spaceDim ); 

    # Prepare for the initial condition

	print("Building initial data " ); 
	initialPressure_Left = 2
	initialPressure_Right =10

	initialVelocity_Left = -3
	initialVelocity_Right = 8

	
	def initialPressure(x):
		if x < discontinuity:
			return initialPressure_Left
		elif discontinuity < x:
			return initialPressure_Right

	def initialVelocityForPb(x): # in order to test th wall boundary cond
		if x < discontinuity:
			return initialVelocity_Left
		elif discontinuity <= x:
			return initialVelocity_Right

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
	fileName = "EulerBarotropicStaggered_1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 120;
	freqSave = 2;
	cfl = 0.5
	maxTime = 20;
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
	def initialPressure(x):
		if x <0:
			return initialPressure_Left
		else :
			return initialPressure_Right

	def ExactPressure(x,t, velocity):
		return initialPressure(x - velocity * t) 

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
	
	#TODO for now pressure law is rho^2 for more general case replac underneath 2 by myEOS.constante("gamma") and 1 by myEOS.constante("p0")
	myProblem.terminate();
	if not os.path.exists(fileName):
		os.mkdir(fileName)
	i=0
	while time[i] < Tmax:
		velocitydata = pd.read_csv(fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		pressuredata = pd.read_csv(fileName + "_Pressure_" + str(i)+ ".csv", sep='\s+')
		pressuredata.columns =['x','pressure', 'index']
		
		#Determine exact solution
		#TODO for now pressure law is rho^2 for more general case replac underneath 2 by myEOS.constante("gamma") and 0 by myEOS.constante("p0")
		#myEOS = myProblem.getStiffenedGasEOS(0)## Needed to retrieve gamma, pinfnity, convert (p,T) to density and (p, rho) to temperature
		#initialDensity_Left  = myEOS.getDensity( initialPressure_Left,  myProblem.getReferenceTemperature() )
		#initialDensity_Right = myEOS.getDensity( initialPressure_Right, myProblem.getReferenceTemperature() )
		#exactDensity, exactVelocity, exactPressure = exact_rs_stiffenedgas.exact_sol_Riemann_problem(xinf, xsup, time, 2 , 0, [ initialDensity_Left, initialVelocity_Left, initialPressure_Left ], [ initialDensity_Right, initialVelocity_Right, initialPressure_Right ], (xinf+xsup)/2, nx)

		exactDensity = np.zeros(nx)
		exactVelocity = np.zeros(nx)
		ExactVelocity = ExactVelocityFromRiemannProblem(initialVelocity_Left,initialVelocity_Right)
		for j in range(nx):
			exactVelocity[j] = ExactVelocity(xinf + j*(xsup - xinf)/nx,time[i])
		for j in range(nx):
			exactDensity[j] = ExactPressure(xinf + j*(xsup - xinf)/nx + (xsup - xinf)/(2*nx),time[i], exactVelocity[j])
		
			
		plt.figure()
		plt.subplot(121)
		plt.plot(pressuredata['x'], exactDensity,  label = "exact pressure")
		plt.plot(pressuredata['x'], pressuredata['pressure'],  label = "pressure results")
		plt.legend()
		plt.subplot(122)
		plt.plot(pressuredata['x'], exactVelocity,label = "exact velocity")
		plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
		plt.legend()
		plt.title("Data at time step"+str(i))
		plt.savefig(fileName + "/Data at time step"+str(i))
		i+= freqSave
	return ok

if __name__ == """__main__""":
    EulerBarotropicStaggered_1DRiemannProblem()
