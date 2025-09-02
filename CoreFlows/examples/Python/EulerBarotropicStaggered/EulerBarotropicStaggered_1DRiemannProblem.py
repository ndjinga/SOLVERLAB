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
	xsup= 1
	nx=400	;
	M=svl.Mesh(xinf,xsup,nx)
	discontinuity=(xinf+xsup)/2 + 0.75/nx

    # set the limit field for each boundary
	a=1.0
	gamma = 2
	myProblem = svl.EulerBarotropicStaggered(svl.GasStaggered, svl.around1bar300K, a, gamma, spaceDim ); 

    # Prepare for the initial condition

	print("Building initial data " ); 
	initialDensity_Left = 3; #12
	initialDensity_Right = 1 #1

	initialVelocity_Left = 5 ;#1.5
	initialVelocity_Right = 8 # -3

	def initialDensity(x):
		if x < discontinuity:
			return initialDensity_Left
		elif discontinuity < x:
			return initialDensity_Right

	def initialVelocity(x): 
		if x < discontinuity:
			return initialVelocity_Left
		elif discontinuity <= x:
			return initialVelocity_Right
	

	Density0 = svl.Field("Density", svl.CELLS, M, 1);
	Momentum0 = svl.Field("Momentum", svl.FACES, M, 1); 
	wallDensityMap = {};
	wallVelocityMap = {}; 
	wallVelocityVector = np.zeros(1)
	vec_normal_sigma = np.zeros(1)

	for j in range( M.getNumberOfFaces() ):
		Fj = M.getFace(j);
		idCells = Fj.getCellsId()
		Ctemp1 = M.getCell(idCells[0])
		for l in range( Ctemp1.getNumberOfFaces()) :
			if (j == Ctemp1.getFacesId()[l]):
				for idim in range(spaceDim):
					vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
		
		Density0[idCells[0]] = initialDensity(Ctemp1.x()) ;
		if(Fj.getNumberOfCells()==2):	
			myProblem.setOrientation(j,vec_normal_sigma)
			myProblem.setInteriorIndex(j);
			Ctemp2 = M.getCell(idCells[1]);
			Density0[idCells[1]] = initialDensity(Ctemp2.x());
			Momentum0[j] = initialVelocity(Fj.x() ) * (initialDensity(Ctemp1.x()) + initialDensity(Ctemp2.x()))/2.0 #TODO times normal sigma
	
		elif (Fj.getNumberOfCells()==1):
			Momentum0[j] = initialVelocity(Fj.x())* (initialDensity(Ctemp1.x()) + initialDensity(Fj.x()))/2.0 #TODO times normal sigma
			# Since we plot the values of the veloccity at the faces every velocity if oriented to the right n_sigma =1
			for idim in range(spaceDim):
				if vec_normal_sigma[idim] < 0:	
					vec_normal_sigma[idim] = -vec_normal_sigma[idim]
			myProblem.setOrientation(j,vec_normal_sigma)
			myProblem.setSteggerBoundIndex(j) 
			wallVelocityMap[j] = initialVelocity(Fj.x())
			wallDensityMap[j] = initialDensity(Fj.x()) ;
			wallVelocityVector[0] = wallVelocityMap[j] 
			myProblem.setboundaryVelocityVector(j, wallVelocityVector)
		
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    # set the numerical method
	myProblem.setTimeScheme(svl.Explicit);
    
    # name of result file
	fileName = "EulerBarotropicStaggered_1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 100000;
	freqSave = 200;
	cfl = 0.8
	maxTime = 0.02
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
	
	
	if not os.path.exists(fileName):
		os.mkdir(fileName)

	i = freqSave
	while time[i] <= Tmax:
		velocitydata = pd.read_csv(fileName + "_Velocity_" + str(i)+ ".csv", sep='\s+')
		velocitydata.columns =['x','velocity', 'index']
		Densitydata = pd.read_csv(fileName + "_Density_" + str(i)+ ".csv", sep='\s+')
		Densitydata.columns =['x','Density', 'index']
		
		myEOS = myProblem.getBarotropicEOS(0)	
		initialPressure_Left  = myEOS.getPressure(initialDensity_Left)
		initialPressure_Right  = myEOS.getPressure(initialDensity_Right)  
		a = myEOS.constante("a");
		gamma = myEOS.constante("gamma");
		exactDensity, exactVelocity = exact_rs_euler_barotropic.exact_sol_Riemann_problem(xinf, xsup, time[i], gamma, a , [initialPressure_Left , initialVelocity_Left ], [ initialPressure_Right, initialVelocity_Right ], (xinf+xsup)/2, nx)

		plt.figure()
		plt.subplot(121)
		plt.plot(Densitydata['x'], exactDensity,  label = "exact Density")
		plt.plot(Densitydata['x'], Densitydata['Density'],  label = "Density results")
		plt.legend()
		plt.subplot(122) 
		plt.plot(Densitydata['x'], exactVelocity,label = "exact velocity")
		plt.plot(velocitydata['x'], velocitydata['velocity'],  label = "velocity results")
		plt.legend()
		plt.title("Data at time step"+str(i)+"t ="+str(time[i]))
		plt.savefig(fileName + "/Data at time step"+str(i))
		i+= freqSave

	# Print only at final time
	timestep_index = -1  
	Density_file = f"{fileName}_Density_{str(len(time)-1)}.csv"
	velocity_file = f"{fileName}_Velocity_{str(len(time) -1)}.csv"
	try:
		velocitydata = pd.read_csv(velocity_file, sep=r'\s+', engine='python', names=['x', 'velocity', 'index'], skiprows=1)
		Densitydata = pd.read_csv(Density_file, sep=r'\s+', engine='python', names=['x', 'Density', 'index'], skiprows=1)
	except FileNotFoundError:
		print(f"Fichiers non trouvés : {velocity_file} ou {Density_file}")
		exit()

	# Trace
	plt.figure(figsize=(12, 5))
	myEOS = myProblem.getBarotropicEOS(0)	
	initialPressure_Left  = myEOS.getPressure(initialDensity_Left)
	initialPressure_Right  = myEOS.getPressure(initialDensity_Right)  
	a = myEOS.constante("a");
	gamma = myEOS.constante("gamma");
	exactDensity, exactVelocity = exact_rs_euler_barotropic.exact_sol_Riemann_problem(xinf, xsup, Tmax, gamma, a , [initialPressure_Left , initialVelocity_Left ], [ initialPressure_Right, initialVelocity_Right ], (xinf+xsup)/2, nx)

	plt.subplot(1, 2, 1)
	
	plt.plot(Densitydata['x'], Densitydata['Density'], label="Density results", color='tab:blue')
	plt.plot(Densitydata['x'], exactDensity, label="Exact density", color='tab:green')
	pyplot.xlabel("x",12,labelpad=1)          
	pyplot.ylabel(r"$\rho$",12,labelpad=1) 
	plt.title("Density profil")
	plt.grid(True)
	plt.legend()

	plt.subplot(1, 2, 2)
	plt.plot(velocitydata['x'], velocitydata['velocity'], label="Velocity results", color='tab:blue')
	plt.plot(Densitydata['x'], exactVelocity, label="Exact velocity", color='tab:green')
	pyplot.xlabel("x",12,labelpad=1)          
	pyplot.ylabel(r"$u$",12,labelpad=1) 
	plt.title("Velocity profil")
	plt.grid(True)
	plt.legend()

	plt.suptitle(f"t = {time[-1]}")
	plt.tight_layout()
	plt.savefig(f"{fileName}_Final_results.png", dpi=300)
	plt.show()


	
	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    EulerBarotropicStaggered_1DRiemannProblem()

