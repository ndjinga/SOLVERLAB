#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl
import matplotlib.pyplot as plt

def IsothermalSinglePhase_1DRiemannProblem():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=100;
	discontinuity=(xinf+xsup)/2
	M=svl.Mesh(xinf,xsup,nx)
	eps=1e-6
	M.setGroupAtPlan(xsup,0,eps,"RightBoundary")
	M.setGroupAtPlan(xinf,0,eps,"LeftBoundary")

    # set the limit field for each boundary
	initialVelocity_Left=1;
	initialPressure_Left=155e5;

	initialVelocity_Right=1;
	initialPressure_Right=155.001e5;

	myProblem = svl.IsothermalSinglePhase(svl.Liquid,svl.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Left  =svl.Vector(nVar)
	VV_Right =svl.Vector(nVar)
	
	# left and right constant vectors		
	VV_Left[0] = initialPressure_Left;
	VV_Left[1] = initialVelocity_Left;

	VV_Right[0] = initialPressure_Right;
	VV_Right[1] = initialVelocity_Right;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary");
	myProblem.setNeumannBoundaryCondition("RightBoundary");

    # set the numerical method
	myProblem.setNumericalScheme(svl.upwind, svl.Explicit);
    
    # name of result file
	fileName = "1DRiemannProblem";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
 
    # evolution
	myProblem.initialize();

    #Postprocessing
	dx=(xsup-xinf)/nx
	x=[ i*dx for i in range(nx)]   # array of cell center (1D mesh)

	myPressureField = myProblem.getPressureField()
	pressureArray=myPressureField.getFieldValues()
	myVelocityField = myProblem.getVelocityField()
	velocityArray=myVelocityField.getFieldValues()

	fig, ([axVelocity, axPressure]) = plt.subplots(1, 2,sharex=True, figsize=(10,10))
	fig.suptitle('Explicit Upwind scheme for isothermal Euler equations')
	axVelocity.plot([xinf+0.5*dx + i*dx for i in range(nx)], velocityArray, label='Velocity time step 0')
	axVelocity.set(xlabel='x (m)', ylabel='Velocity')
	axVelocity.set_xlim(xinf,xsup)
	axVelocity.set_ylim( 0.999*min(initialVelocity_Left, initialVelocity_Right), 1.001*max(initialVelocity_Left, initialVelocity_Right) )
	axVelocity.legend()
	axPressure.plot([xinf+0.5*dx + i*dx for i in range(nx)], pressureArray, label='Pressure  time step 0')
	axPressure.set(xlabel='x (m)', ylabel='Pressure')
	axPressure.set_xlim(xinf,xsup)
	axPressure.set_ylim(0.999999*min(initialPressure_Left, initialPressure_Right), 1.000001*max(initialPressure_Left, initialPressure_Right) )
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axPressure.legend()

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );

	myPressureField = myProblem.getPressureField()
	pressureArray=myPressureField.getFieldValues()
	myVelocityField = myProblem.getVelocityField()
	velocityArray=myVelocityField.getFieldValues()
	timeStep=myProblem.getNbTimeStep()#Final time step
	axPressure.plot(x, pressureArray,  label='Pressure time step '+str(timeStep))
	axVelocity.plot(x, velocityArray,  label='Velocity time step '+str(timeStep))
	axPressure.legend()
	axVelocity.legend()
	plt.savefig(fileName+".png")

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    IsothermalSinglePhase_1DRiemannProblem()