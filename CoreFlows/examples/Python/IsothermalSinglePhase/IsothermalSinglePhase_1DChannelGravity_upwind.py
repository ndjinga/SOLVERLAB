#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl
import matplotlib.pyplot as plt

def IsothermalSinglePhase_1DChannelGravity():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1;
	nx=10;


    # physical parameters
	gravite=[10];

	myProblem = svl.IsothermalSinglePhase(svl.Liquid,svl.around1bar300K,spaceDim);#,False
	nVar =  myProblem.getNumberOfVariables();

	# Prepare for the boundary conditions
	inletVelocity = 1
	outletPressure = 1.e5

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;
	initialVelocity = 1
	initialPressure = 1.e5

	# constant vector
	VV_Constant[0] = initialPressure ;
	VV_Constant[1] = initialVelocity;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Inlet","Outlet");

    # set the boundary conditions
	#myProblem.setNeumannBoundaryCondition("Inlet")
	#myProblem.setNeumannBoundaryCondition("Outlet");
	#myProblem.setWallBoundaryCondition("Inlet")
	#myProblem.setWallBoundaryCondition("Outlet");
	myProblem.setInletBoundaryCondition("Inlet",[inletVelocity])
	myProblem.setOutletBoundaryCondition("Outlet",outletPressure);

    # set physical parameters
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(svl.upwind, svl.Implicit);
	#myProblem.setLinearSolver(svl.GMRES, svl.NOPC);
    
    # name of result file
	fileName = "1DChannelGravityUpwind";

    # simulation parameters 
	MaxNbOfTimeStep = 100 ;
	freqSave = 100;
	cfl = 10;
	maxTime = 500000000;
	precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
 
    # evolution
	myProblem.initialize();

    #Postprocessing
	dx = (xsup-xinf)/nx
	x  = [ i*dx for i in range(nx)]   # array of cell center (1D mesh)

	myPressureField = myProblem.getPressureField()
	pressureArray=myPressureField.getFieldValues()
	myVelocityField = myProblem.getVelocityField()
	velocityArray=myVelocityField.getFieldValues()

	fig, ([axVelocity, axPressure]) = plt.subplots(1, 2,sharex=True, figsize=(10,10))
	fig.suptitle('Implicit upwind scheme for isothermal Euler equations')
	#axVelocity.plot([xinf+0.5*dx + i*dx for i in range(nx)], velocityArray, label='Initial velocity time step 0')
	axVelocity.set(xlabel='x (m)', ylabel='Velocity')
	axVelocity.set_xlim(xinf,xsup)
	axVelocity.set_ylim( 0.999*initialVelocity, 1.001*initialVelocity )
	axVelocity.legend()
	#axPressure.plot([xinf+0.5*dx + i*dx for i in range(nx)], pressureArray, label='Initial pressure time step 0')
	axPressure.set(xlabel='x (m)', ylabel='Pressure')
	axPressure.set_xlim(xinf,xsup)
	axPressure.set_ylim(0.999999*initialPressure, 1.000001*initialPressure )
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
	axPressure.set_ylim(min(pressureArray), max(pressureArray) )
	axPressure.plot(x, pressureArray,  label='Numerical pressure at time step '+str(timeStep))
	axVelocity.set_ylim(min(velocityArray), max(velocityArray) )
	axVelocity.plot(x, velocityArray,  label='Numerical velocity at time step '+str(timeStep))
	axPressure.legend()
	axVelocity.legend()
	plt.savefig(fileName+".png")

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    IsothermalSinglePhase_1DChannelGravity()
