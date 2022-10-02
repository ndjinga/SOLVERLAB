#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl
import matplotlib.pyplot as plt

def IsothermalSinglePhase_1DChannelGravity():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh" );
	xinf = 0 ;
	xsup = 1;
	nx=100


    # physical parameters
	gravite=[10];

	myProblem = svl.IsothermalSinglePhase(svl.Liquid,svl.around155bars600K,spaceDim,False);
	nVar =  myProblem.getNumberOfVariables();

	# Prepare for the boundary conditions
	inletVelocity = 1
	outletPressure = 155.e5

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;
	initialVelocity = 1
	initialPressure = 155.e5

	# constant vector
	VV_Constant[0] = initialPressure ;
	VV_Constant[1] = initialVelocity;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Inlet","Outlet");

    # set the boundary conditions
	myProblem.setInletBoundaryCondition( "Inlet",[inletVelocity])
	myProblem.setOutletBoundaryCondition("Outlet",outletPressure,[xsup,0,0]);

    # set physical parameters
	myProblem.setGravity(gravite);

    # set the numerical method
	myProblem.setNumericalScheme(svl.staggered, svl.Implicit);
    
    # name of result file
	fileName = "1DChannelGravityStaggered_Incompressible";

    # simulation parameters 
	MaxNbOfTimeStep = 100 ;
	freqSave = 100;
	cfl = 1;
	maxTime = 5000;
	precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
 
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

    #Postprocessing
	dx = (xsup-xinf)/nx
	x  = [ xinf+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)

	fig, ([axVelocity, axPressure]) = plt.subplots(1, 2,sharex=True, figsize=(10,10))
	fig.suptitle('Implicit pstaggered scheme for isothermal Euler equations \n 1D channel under gravity')

	# Extract density
	myEOS = myProblem.getIncompressibleEOS(0)## Needed to retrieve density
	density  = myEOS.getDensity( outletPressure,  myProblem.getReferenceTemperature() )

	# Build exact solution
	ExactVelocityArray = [inletVelocity]*nx
	ExactPressureArray = [outletPressure + density*gravite[0]*(x[i] - xsup) for i in range(nx)]

	axVelocity.plot( x , ExactVelocityArray, label='Exact stationary velocity')
	axVelocity.set(xlabel='x (m)', ylabel='Velocity')
	axVelocity.set_xlim(xinf,xsup)
	axPressure.plot(x, ExactPressureArray, label='Exact stationary pressure')
	axPressure.set(xlabel='x (m)', ylabel='Pressure')
	axPressure.set_xlim(xinf,xsup)
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

	myPressureField = myProblem.getPressureField()
	pressureArray=myPressureField.getFieldValues()
	myVelocityField = myProblem.getVelocityField()
	velocityArray=myVelocityField.getFieldValues()

	timeStep=myProblem.getNbTimeStep()#Final time step
	axPressure.set_ylim(0.999*min(pressureArray), 1.001*max(pressureArray) )
	axPressure.plot(x, pressureArray,  label='Numerical pressure at time step '+str(timeStep))
	axVelocity.set_ylim(0.999*min(velocityArray), 1.001*max(velocityArray) )
	axVelocity.plot(x, velocityArray,  label='Numerical velocity at time step '+str(timeStep))
	axPressure.legend()
	axVelocity.legend()
	plt.savefig(fileName+".png")

	error_velocity = max( [ abs(ExactVelocityArray[i] - velocityArray[i]) for i in range(nx)] )
	max_velocity = max( [ abs(ExactVelocityArray[i] ) for i in range(nx)] )
	error_pressure = max( [ abs(ExactPressureArray[i] - pressureArray[i]) for i in range(nx)] )
	max_pressure = max( [ abs(ExactPressureArray[i] ) for i in range(nx)] )
	print("Relative error velocity ", error_velocity/max_velocity)
	print("Relative error pressure ", error_pressure/max_pressure)
	assert error_velocity/max_velocity<1e-5
	assert error_pressure/max_pressure<1e-5

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    IsothermalSinglePhase_1DChannelGravity()
