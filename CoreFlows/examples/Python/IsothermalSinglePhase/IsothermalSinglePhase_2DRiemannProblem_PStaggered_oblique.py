#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl
import matplotlib.pyplot as plt
import exact_rs_stiffenedgas
import VTK_routines
from math import sqrt

def IsothermalSinglePhase_2DRiemannProblem_Staggered_symmetric():

	spaceDim = 2;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1;
	yinf = 0 ;
	ysup=1;
	nx=20;
	ny=20;
	discontinuity=(xinf+xsup)/2
	M=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	eps=1e-6
	M.setGroupAtPlan(xsup,0,eps,"RightBoundary")
	M.setGroupAtPlan(xinf,0,eps,"LeftBoundary")
	M.setGroupAtPlan(ysup,1,eps,"TopBoundary")
	M.setGroupAtPlan(yinf,1,eps,"BottomBoundary")

    # Prepare initial data
	initialVelocityX_Left=-1;
	initialVelocityY_Left= 1;
	initialPressure_Left=1e5;
	initialVelocityX_Right=1;
	initialVelocityY_Right=1;
	initialPressure_Right=1e5;

	myProblem = svl.IsothermalSinglePhase(svl.Gas,svl.around1bar300K,spaceDim,True);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Left  = svl.Vector(nVar)
	VV_Right = svl.Vector(nVar)
	
	# left and right constant vectors		
	VV_Left[0] = initialPressure_Left;
	VV_Left[1] = initialVelocityX_Left;
	VV_Left[2] = initialVelocityY_Left;
	VV_Right[0] = initialPressure_Right;
	VV_Right[1] = initialVelocityX_Right;
	VV_Right[2] = initialVelocityY_Right;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary");
	myProblem.setNeumannBoundaryCondition("RightBoundary");
	myProblem.setNeumannBoundaryCondition("TopBoundary");
	myProblem.setNeumannBoundaryCondition("BottomBoundary");

    # set the numerical method
	myProblem.setNumericalScheme(svl.staggered, svl.Implicit);
    
    # name of result file
	fileName = "2DRiemannProblem_staggered_oblique";

    # simulation parameters 
	MaxNbOfTimeStep = 5 ;
	freqSave = 1;
	cfl = 1;
	maxTime = 500;
	precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass
	
	myProblem.setLinearSolver(svl.GMRES, svl.LU)
	myProblem.setNewtonSolver(precision,50, svl.Newton_SOLVERLAB)
 
    # Initialisation
	myProblem.initialize();

    # Postprocessing
	dx=(xsup-xinf)/nx
	x=[xinf+0.5*dx + i*dx for i in range(nx+1)]   # array of cell center (1D mesh)

	myPressureField = myProblem.getPressureField()
	pressureArray=VTK_routines.Extract_field_data_over_line_to_numpyArray(myPressureField,[xinf,(yinf+ysup)/2,0],[xsup,(yinf+ysup)/2,0], nx)
	myVelocityField = myProblem.getVelocityXField()
	velocityArray=VTK_routines.Extract_field_data_over_line_to_numpyArray(myVelocityField,[xinf,(yinf+ysup)/2,0],[xsup,(yinf+ysup)/2,0], nx)

	fig, ([axVelocity, axPressure]) = plt.subplots(1, 2,sharex=True, figsize=(10,10))
	fig.suptitle('PStaggered scheme for oblique Riemmann Problem \n 2D isothermal Euler equations')
	axVelocity.plot(x, velocityArray, label='Initial velocity time step 0')
	axVelocity.set(xlabel='x (m)', ylabel='Velocity')
	axVelocity.set_xlim(xinf,xsup)
	axVelocity.set_ylim( 0.999*min(initialVelocityX_Left, initialVelocityX_Right), 1.001*max(initialVelocityX_Left, initialVelocityX_Right) )
	axVelocity.legend()
	axPressure.plot(x, pressureArray, label='Initial pressure time step 0')
	axPressure.set(xlabel='x (m)', ylabel='Pressure')
	axPressure.set_xlim(xinf,xsup)
	axPressure.set_ylim(0.95*min(initialPressure_Left, initialPressure_Right), 1.05*max(initialPressure_Left, initialPressure_Right) )
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axPressure.legend()


    # Evolution
	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );

	# Extract EOS and Riemann problem parameters
	myEOS = myProblem.getStiffenedGasEOS(0)## Needed to retrieve gamma, pinfnity, convert (p,T) to density and (p, rho) to temperature
	initialDensity_Left  = myEOS.getDensity( initialPressure_Left,  myProblem.getReferenceTemperature() )
	initialDensity_Right = myEOS.getDensity( initialPressure_Right, myProblem.getReferenceTemperature() )

	#Determine exact solution
	normVelocityLeft = sqrt( initialVelocityX_Left*initialVelocityX_Left + initialVelocityY_Left*initialVelocityY_Left )
	normVelocityRight = sqrt( initialVelocityX_Right*initialVelocityX_Right + initialVelocityY_Right*initialVelocityY_Right )
	exactDensity, exactVelocity, exactPressure = exact_rs_stiffenedgas.exact_sol_Riemann_problem(xinf, xsup, myProblem.presentTime(), myEOS.constante("gamma"), myEOS.constante("p0"), [ initialDensity_Left, initialVelocityX_Left/normVelocityLeft, initialPressure_Left ], [ initialDensity_Right, initialVelocityX_Right/normVelocityRight, initialPressure_Right ], (xinf+xsup)/2, nx+1)

	myPressureField = myProblem.getPressureField()
	pressureArray=VTK_routines.Extract_field_data_over_line_to_numpyArray(myPressureField,[xinf,(yinf+ysup)/2,0],[xsup,(yinf+ysup)/2,0], nx)
	myVelocityField = myProblem.getVelocityXField()
	velocityArray=VTK_routines.Extract_field_data_over_line_to_numpyArray(myVelocityField,[xinf,(yinf+ysup)/2,0],[xsup,(yinf+ysup)/2,0], nx)

	timeStep=myProblem.getNbTimeStep()#Final time step
	axPressure.plot(x, exactPressure,  label='Exact pressure at time step '+str(timeStep))
	axPressure.plot(x, pressureArray,  label='Numerical pressure at time step '+str(timeStep))
	axVelocity.plot(x, exactVelocity,  label='Exact velocity at time step '+str(timeStep))
	axVelocity.plot(x, velocityArray,  label='Numerical velocity at time step '+str(timeStep))
	axPressure.legend()
	axVelocity.legend()
	plt.savefig(fileName+".png")


	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	IsothermalSinglePhase_2DRiemannProblem_Staggered_symmetric()
