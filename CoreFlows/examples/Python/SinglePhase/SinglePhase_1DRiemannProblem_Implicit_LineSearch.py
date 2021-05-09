#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm
import matplotlib.pyplot as plt
import VTK_routines

def SinglePhase_1DRiemannProblem_Implicit_LineSearch():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;
	discontinuity=(xinf+xsup)/2
	M=cm.Mesh(xinf,xsup,nx)
	eps=1e-6
	M.setGroupAtPlan(xsup,0,eps,"RightBoundary")
	M.setGroupAtPlan(xinf,0,eps,"LeftBoundary")

    # Prepare initial data
	initialVelocity_Left=1;
	initialTemperature_Left=565;
	initialPressure_Left=155e5;
	initialVelocity_Right=1;
	initialTemperature_Right=565;
	initialPressure_Right=1e5;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Left  = cm.Vector(nVar)
	VV_Right = cm.Vector(nVar)
	
	# left and right constant vectors		
	VV_Left[0] = initialPressure_Left;
	VV_Left[1] = initialVelocity_Left;
	VV_Left[2] = initialTemperature_Left ;
	VV_Right[0] = initialPressure_Right;
	VV_Right[1] = initialVelocity_Right;
	VV_Right[2] = initialTemperature_Right ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

    #Postprocessing
	plt.xlabel('x')
	plt.ylabel('Pressure')
	plt.xlim(xinf,xsup)
	plt.ylim( initialPressure_Right - 0.1*(initialPressure_Left-initialPressure_Right), initialPressure_Left +  0.5*(initialPressure_Left-initialPressure_Right) )
	plt.title('Solving Riemann problem for Euler equations\n with Finite volume schemes method')
	dx=(xsup-xinf)/nx
	x=[ i*dx for i in range(nx+1)]   # array of cell center (1D mesh)
	myDensityField = myProblem.getDensityField()
	#densityArray=VTK_routines.Extract_field_data_over_line_to_numpyArray(myDensityField(), [xinf,0,0], [xsup,0,0],nx)
	#line1, = plt.plot(x, densityArray,  label='Time step 0')

    # set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary");
	myProblem.setNeumannBoundaryCondition("RightBoundary");

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
    
    # name of result file
	fileName = "1DRiemannProblem_Implicit_LineSearch";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 1;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(cf.CSV)
	myProblem.saveConservativeField(True);
	
	myProblem.setLinearSolver(cf.GMRES, cf.ILU)
	myProblem.setNewtonSolver(precision,20, cf.Newton_PETSC_LINESEARCH)
	myProblem.usePrimitiveVarsInNewton(False)
 
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

	momentumField=myProblem.getMomentumField()
	print("Momentum X in first cell ", momentumField[0,0])
	print("Momentum Y in first cell ", momentumField[0,1])

	densityField=myProblem.getDensityField()
	print("density X in first cell ", densityField[0])

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_1DRiemannProblem_Implicit_LineSearch()
