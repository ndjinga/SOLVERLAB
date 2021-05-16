#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm
import matplotlib.pyplot as plt
from numpy.linalg import norm
import VTK_routines
import exact_rs_stiffenedgas

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

    ### Postprocessing initial data
	dx=(xsup-xinf)/nx
	x=[ i*dx for i in range(nx+1)]   # array of cell center (1D mesh)
	fig, ([axDensity, axPressure], [axVelocity, axTemperature]) = plt.subplots(2, 2,sharex=True, figsize=(10,10))
	plt.gcf().subplots_adjust(wspace = 0.5)

	myEOS = myProblem.getFluidEOS()## Needed to retrieve gamma, pinfnity, convert (p,T) to density and (p, rho) to temperature

	axPressure.set(xlabel='x (m)', ylabel='Pressure (bar)')
	axPressure.set_xlim(xinf,xsup)
	axPressure.set_ylim( initialPressure_Right - 0.1*(initialPressure_Left-initialPressure_Right), initialPressure_Left +  0.5*(initialPressure_Left-initialPressure_Right) )

	myPressureField = myProblem.getPressureField()
	myPressureField.writeVTK("PressureField")
	pressureArray=VTK_routines.Extract_VTK_data_over_line_to_numpyArray("PressureField"+"_0.vtu", [xinf,0,0], [xsup,0,0],nx)
	axPressure.plot(x, pressureArray,  label='Pressure time step 0')

	initialDensity_Left = myEOS.getDensity( initialPressure_Left, initialTemperature_Left)
	initialDensity_Right = myEOS.getDensity( initialPressure_Right, initialTemperature_Right)

	axDensity.set(xlabel='x (m)', ylabel='Density (Kg/m3)')
	axDensity.set_xlim(xinf,xsup)
	axDensity.set_ylim( initialDensity_Right - 0.1*(initialDensity_Left-initialDensity_Right), initialDensity_Left +  0.5*(initialDensity_Left-initialDensity_Right) )

	myDensityField = myProblem.getDensityField()
	myDensityField.writeVTK("DensityField")
	densityArray=VTK_routines.Extract_VTK_data_over_line_to_numpyArray("DensityField"+"_0.vtu", [xinf,0,0], [xsup,0,0],nx)
	axDensity.plot(x, densityArray,  label='Density time step 0')

	axVelocity.set(xlabel='x (m)', ylabel='Velocity (m/s)')
	axVelocity.set_xlim(xinf,xsup)
	axVelocity.set_ylim( 0.9*initialVelocity_Right - 0.1*(initialVelocity_Left-initialVelocity_Right), 22*initialVelocity_Left +  0.5*(initialVelocity_Left-initialVelocity_Right) )

	myVelocityField = myProblem.getVelocityXField()
	myVelocityField.writeVTK("VelocityField")
	velocityArray=VTK_routines.Extract_VTK_data_over_line_to_numpyArray("VelocityField"+"_0.vtu", [xinf,0,0], [xsup,0,0],nx)
	axVelocity.plot(x, velocityArray,  label='Velocity time step 0')
	
	axTemperature.set(xlabel='x (m)', ylabel='Temperature (K)')
	axTemperature.set_xlim(xinf,xsup)
	axTemperature.set_ylim( 0.999*initialTemperature_Right - 0.1*(initialTemperature_Left-initialTemperature_Right), 1.001*initialTemperature_Left +  0.5*(initialTemperature_Left-initialTemperature_Right) )

	myTemperatureField = myProblem.getTemperatureField()
	myTemperatureField.writeVTK("TemperatureField")
	temperatureArray=VTK_routines.Extract_VTK_data_over_line_to_numpyArray("TemperatureField"+"_0.vtu", [xinf,0,0], [xsup,0,0],nx)
	axTemperature.plot(x, temperatureArray,  label='Temperature time step 0')
	
	### run simulation
	ok = myProblem.run();

	#Determine exact solution
	exactDensity, exactVelocity, exactPressure = exact_rs_stiffenedgas.exact_sol_Riemann_problem(xinf, xsup, myProblem.presentTime(), myEOS.constante("gamma"), myEOS.constante("p0"), [ initialDensity_Left, initialVelocity_Left, initialPressure_Left ], [ initialDensity_Right, initialVelocity_Right, initialPressure_Right ], (xinf+xsup)/2, nx+1)

	### Plot curves
	axPressure.plot(x, exactPressure,  label='Exact Pressure ')
	myPressureField = myProblem.getPressureField()
	myPressureField.writeVTK("PressureField")
	pressureArray=VTK_routines. Extract_VTK_data_over_line_to_numpyArray("PressureField_"+str(myProblem.getNbTimeStep())+".vtu", [xinf,0,0], [xsup,0,0],nx)
	axPressure.plot(x, pressureArray,  label='Pressure time step '+str(myProblem.getNbTimeStep()))
	axPressure.legend()

	axDensity.plot(x, exactDensity,  label='Exact Density ')
	myDensityField = myProblem.getDensityField()
	myDensityField.writeVTK("DensityField")
	densityArray=VTK_routines. Extract_VTK_data_over_line_to_numpyArray("DensityField_"+str(myProblem.getNbTimeStep())+".vtu", [xinf,0,0], [xsup,0,0],nx)
	axDensity.plot(x, densityArray,  label='Density time step '+str(myProblem.getNbTimeStep()))
	axDensity.legend()

	axVelocity.plot(x, exactVelocity,  label='Exact Velocity ')
	myVelocityField = myProblem.getVelocityXField()
	myVelocityField.writeVTK("VelocityField")
	velocityArray=VTK_routines. Extract_VTK_data_over_line_to_numpyArray("VelocityField_"+str(myProblem.getNbTimeStep())+".vtu", [xinf,0,0], [xsup,0,0],nx)
	axVelocity.plot(x, velocityArray,  label='Velocity time step '+str(myProblem.getNbTimeStep()))
	axVelocity.legend()

	exactTemperature = [0.]*(nx+1)
	for i in range(nx+1):
		exactTemperature[i] = myEOS.getTemperatureFromPressure(exactPressure[i], exactDensity[i])

	axTemperature.plot(x, exactTemperature,  label='Exact Temperature ')
	myTemperatureField = myProblem.getTemperatureField()
	myTemperatureField.writeVTK("TemperatureField")
	temperatureArray=VTK_routines. Extract_VTK_data_over_line_to_numpyArray("TemperatureField_"+str(myProblem.getNbTimeStep())+".vtu", [xinf,0,0], [xsup,0,0],nx)
	axTemperature.plot(x, temperatureArray,  label='Temperature time step '+str(myProblem.getNbTimeStep()))
	axTemperature.legend()

	#plt.title('Solving Riemann problem for Euler equations\n with Finite volume schemes method')
	plt.savefig(fileName+".png")
	
	#Compute numerical error
	error_pressure = norm( exactPressure - pressureArray )/norm( exactPressure )
	print('relative error on pressure = ', error_pressure )
	
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );


	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_1DRiemannProblem_Implicit_LineSearch()
