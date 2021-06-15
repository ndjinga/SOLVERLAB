#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm
import VTK_routines
import exact_rs_stiffenedgas
import matplotlib.pyplot as plt
from numpy.linalg import norm
import time

test_desc={}
test_desc["Initial_data"]="Riemann problem"
test_desc["PDE_model"]="Euler"
test_desc["PDE_is_stationary"]=False
test_desc["PDE_search_for_stationary_solution"]=False
test_desc["Mesh_is_unstructured"]=False
test_desc["Part_of_mesh_convergence_analysis"]=True
test_desc["Numerical_method_name"]="Upwind"
test_desc["Numerical_method_space_discretization"]="Finite volumes"
test_desc["Boundary_conditions"]="Neumann"
test_desc["Geometry"]="Segment"


def solve(xinf,xsup,nx,cfl,isExplicit,scheme):
	start = time.time()

	if(isExplicit):
		ExplicitOrImplicit="Explicit"
	else:
		ExplicitOrImplicit="Implicit"

	test_desc["Numerical_method_time_discretization"]=ExplicitOrImplicit

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	discontinuity=(xinf+xsup)/2
	M=cm.Mesh(xinf,xsup,nx)
	eps=1e-6
	M.setGroupAtPlan(xsup,0,eps,"RightBoundary")
	M.setGroupAtPlan(xinf,0,eps,"LeftBoundary")

	test_desc["Space_dimension"]=M.getSpaceDimension()
	test_desc["Mesh_dimension"]=M.getMeshDimension()
	test_desc["Mesh_number_of_elements"]=M.getNumberOfCells()
	test_desc["Mesh_cell_type"]=M.getElementTypesNames()
		
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
	if(isExplicit):
		cf_ExplicitOrImplicit=cf.Implicit
	else:
		cf_ExplicitOrImplicit=cf.Explicit
			
	if(scheme=="Upwind"):
		cf_Scheme=cf.upwind
	elif(scheme=="Centered"):
		cf_Scheme=cf.centered

	myProblem.setNumericalScheme(cf_Scheme, cf_ExplicitOrImplicit);
    
    # name of result file
	fileName = "1DRiemannProblem_"+ExplicitOrImplicit+scheme;

    # simulation parameters 
	MaxNbOfTimeStep = nx ;
	freqSave = 10;
	maxTime = (xsup-xinf)/4./max(myProblem.getFluidEOS().vitesseSonPressure(initialPressure_Left,initialTemperature_Left), myProblem.getFluidEOS().vitesseSonPressure(initialPressure_Right,initialTemperature_Right) )
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveConservativeField(True);
	
	myProblem.setLinearSolver(cf.GMRES, cf.NOPC)
	myProblem.setNewtonSolver(precision,20, cf.Newton_SOLVERLAB)
	myProblem.usePrimitiveVarsInNewton(False)
 
    # evolution
	myProblem.initialize();

	ok = myProblem.run();

	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );
		####################### Postprocessing #########################

		dx=(xsup-xinf)/nx
		x=[ i*dx for i in range(nx+1)]   # array of cell center (1D mesh)
		fig, ([axDensity, axPressure], [axVelocity, axTemperature]) = plt.subplots(2, 2,sharex=True, figsize=(10,10))
		plt.gcf().subplots_adjust(wspace = 0.5)
	
		myEOS = myProblem.getFluidEOS()## Needed to retrieve gamma, pinfnity, convert (p,T) to density and (p, rho) to temperature
		initialDensity_Left = myEOS.getDensity( initialPressure_Left, initialTemperature_Left)
		initialDensity_Right = myEOS.getDensity( initialPressure_Right, initialTemperature_Right)

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
		error_velocity = norm( exactVelocity - velocityArray )/norm( exactVelocity )
		error_temperature = norm( exactTemperature - temperatureArray )/norm( exactTemperature )

		print("Absolute error = ", error_pressure, " (pressure), ", error_velocity, " (velocity), ", error_temperature, " (temperature), " )
			
		assert error_pressure <1.
		assert error_velocity <1.
		assert error_temperature <1.

	myProblem.terminate();

	end = time.time()

	test_desc["Computational_time_taken_by_run"]=end-start

	return pressureArray, velocityArray, temperatureArray, error_pressure, error_velocity, error_temperature, x, end - start 

if __name__ == """__main__""":
    solve(0.99,True,Upwind)
