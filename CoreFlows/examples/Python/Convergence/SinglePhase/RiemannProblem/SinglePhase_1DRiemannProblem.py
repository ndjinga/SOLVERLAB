#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

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


def SinglePhase_1DRiemannProblem_Implicit(xinf,xsup,nx,cfl,isExplicit,scheme)):
	start = time.time()
    if(isExplicit):
        test_desc["Numerical_method_time_discretization"]="Explicit"
    else:
        test_desc["Numerical_method_time_discretization"]="Implicit"

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
			
	if(scheme="Upwind"):
		cf_Scheme=cf.upwind
	elif(scheme="Centered"):
		cf_Scheme=cf.centered

	myProblem.setNumericalScheme(cf_Scheme, cf_ExplicitOrImplicit);
    
    # name of result file
	fileName = "1DRiemannProblem_"+ExplicitOrImplicit+scheme;

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	maxTime = 500;
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
		my_ResultField = myProblem.getOutputPressureField()
		#The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi
		max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(spaceDim*pi*pi)
		max_sol_num=my_ResultField.max()
		min_sol_num=my_ResultField.min()
		erreur_abs=0
		if method =='FE':
			for i in range(my_mesh.getNumberOfNodes()) :
				if  erreur_abs < abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) :
					erreur_abs = abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i])
		else:
			for i in range(my_mesh.getNumberOfCells()) :
				if  erreur_abs < abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) :
					erreur_abs = abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) 				
		print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
		print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
			
		assert erreur_abs/max_abs_sol_exacte <1.

	myProblem.terminate();

	end = time.time()

	test_desc["Absolute_error"]=erreur_abs
	test_desc["Relative_error"]=erreur_abs/max_abs_sol_exacte
	test_desc["Computational_time_taken_by_run"]=end-start

	return nx, end - start 

if __name__ == """__main__""":
    SinglePhase_1DRiemannProblem_Implicit(0.99,True,Upwind)
