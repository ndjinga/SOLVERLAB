#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm
import math 

def main():

	spaceDim = 1;
    # Prepare for the mesh
	xinf = 0 ;
	xsup=4.2;
	eps=1.E-8
	M=cm.Mesh("../CoreFlows_src/examples/ressources/BifurcatingFlow2BranchesEqualSections.med")
	M.getFace(0).setGroupName("Inlet")
	M.getFace(31).setGroupName("Outlet")
	#M.setGroupAtPlan(xsup,0,eps,"Outlet");
	#M.setGroupAtPlan(xinf,0,eps,"Inlet");

    # set the initial field

	initialPressure=155e5;
	initialVelocityX=5;
	initialTemperature=573;

   # set the limit field for each boundary

	inletVelocityX=5;
	inletTemperature=573;
	outletPressure=155e5

 	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = initialPressure ;
	VV_Constant[1] = initialVelocityX;
	VV_Constant[2] = initialTemperature ;


    #Initial field creation
	print("Building initial data" ); 
	myProblem.setInitialFieldConstant( M, VV_Constant);

    # set the boundary conditions
	myProblem.setInletBoundaryCondition("Inlet", inletTemperature, inletVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet",outletPressure);

	#set porosity, heat and gravity source
	Sections=cm.Field("../CoreFlows_src/examples/ressources/BifurcatingFlow2BranchesEqualSections", cm.CELLS,"Section area");
	heatPowerField=cm.Field("../CoreFlows_src/examples/ressources/BifurcatingFlow2BranchesEqualSections", cm.CELLS,"Heat power");
	myProblem.setSectionField(Sections);
	myProblem.setHeatPowerField(heatPowerField)
	gravite=[-10]
	myProblem.setGravity(gravite)
    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setWellBalancedCorrection(True)    

    # name of result file
	fileName = "2BranchesHeatedChannels";

    # simulation parameters 
	MaxNbOfTimeStep = 100000 ;
	freqSave = 50;
	cfl = 0.5;
	maxTime = 500;
	precision = 1e-5;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	#myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass
 
    # evolution
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.run();
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
    main()
