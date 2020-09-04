#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm
 

def DriftModel_2BranchesBoilingChannels():

	spaceDim = 1;
    # Prepare for the mesh
	M=cm.Mesh("../resources/BifurcatingFlow2BranchesEqualSections.med")
	M.getFace(0).setGroupName("Inlet")#z=0
	M.getFace(31).setGroupName("Outlet")#z=4.2

    # set the initial field
	initialConc=0;
	initialPressure=155e5;
	initialVelocityX=5;
	initialTemperature=573;

   # set the limit field for each boundary
	inletConc=0;
	inletVelocityX=5;
	inletTemperature=573;
	outletPressure=155e5

 	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = initialConc ;
	VV_Constant[1] = initialPressure ;
	VV_Constant[2] = initialVelocityX;
	VV_Constant[3] = initialTemperature ;


    #Initial field creation
	print("Building initial data" ); 
	myProblem.setInitialFieldConstant( M, VV_Constant);

    # set the boundary conditions
	myProblem.setInletBoundaryCondition("Inlet", inletTemperature, inletConc,inletVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet",outletPressure);

	#set porosity, heat and gravity source
	Sections=cm.Field("../resources/BifurcatingFlow2BranchesEqualSections", cm.CELLS,"Section area");
	heatPowerField=cm.Field("../resources/BifurcatingFlow2BranchesEqualSections", cm.CELLS,"Heat power");
	myProblem.setSectionField(Sections);
	myProblem.setHeatPowerField(heatPowerField)
	gravite=[-10]
	myProblem.setGravity(gravite)
    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setWellBalancedCorrection(True)    
	myProblem.setNonLinearFormulation(cf.VFFC) 

    # name of result file
	fileName = "2BranchesBoilingChannels";

    # simulation parameters 
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveConservativeField(True);
 
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

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    DriftModel_2BranchesBoilingChannels()
