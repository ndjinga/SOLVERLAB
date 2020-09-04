#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def FiveEqsTwoFluid_1DBoilingAssembly():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;
	M=cm.Mesh(xinf,xsup,nx)
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"outlet")
	M.setGroupAtPlan(xinf,0,eps,"inlet")

    # set the limit field for each boundary

	inletVoidFraction=0;
	inletVelocityX=[1]*2;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	heatPowerField=cm.Field("heatPowerField", cm.CELLS, M, 1);

	nbCells=M.getNumberOfCells();
	
	for i in range (nbCells):
		x=M.getCell(i).x();
		if (x> (xsup-xinf)/4) and (x< (xsup-xinf)*3/4):
			heatPowerField[i]=1e8
		else:
			heatPowerField[i]=0
	heatPowerField.writeVTK("heatPowerField",True)		
		

	myProblem = cf.FiveEqsTwoFluid(cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = inletVoidFraction;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX[0];
	VV_Constant[3] = inletVelocityX[1];
	VV_Constant[4] = inletTemperature ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant(M,VV_Constant)

    # set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet",inletVoidFraction,inletTemperature,inletVelocityX)
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);

    # set physical parameters
	myProblem.setHeatPowerField(heatPowerField)

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setEntropicCorrection(True);
	myProblem.setWellBalancedCorrection(True);  
    
    # name of result file
	fileName = "1DBoilingAssembly";

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
    FiveEqsTwoFluid_1DBoilingAssembly()
