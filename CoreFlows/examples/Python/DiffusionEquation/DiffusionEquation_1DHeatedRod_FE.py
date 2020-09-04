#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def DiffusionEquation_1DHeatedRod():
	spaceDim = 1;
    # Prepare for the mesh
	xinf = 0 ;
	xsup=4.2;
	nx=10;
	
	print ("Building of a 1D mesh with ", nx ," cells")

    # set the limit field for each boundary
	Temperature=623;
	cp_ur=300
	rho_ur=10000
	lambda_ur=5

	FECalculation=True    
	myProblem = cf.DiffusionEquation(spaceDim,FECalculation,rho_ur,cp_ur,lambda_ur);
	nVar = myProblem.getNumberOfVariables();

     #Set heat exchanges
	fluidTemp=573.;#fluid mean temperature
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient
	phi=1e5;#heat power ddensity
	myProblem.setFluidTemperature(fluidTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	myProblem.setHeatSource(phi);

	# constant vector
	VV_Constant = [623];
	
    #Initial field creation
	print("Building initial data" );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"Neumann","Neumann");

    # the boundary conditions
	myProblem.setNeumannBoundaryCondition("Neumann");

    # set the numerical method
	myProblem.setTimeScheme( cf.Explicit);
	# myProblem.setLinearSolver(GMRES,ILU,True);

    # name of result file
	fileName = "1DRodTemperature_FE";

    # computation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 100000000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

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
    DiffusionEquation_1DHeatedRod()
