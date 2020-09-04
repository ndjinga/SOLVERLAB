#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf

def SinglePhase_2DLidDrivenCavity_unstructured():
	spaceDim = 2;

	print( "Loading unstructured mesh " );
	inputfile="../resources/BoxWithMeshWithTriangularCells.med";	
 
    # set the limit field for each boundary

	fixedWallVelocityX=0;
	fixedWallVelocityY=0;
	fixedWallTemperature=273;

	movingWallVelocityX=1;
	movingWallVelocityY=0;
	movingWallTemperature=273;
	
    # physical constants

	viscosite=[0.025];

	myProblem = cf.SinglePhase(cf.Gas,cf.around1bar300K,spaceDim);
	nVar = myProblem.getNumberOfVariables();

	#Initial field creation
	print("Building initial data " ); 
	
    # Prepare for the initial condition
    
	VV_Constant = [0] * nVar

	# constant vector
	VV_Constant[0] = 1e5;
	VV_Constant[1] = 0 ;
	VV_Constant[2] = 0;
	VV_Constant[3] = 273;

    #Initial field creation
	print("Setting mesh and initial data" );
	myProblem.setInitialFieldConstant(inputfile,VV_Constant);

    # Set the boundary conditions
	myProblem.setWallBoundaryCondition("BAS", fixedWallTemperature, fixedWallVelocityX, fixedWallVelocityY);
	myProblem.setWallBoundaryCondition("GAUCHE", fixedWallTemperature, fixedWallVelocityX, fixedWallVelocityY);
	myProblem.setWallBoundaryCondition("DROITE", fixedWallTemperature, fixedWallVelocityX, fixedWallVelocityY);
	myProblem.setWallBoundaryCondition("HAUT", movingWallTemperature, movingWallVelocityX, movingWallVelocityY);

    # set physical parameters
	myProblem.setViscosity(viscosite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.pressureCorrection, cf.Implicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
   
    # name file save
	fileName = "2DLidDrivenCavityUnstructured";

    # simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 10;
	maxTime = 50;
	precision = 1e-9;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass

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
    SinglePhase_2DLidDrivenCavity_unstructured()
