#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath


def SinglePhase_3DSphericalExplosion_unstructured():

	inputfile="../../resources/meshCube.med";

	xinf=0;
	xsup=1;
	yinf=0;
	ysup=1;
	zinf=0;
	zsup=1;
	M=cdmath.Mesh(inputfile);
	eps=1.E-6;
	M.setGroupAtPlan(xinf,0,eps,"GAUCHE");
	M.setGroupAtPlan(xsup,0,eps,"DROITE");
	M.setGroupAtPlan(yinf,1,eps,"ARRIERE");
	M.setGroupAtPlan(ysup,1,eps,"AVANT");
	M.setGroupAtPlan(zinf,2,eps,"BAS");
	M.setGroupAtPlan(zsup,2,eps,"HAUT");

	# Initial field data
	spaceDim = 3;
	nVar=2+spaceDim;
	radius=0.5;
	Center=cdmath.Vector(3);#default value is (0,0,0)
	Vout=cdmath.Vector(nVar)
	Vin =cdmath.Vector(nVar)
	Vin[0]=1.1;
	Vin[1]=0;
	Vin[2]=0;
	Vin[3]=0;
	Vin[4]=300;
	Vout[0]=1;
	Vout[1]=0;
	Vout[2]=0;
	Vout[3]=0;
	Vout[4]=300;
	
	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar=myProblem.getNumberOfVariables();

	# Initial field creation
	print ("Setting mesh and initial data" ) ;
	myProblem.setInitialFieldSphericalStepFunction( M, Vout, Vin, radius, Center);

	# set the boundary conditions
	wallVelocityX=0;
	wallVelocityY=0;
	wallVelocityZ=0;
	wallTemperature=563;

	myProblem.setWallBoundaryCondition("GAUCHE", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("DROITE", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("HAUT", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("BAS", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("AVANT", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("ARRIERE", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
 	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
	myProblem.setEntropicCorrection(False);
	myProblem.setWellBalancedCorrection(False);
    
	# name file save
	fileName = "3DSphericalExplosion_unstructured";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5;
	precision = 1e-6;

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
    SinglePhase_3DSphericalExplosion_unstructured()
