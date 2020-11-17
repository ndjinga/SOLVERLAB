#!/usr/bin/env python3
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath


def SinglePhase_2DSphericalExplosion_HEXAGON():

	inputfile="./resources/meshHexagonWithTriangles10.med";
	my_mesh=cdmath.Mesh(inputfile);
	spaceDim=2
	
	# Initial field data
	nVar=2+spaceDim;
	radius=0.5;
	Center=cdmath.Vector(spaceDim);#default value is (0,0,0)
	Vout=cdmath.Vector(nVar)
	Vin =cdmath.Vector(nVar)
	Vin[0]=155e5;
	Vin[1]=0;
	Vin[2]=0;
	Vin[3]=563;
	Vout[0]=154e5;
	Vout[1]=0;
	Vout[2]=0;
	Vout[3]=563;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);

	# Initial field creation
	print ("Setting mesh and initial data" ) ;
	myProblem.setInitialFieldSphericalStepFunction( my_mesh, Vout, Vin, radius, Center);

	# set the boundary conditions
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=563;

	myProblem.setWallBoundaryCondition("boundaries", wallTemperature, wallVelocityX, wallVelocityY);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    
	# name file save
	fileName = "2DSphericalExplosion_HEXAGON";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 5;
	cfl = 0.49;
	maxTime = 5;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveConservativeField(False);
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
    SinglePhase_2DSphericalExplosion_HEXAGON()
