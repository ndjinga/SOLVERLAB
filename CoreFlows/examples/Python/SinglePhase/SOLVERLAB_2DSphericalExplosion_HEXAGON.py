#!/usr/bin/env python3
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath
import sys

def SinglePhase_2DSphericalExplosion_HEXAGON(inputfile):

	my_mesh = cdmath.Mesh(inputfile)
	if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
		raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")
	spaceDim=2
	
	# Initial field data
	nVar=2+spaceDim;
	radius=0.5;
	Center=cdmath.Vector(spaceDim);#default value is (0,0,0)
	Vout = cdmath.Vector(nVar)
	Vin  = cdmath.Vector(nVar)
	
	Pmin=100e5
	Pmax=155e5
	InitialTemperature = 563

	Vin[0]=Pmax
	Vin[1]=0;
	Vin[2]=0;
	Vin[3]=InitialTemperature
	Vout[0]=Pmin;
	Vout[1]=0;
	Vout[2]=0;
	Vout[3]=InitialTemperature

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);

	# Initial field creation
	print ("Setting mesh and initial data" ) ;
	myProblem.setInitialFieldSphericalStepFunction( my_mesh, Vin, Vout, radius, Center);

	# set the boundary conditions
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=InitialTemperature;

	myProblem.setWallBoundaryCondition("BoundaryFaces", wallTemperature, wallVelocityX, wallVelocityY);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    
	# name file save
	fileName = "2DSphericalExplosion_HEXAGON";

	# parameters calculation
	MaxNbOfTimeStep = 400 ;
	freqSave = 25;
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
	myProblem.setSaveFileFormat(cf.MED)


	# evolution
	myProblem.initialize();

	masseInitiale = abs( myProblem.getDensityField().integral()[0] )
	initialTotalEnergy=abs(myProblem.getTotalEnergyField().integral()[0] )

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );

	# Control of the results
	PressureField=myProblem.getPressureField()
	print( "La pression est bornée par le maximum de pression initiale : PressureField.max() < InitialPressure.max()" )
	assert PressureField.max() < Pmax

	temperatureField=myProblem.getTemperatureField()
	print("La température est constante à 0.1 % près : erreur relative = ", max(abs(temperatureField.max() - InitialTemperature),abs(temperatureField.min() - InitialTemperature))/InitialTemperature )
	assert abs(temperatureField.max() - InitialTemperature)/InitialTemperature < 0.001
	assert abs(temperatureField.min() - InitialTemperature)/InitialTemperature < 0.001

	densityField=myProblem.getDensityField()
	print("La masse totale est conservée au zero machine près : erreur relative = ", abs( (abs(densityField.integral()[0]) - masseInitiale)/masseInitiale ) )
	assert abs( (abs(densityField.integral()[0]) - masseInitiale)/masseInitiale ) < 1e-14
	
	totalEnergyField=myProblem.getTotalEnergyField()
	print("L'énergie totale est conservée au zero machine près : erreur relative = ", abs( (abs(totalEnergyField.integral()[0]) - initialTotalEnergy)/initialTotalEnergy ) )
	assert abs( (abs(totalEnergyField.integral()[0]) - initialTotalEnergy)/initialTotalEnergy ) < 1e-13

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
	if len(sys.argv) > 1 :
		inputfile=sys.argv[1]
		SinglePhase_2DSphericalExplosion_HEXAGON(inputfile)
	else :
		raise ValueError("Error : Expecting a mesh file name argument")
