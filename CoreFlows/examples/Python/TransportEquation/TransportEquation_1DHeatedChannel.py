#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf

def TransportEquation_1DHeatedChannel():

	spaceDim = 1;
    # Prepare for the mesh
	xinf = 0 ;
	xsup=4.2;
	nx=10;

    # set the limit field for each boundary
	inletEnthalpy=1.3e6;

    # Set the transport velocity
	transportVelocity=[5];

	myProblem = cf.TransportEquation(cf.LiquidPhase,cf.around155bars600KTransport,transportVelocity);
	nVar = myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant = [1.3e6]; #initial enthalpy

	#Set rod temperature and heat exchamge coefficient
	rodTemp=623;#Rod clad temperature 
	heatTransfertCoeff=1000;#fluid/solid heat exchange coefficient
	myProblem.setRodTemperature(rodTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);

    #Initial field creation
	print("Building mesh and initial data " );
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"inlet","neumann");
 
    # Set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet", inletEnthalpy);
	myProblem.setNeumannBoundaryCondition("neumann")

    # Set the numerical method
	myProblem.setTimeScheme( cf.Explicit);

    # name file save
	fileName = "1DFluidEnthalpy";

    # parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 5;
	cfl = 0.95;
	maxTime = 5;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);


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
    TransportEquation_1DHeatedChannel()
