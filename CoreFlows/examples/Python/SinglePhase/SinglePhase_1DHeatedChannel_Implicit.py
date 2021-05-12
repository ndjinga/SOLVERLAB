#!/usr/bin/env python3
# -*-coding:utf-8 -*

import CoreFlows as cf
import matplotlib.pyplot as plt
import cdmath as cm
import VTK_routines

def SinglePhase_1DHeatedChannel_Implicit():

	spaceDim = 1;
    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=4.2;
	nx=50;

    # set the limit field for each boundary

	inletVelocityX=5;
	inletTemperature=565;
	outletPressure=155e5;

    # physical parameters
	heatPower=1e8;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar;

	# constant vector
	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletTemperature ;


    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"inlet","outlet");

    # set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet",inletTemperature,inletVelocityX)
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup]);

    # set physical parameters
	myProblem.setHeatSource(heatPower);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
    
    # name of result file
	fileName = "1DHeatedChannelUpwind_Implicit";

    # simulation parameters 
	MaxNbOfTimeStep = 1000 ;
	freqSave = 100;
	cfl = 100;
	maxTime = 500;
	precision = 1e-7;

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
 
	myProblem.setLinearSolver(cf.GMRES, cf.ILU)
	myProblem.setNewtonSolver(1e-3, 50, cf.Newton_PETSC_LINESEARCH)

    # evolution
	myProblem.initialize();

    #Postprocessing
	plt.xlabel('x')
	plt.ylabel('Pressure')
	plt.xlim(xinf,xsup)
	plt.ylim( 0.999*outletPressure, 1.001*outletPressure )
	plt.title('Solving Riemann problem for Euler equations\n with Finite volume schemes method')
	dx=(xsup-xinf)/nx
	x=[ i*dx for i in range(nx)]   # array of cell center (1D mesh)

	myPressureField = myProblem.getPressureField()
	pressureArray=myPressureField.getFieldValues()
	line_pressure, = plt.plot(x, pressureArray,  label='Pressure time step 0')
	plt.legend()
	plt.savefig(fileName+".png")

	ok = myProblem.run();

	myPressureField = myProblem.getPressureField()
	timeStep=myProblem.getNbTimeStep()#Final time step
	pressureArray=myPressureField.getFieldValues()
	line_pressure, = plt.plot(x, pressureArray,  label='Pressure time step '+str(timeStep))
	plt.legend()
	plt.savefig(fileName+".png")

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
    SinglePhase_1DHeatedChannel_Implicit()
