#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_2DPoiseuilleFlow():
	spaceDim = 2;

	print("Building the mesh" );
    # Prepare the mesh data
	xinf = 0 ;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 4;
	nx   = 10;
	ny   = 40; 

	my_mesh=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	# set the boundary names for each boundary
	eps=1e-6;
	my_mesh.setGroupAtPlan(xsup,0,eps,"wall")
	my_mesh.setGroupAtPlan(xinf,0,eps,"wall")
	my_mesh.setGroupAtPlan(ysup,1,eps,"neumann")
	my_mesh.setGroupAtPlan(yinf,1,eps,"neumann")

    # physical constants
	viscosity=0.025
	viscosite=[viscosity];
	Vy_max             = 1.5
	a = -8*viscosity*Vy_max / ((ysup-yinf)*(ysup-yinf) ) #pressure slope
       
	initialTemperature = 573
	outletPressure     = 155e5

	initial_field=cm.Field("Initial field", cm.CELLS, my_mesh, 4)
	for i in range( 0 , my_mesh.getNumberOfCells() ):
		Ci=my_mesh.getCell(i)
		x=Ci.x()
		y=Ci.y()
		initial_field[i,0] =  outletPressure + a*(y - ysup )
		initial_field[i,1] =  0  #x component of the velocity
		initial_field[i,2] =  a/(2*viscosity)*( (x-(xsup+xinf)/2)*(x-(xsup+xinf)/2) - (xsup-xinf)*(xsup-xinf)/4)  #y component of the velocity
		initial_field[i,3] =  initialTemperature  
        
        
    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=573;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    #Initial field creation
	print("Setting initial data" );
	myProblem.setInitialField(initial_field)

    # the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setNeumannBoundaryCondition("neumann");

    # set physical parameters
	myProblem.setViscosity(viscosite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
	#myProblem.setLinearSolver(cf.GMRES,cf.LU,True);
    
	# name file save
	fileName = "2DPoiseuilleFlow";

	# parameters calculation
	MaxNbOfTimeStep = 10000 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setNewtonSolver(float('inf'),20);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
    SinglePhase_2DPoiseuilleFlow()
