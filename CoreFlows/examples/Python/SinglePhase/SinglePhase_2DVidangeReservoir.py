#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_2DVidangeReservoir():

	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=1.0;
	yinf=0.0;
	ysup=1.0;
	nx=50;
	ny=50; 
	diametreSortie=(ysup-yinf)/10.#10 percent of the height
	nsortie=ny*diametreSortie/(ysup-yinf)
	M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xinf,0,eps,"wall")
	M.setGroupAtPlan(ysup,1,eps,"inlet")
	M.setGroupAtPlan(yinf,1,eps,"wall")
	dy=(ysup-yinf)/ny
	i=0
        while i < nsortie:	
		M.setGroupAtFaceByCoords(xsup,yinf+(i+0.5)*dy,0,eps,"outlet")
                i=i+1
        while i < ny:	
		M.setGroupAtFaceByCoords(xsup,yinf+(i+0.5)*dy,0,eps,"wall")
                i=i+1
	
	

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=300;
	inletTemperature=300;
	outletPressure=1e5;

    # set the limit field for each boundary
	initialVelocityX=0;
	initialVelocityY=0;
	initialTemperature=300;
	initialPressure=1e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[0]=0;
	gravite[1]=-10;

	myProblem = cf.SinglePhase(cf.Liquid,cf.around1bar300K,spaceDim);
	nVar =myProblem.getNumberOfVariables();
    # Prepare for the initial condition
	VV_constant =cm.Vector(nVar)

	# constant vector
	VV_constant[0] = initialPressure ;
	VV_constant[1] = initialVelocityX;
	VV_constant[2] = initialVelocityY;
	VV_constant[3] = initialTemperature ;


    #Initial field creation
	print("Building initial data" );
        myProblem.setInitialFieldConstant( M, VV_constant)

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletPressureBoundaryCondition("inlet", outletPressure, inletTemperature,[xsup,yinf]);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);

	# name file save
	fileName = "2DVidangeReservoir";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = .05;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveVelocity();

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
    SinglePhase_2DVidangeReservoir()
