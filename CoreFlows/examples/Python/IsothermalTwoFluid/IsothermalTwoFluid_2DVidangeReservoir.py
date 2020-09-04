#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def IsothermalTwoFluid_2DVidangeReservoir():

	spaceDim = 2;
	print("Building mesh " );
	# Prepare for the mesh
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
		M.setGroupAtFaceByCoords(xsup,yinf+(i+0.5)*dy,0,eps,"wall")
                i=i+1
        while i < ny:	
		M.setGroupAtFaceByCoords(xsup,yinf+(i+0.5)*dy,0,eps,"wall")
                i=i+1
	
	

    # set the limit field for each boundary
	wallVelocityX=[0]*2;
	wallVelocityY=[0]*2;
	inletAlpha=1;
	outletPressure=1e5;

    # set the limit field for each boundary
	initialAlphaTop=1;
	initialAlphaBottom=0.;
	initialVelocityX=[0]*2;
	initialVelocityY=[0]*2;
	initialPressure=1e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[0]=0;
	gravite[1]=-10;

	myProblem = cf.IsothermalTwoFluid(cf.around1bar300K,spaceDim);
	nVar =myProblem.getNumberOfVariables();
    # Prepare for the initial condition
	VV_top =cm.Vector(nVar)
	VV_bottom =cm.Vector(nVar)

	# top and bottom vectors
	VV_top[0] = initialAlphaTop ;
	VV_top[1] = initialPressure ;
	VV_top[2] = initialVelocityX[0];
	VV_top[3] = initialVelocityX[1];
	VV_top[4] = initialVelocityY[0];
	VV_top[5] = initialVelocityY[1];

	VV_bottom[0] = initialAlphaBottom ;
	VV_bottom[1] = initialPressure ;
	VV_bottom[2] = initialVelocityX[0];
	VV_bottom[3] = initialVelocityX[1];
	VV_bottom[4] = initialVelocityY[0];
	VV_bottom[5] = initialVelocityY[1];

    #Initial field creation
	print("Building initial data" );
	myProblem.setInitialFieldStepFunction( M, VV_bottom, VV_top, .8, 1);

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);
	myProblem.setNeumannBoundaryCondition("inlet");
	myProblem.setWallBoundaryCondition("wall",wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setNonLinearFormulation(cf.VFFC) 
	myProblem.setEntropicCorrection(True);

	# name file save
	fileName = "2DVidangeReservoir";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = .1;
	maxTime = 500;
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
    IsothermalTwoFluid_2DVidangeReservoir()
