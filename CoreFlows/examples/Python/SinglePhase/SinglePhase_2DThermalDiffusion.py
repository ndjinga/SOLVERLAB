#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf
import cdmath as cm
import math 

def SinglePhase_2DThermalDiffusion():
	spaceDim = 2;

    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=0.2;
	yinf=0.0;
	ysup=0.4;
	nx=10;
	ny=20; 
	discontinuity=(xinf+xsup)/2
	M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"RightWall")
	M.setGroupAtPlan(xinf,0,eps,"LeftWall")
	M.setGroupAtPlan(ysup,1,eps,"outlet")
	dx=(xsup-xinf)/nx
	ndis=(discontinuity-xinf)/dx
	print("ndis=",math.floor(ndis) );
	i=0	
	while i<= ndis:
		M.setGroupAtFaceByCoords(xinf+(i+0.5)*dx,yinf,0,eps,"inlet_1")
		i=i+1
	while i<= nx:
		M.setGroupAtFaceByCoords(xinf+(i+0.5)*dx,yinf,0,eps,"inlet_2")
		i=i+1
	

        # set the limit field for each boundary
	inletVelocityX=0;
	inletVelocityY=1;
	inletTemperature1=563;
	inletTemperature2=593;
	outletPressure=155e5;

	viscosite= [1.5];
	conductivite=[5000];
	
	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

        # Prepare for the initial condition
	VV_Constant_Left =cm.Vector(nVar)
	VV_Constant_Right =cm.Vector(nVar)
	
	# constant vectors		
	VV_Constant_Left[0] = outletPressure ;
	VV_Constant_Left[1] = inletVelocityX;
	VV_Constant_Left[2] = inletVelocityY;
	VV_Constant_Left[3] = inletTemperature1 ;

	VV_Constant_Right[0] = outletPressure ;
	VV_Constant_Right[1] = inletVelocityX;
	VV_Constant_Right[2] = inletVelocityY;
	VV_Constant_Right[3] = inletTemperature2 ;
	
	print("Building initial data" );
	myProblem.setInitialFieldStepFunction(M,VV_Constant_Left,VV_Constant_Right,discontinuity)

    
	# the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);
	myProblem.setInletBoundaryCondition("inlet_1", inletTemperature1, inletVelocityX, inletVelocityY);
	myProblem.setInletBoundaryCondition("inlet_2", inletTemperature2, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("LeftWall", inletTemperature1, 0,0);
	myProblem.setWallBoundaryCondition("RightWall", inletTemperature2, 0,0);


	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);

	# name file save
	fileName = "2DThermalDiffusion";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = .01;
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
    SinglePhase_2DThermalDiffusion()
