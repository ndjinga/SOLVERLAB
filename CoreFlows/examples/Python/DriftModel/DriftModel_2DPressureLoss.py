#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf
import cdmath as cm
import math 

def DriftModel_2DPressureLoss():
	spaceDim = 2;

    # Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=0.42;
	yinf=0.0;
	ysup=4.2;
	nx=2;
	ny=100; 
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
#	while i<= ndis:
	M.setGroupAtFaceByCoords(xinf+(i+0.5)*dx,yinf,0,eps,"inlet_1")
	i=1
#	while i<= nx:
	M.setGroupAtFaceByCoords(xinf+(i+0.5)*dx,yinf,0,eps,"inlet_2")
	
	
	
    # set the limit field for each boundary
	inletConc=0;    
	inletVelocityX1=0. ;
	inletVelocityY1=2.7 ;
	inletVelocityX2=0 ;
	inletVelocityY2=1.5 ;
	inletTemperature=563 ;
	outletPressure=155e5 ;


    # physical parameters
	pressureLossField=cm.Field("PressureLossCoeff", cm.FACES, M, 1);
	nbFaces=M.getNumberOfFaces();
	
	for i in range (nbFaces):
		if abs(M.getFace(i).y() - 1.0500)<eps :
			pressureLossField[i]=50 
			print("Premiere perte de charge ok")
		elif abs(M.getFace(i).y() - 1.680)<eps :
			pressureLossField[i]=50 
			print("Deuxieme perte de charge ok")
		elif M.getFace(i).y() == 2.7300:
			pressureLossField[i]=50 
			print("Troisieme perte de charge ok")
		else:
			pressureLossField[i]=0;
	
	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

        # Prepare for the initial condition
	VV_Constant_Left =cm.Vector(nVar)
	VV_Constant_Right =cm.Vector(nVar)
	
	# constant vectorsetOutletBoundaryCondition		
	VV_Constant_Left[0] = inletConc;
	VV_Constant_Left[1] = outletPressure ;
	VV_Constant_Left[2] = inletVelocityX1;
	VV_Constant_Left[3] = inletVelocityY1;
	VV_Constant_Left[4] = inletTemperature ;

	VV_Constant_Right[0] = inletConc;
	VV_Constant_Right[1] = outletPressure ;
	VV_Constant_Right[2] = inletVelocityX2;
	VV_Constant_Right[3] = inletVelocityY2;
	VV_Constant_Right[4] = inletTemperature ;

	
    # set physical parameters
	myProblem.setPressureLossField(pressureLossField);
	
	print("Building initial data" );
	myProblem.setInitialFieldStepFunction(M,VV_Constant_Left,VV_Constant_Right,discontinuity)

    
	# the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);
	myProblem.setInletBoundaryCondition("inlet_1", inletTemperature,inletConc, inletVelocityX1, inletVelocityY1);
	myProblem.setInletBoundaryCondition("inlet_2", inletTemperature,inletConc, inletVelocityX2, inletVelocityY2);
	myProblem.setWallBoundaryCondition("LeftWall", inletTemperature, 0,0);
	myProblem.setWallBoundaryCondition("RightWall", inletTemperature, 0,0);


	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setWellBalancedCorrection(True);  
	myProblem.setNonLinearFormulation(cf.VFFC) 

	# name file save
	fileName = "2DPressureLoss";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5000;
	precision = 1e-4;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(1e10,20);
	myProblem.saveVelocity();

	# evolution
	myProblem.initialize();
	print("Running python "+ fileName );

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
    DriftModel_2DPressureLoss()
