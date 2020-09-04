#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def DriftModel_2DInclinedBoilingChannel():
	spaceDim = 2;

	# Prepare for the mesh
	xinf = 0 ;
	xsup=.6;
	yinf=0.0;
	ysup=2.0;
	nx=3;
	ny=100; 
	M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)

	#Set the barriers
	xcloison1=xinf+(xsup-xinf)/3
	xcloison2=xinf+2*(xsup-xinf)/3
	barrierField=cm.Field("Barrier Field", cm.FACES, M, 1);
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"wall")
	M.setGroupAtPlan(xinf,0,eps,"wall")
	M.setGroupAtPlan(ysup,1,eps,"outlet")
	M.setGroupAtPlan(yinf,1,eps,"inlet")
	dy=(ysup-yinf)/ny
	ncloison=3*ny/4
	i=0	
	while i<= ncloison+1:
		M.setGroupAtFaceByCoords(xcloison1,yinf+((ysup-yinf)/4)+(i+0.5)*dy,0,eps,"wall")
		M.setGroupAtFaceByCoords(xcloison2,yinf+((ysup-yinf)/4)+(i+0.5)*dy,0,eps,"wall")
		i=i+1

	nbFaces=M.getNumberOfFaces();
	for i in range (nbFaces):
		x=M.getFace(i).x();
		y=M.getFace(i).y();
		if ((y> yinf+(ysup-yinf)/4) and (abs(x-xcloison1)< eps or abs(x-xcloison2)< eps)) or abs(x-xinf)< eps or abs(x-xsup)< eps :
			barrierField[i]=1
		else:
			barrierField[i]=0
	barrierField.writeVTK("barrierField",True)		

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=563;
	inletConcentration=0;
	inletVelocityX=0;
	inletVelocityY=1;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[1]=-8.5;
	gravite[0]=5;

	heatPower=0e8;

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = inletConcentration;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletVelocityY;
	VV_Constant[4] = inletTemperature ;

    #Initial field creation
	print("Building mesh and initial data " );
	myProblem.setInitialFieldConstant(M,VV_Constant)

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConcentration, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.staggered, cf.Implicit);
	myProblem.setWellBalancedCorrection(True);    
	myProblem.setNonLinearFormulation(cf.VFFC) 
    
	# name of result file
	fileName = "2DInclinedChannelVFFCStaggeredWB3x100CFL100";

	# simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1000;
	cfl = 100;
	maxTime = 500;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.usePrimitiveVarsInNewton(True)

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
    DriftModel_2DInclinedBoilingChannel()
