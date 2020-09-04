#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf
import cdmath as cm


def DriftModel_2DInclinedBoilingChannelBarrier():

	spaceDim = 2;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=2.0;
	yinf=0.0;
	ysup=4.0;
	nx=20;
	ny=40; 
	cloison=(xinf+xsup)/2
	M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"wall")
	M.setGroupAtPlan(xinf,0,eps,"wall")
	M.setGroupAtPlan(ysup,1,eps,"outlet")
	M.setGroupAtPlan(yinf,1,eps,"inlet")
	dy=(ysup-yinf)/ny
	ncloison=ny/2
	i=0	
	while i<= ncloison+1:
		M.setGroupAtFaceByCoords(cloison,((ysup-yinf)/4)+(i+0.5)*dy,0,eps,"wall")
		i=i+1
	

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=573;
	inletConc=0;
	inletVelocityX=0;
	inletVelocityY=1;
	inletTemperature=563;
	outletPressure=155e5;

    # source terms
	gravite = [0] * spaceDim
    
	gravite[0]=7;
	gravite[1]=-7;

	heatPower=1e8;
	heatPower_nul=0;
	heatPowerField=cm.Field("heatPowerField", cm.CELLS, M, 1);
	
	nbCells=M.getNumberOfCells();
	
	for i in range (nbCells):
		x=M.getCell(i).x();
		y=M.getCell(i).y();
		if (y> (ysup-yinf)/4) and (y< (ysup-yinf)*3/4) and (x<cloison):
			heatPowerField[i]=heatPower
		else:
			heatPowerField[i]=heatPower_nul
	heatPowerField.writeVTK("heatPowerField",True)		
		

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();
    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = inletConc ;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletVelocityY;
	VV_Constant[4] = inletTemperature ;

    #Initial field creation
	print("Building initial data" ); 
	myProblem.setInitialFieldConstant(M,VV_Constant)

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConc, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setHeatPowerField(heatPowerField)
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setWellBalancedCorrection(True);
	myProblem.setNonLinearFormulation(cf.VFFC) 

	# name file save
	fileName = "2DBInclinedoilingChannelBarrier";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = .5;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
    DriftModel_2DInclinedBoilingChannelBarrier()
