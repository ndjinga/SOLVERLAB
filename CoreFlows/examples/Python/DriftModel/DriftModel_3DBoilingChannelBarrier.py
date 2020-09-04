#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf
import cdmath as cm


def DriftModel_3DBoilingChannelBarrier():

	spaceDim = 3;
	# Prepare for the mesh
	print("Building mesh " );
	xinf = 0 ;
	xsup=2.0;
	yinf=0.0;
	ysup=2.0;
	zinf=0.0;
	zsup=4.0;
	nx=10;
	ny=nx; 
	nz=20; 
	xcloison=(xinf+xsup)/2
	ycloison=(yinf+ysup)/2
	zcloisonmin=1
	zcloisonmax=3
	M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny,zinf,zsup,nz)
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"wall")
	M.setGroupAtPlan(xinf,0,eps,"wall")
	M.setGroupAtPlan(ysup,1,eps,"wall")
	M.setGroupAtPlan(yinf,1,eps,"wall")
	M.setGroupAtPlan(zsup,2,eps,"outlet")
	M.setGroupAtPlan(zinf,2,eps,"inlet")
	dx=(xsup-xinf)/nx
	dy=(ysup-yinf)/ny
	dz=(zsup-zinf)/nz
	ncloison=nz*(zcloisonmax-zcloisonmin)/(zsup-zinf)
	i=0	
	j=0
	while i< ncloison:
		while j< ny:
			M.setGroupAtFaceByCoords(xcloison,(j+0.5)*dy,zcloisonmin+(i+0.5)*dz,eps,"wall")
			M.setGroupAtFaceByCoords((j+0.5)*dx,ycloison,zcloisonmin+(i+0.5)*dz,eps,"wall")
			j=j+1
		i=i+1

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallVelocityZ=0;
	wallTemperature=573;
	inletConc=0;
	inletVelocityX=0;
	inletVelocityY=0;
	inletVelocityZ=1;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	gravite = [0] * spaceDim
    
	gravite[0]=0;
	gravite[1]=0;
	gravite[2]=-10;

	heatPower1=0;
	heatPower2=0.25e8;
	heatPower3=0.5e8;
	heatPower4=1e8;

	myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();
	heatPowerField=cm.Field("heatPowerField", cm.CELLS, M, 1);
	
	nbCells=M.getNumberOfCells();
	
	for i in range (nbCells):
		x=M.getCell(i).x();
		y=M.getCell(i).y();
		z=M.getCell(i).z();
		if (z> zcloisonmin) and (z< zcloisonmax) :
			if (y<ycloison) and (x<xcloison):
				heatPowerField[i]=heatPower1
			if (y<ycloison) and (x>xcloison):
				heatPowerField[i]=heatPower2
			if (y>ycloison) and (x<xcloison):
				heatPowerField[i]=heatPower3
			if (y>ycloison) and (x>xcloison):
				heatPowerField[i]=heatPower4
		else:
			heatPowerField[i]=0
			
	heatPowerField.writeVTK("heatPowerField",True)		
		
    # Prepare for the initial condition
	VV_Constant =[0]*nVar

	# constant vector
	VV_Constant[0] = inletConc ;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletVelocityY;
	VV_Constant[4] = inletVelocityZ;
	VV_Constant[5] = inletTemperature ;

    #Initial field creation
	print("Building initial data " ); 
	myProblem.setInitialFieldConstant(M,VV_Constant)

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,[xsup,ysup,zsup]);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConc, inletVelocityX, inletVelocityY, inletVelocityZ);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY,wallVelocityZ);

    # set physical parameters
	myProblem.setHeatPowerField(heatPowerField)
	myProblem.setGravity(gravite);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
	myProblem.setWellBalancedCorrection(True);

	# name file save
	fileName = "3DBoilingChannelBarrier";

	# parameters calculation
	MaxNbOfTimeStep = 3;
	freqSave = 100;
	cfl = .3;
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
    DriftModel_3DBoilingChannelBarrier()
