#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def SinglePhase_2DWallHeatedChannel_ChangeSect():

	#import themesh
	print("Reading a mesh with sudden cross-section change for test SinglePhase_2DWallHeatedChannel_ChangeSect()")
	M=cm.Mesh("../resources/VaryingSectionDuct.med")

    # Prepare the mesh boundaries
	xinf=0.0;
	xsup=0.01;
	yinf=0.0;
	ysup=0.01;
	eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Wall");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	M.setGroupAtPlan(yinf,1,eps,"Inlet");
	M.setGroupAtPlan(ysup,1,eps,"Outlet");

	#Nombre de cellules utilisees au depart dans Salome ou Alamos	
	nx=60
	ny=60
	#taille d'une cellule
	dx = (xsup-xinf)/nx
	dy = (ysup-yinf)/ny;
	for i in range(ny/2):
		 M.setGroupAtFaceByCoords((xsup-xinf)/4,(ysup-yinf)/4+(i+0.5)*dy,0,eps,"Wall");#Paroi verticale intérieure gauche
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4,(ysup-yinf)/4+(i+0.5)*dy,0,eps,"Wall");#Paroi verticale intérieure droitee
	
	for i in range(nx/4):
		 M.setGroupAtFaceByCoords((i+0.5)*dx,(ysup-yinf)/4,0,eps,"Wall");#paroi horizontale en bas à gauche
		 M.setGroupAtFaceByCoords((i+0.5)*dx,(ysup-yinf)*3/4,0,eps,"Wall");#paroi horizontale en haut à gauche
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4+(i+0.5)*dx,(ysup-yinf)/4,0,eps,"Wall");#paroi horizontale en bas à droite
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4+(i+0.5)*dx,(ysup-yinf)*3/4,0,eps,"Wall");#paroi horizontale en haut à droite

	spaceDim = M.getSpaceDimension();

    # set the limit field for each boundary
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=623;
	inletVelocityX=0;
	inletVelocityY=2.5;
	inletTemperature=563;
	outletPressure=155e5;

    # physical constants
	viscosite=[8.85e-5]
	conductivite=[1000]

	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar =myProblem.getNumberOfVariables();

    # Prepare for the initial condition
	VV_Constant =cm.Vector(nVar)

	# constant vector
	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletVelocityY;
	VV_Constant[3] = inletTemperature ;

    #Initial field creation
	print("Building mesh and initial data" );
	myProblem.setInitialFieldConstant(M,VV_Constant);

    # the boundary conditions
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure,[xsup,ysup]);
	myProblem.setInletBoundaryCondition("Inlet", inletTemperature, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("Wall", wallTemperature, wallVelocityX, wallVelocityY);

    # set physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);


	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Implicit);
	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
    
	# name file save
	fileName = "2DWallHeatedChannel_ChangeSect";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.5;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(float('inf'),20);#newton precision should be infinite for staggered scheme!!!
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
    SinglePhase_2DWallHeatedChannel_ChangeSect()
