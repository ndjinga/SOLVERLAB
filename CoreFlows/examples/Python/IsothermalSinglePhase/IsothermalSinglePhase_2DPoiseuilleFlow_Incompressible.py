#!/usr/bin/env python3
# -*-coding:utf-8 -*

import solverlab as svl

def IsothermalSinglePhase_2DPoiseuilleFlow():
	spaceDim = 2;

	print("Building the mesh" );
    # Prepare the mesh data
	xinf = 0 ;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 4;
	nx   = 10;
	ny   = 40; 

	my_mesh=svl.Mesh(xinf,xsup,nx,yinf,ysup,ny)
	# set the boundary names for each boundary
	eps=1e-6;
	my_mesh.setGroupAtPlan(xsup,0,eps,"Wall")
	my_mesh.setGroupAtPlan(xinf,0,eps,"Wall")
	my_mesh.setGroupAtPlan(ysup,1,eps,"Neumann")
	my_mesh.setGroupAtPlan(yinf,1,eps,"Neumann")
        
	myProblem = svl.IsothermalSinglePhase(svl.Liquid,svl.around155bars600K,spaceDim,False);
	nVar =myProblem.getNumberOfVariables();

    # physical constants
	viscosity=0.0025
	viscosite=[viscosity];
	Vy_max             = 1.5
	a = -8*viscosity*Vy_max / ((ysup-yinf)*(ysup-yinf) ) #pressure slope
       
	outletPressure     = 155e5

	print("Setting initial data" );
	initial_field=svl.Field("Initial field", svl.CELLS, my_mesh,3)
	for i in range( 0 , my_mesh.getNumberOfCells() ):
		Ci=my_mesh.getCell(i)
		x=Ci.x()
		y=Ci.y()
		initial_field[i,0] =  outletPressure + a*(y - ysup )
		initial_field[i,1] =  0  #x component of the velocity
		initial_field[i,2] =  a/(2*viscosity)*( (x-(xsup+xinf)/2)*(x-(xsup+xinf)/2) - (xsup-xinf)*(xsup-xinf)/4)  #y component of the velocity
        
	myProblem.setInitialField(initial_field)

    # the boundary conditions
	myProblem.setWallBoundaryCondition("Wall");
	myProblem.setNeumannBoundaryCondition("Neumann");

    # set physical parameters
	myProblem.setViscosity(viscosite);

	# set the numerical method
	myProblem.setNumericalScheme(svl.staggered, svl.Implicit);
	myProblem.setLinearSolver(svl.GMRES,svl.LU);
    
	# name file save
	fileName = "2DPoiseuilleFlow_Incompressible_staggered_CFL10";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 10;
	maxTime = 5000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	#myProblem.setNewtonSolver(float('inf'),20);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
	
	message=myProblem.getOutputFieldsNames()
	numberOfFields=len(message)
	
	for i in range(numberOfFields):
		print( message[i])
	
	pressureField=myProblem.getOutputField("Pressure")
	pressureField.writeMED("pressureField")
	pressureField.writeVTK("pressureField")
	print("Pressure in first cell ", pressureField[0])
	
	PressureField=myProblem.getPressureField()
	PressureField.writeMED("PressureField2")
	PressureField.writeVTK("PressureField2")
	print("Pressure in first cell ", pressureField[0])

	velocityField=myProblem.getVelocityField()
	velocityField.writeMED("velocityField")
	velocityField.writeVTK("velocityField")
	print("Velocity X in first cell ", velocityField[0,0])
	print("Velocity Y in first cell ", velocityField[0,1])

	velocityField=myProblem.getOutputField("Velocity")
	velocityField.writeMED("velocityField2")
	velocityField.writeVTK("velocityField2")
	print("Velocity X in first cell ", velocityField[0,0])
	print("Velocity Y in first cell ", velocityField[0,1])

	densityField=myProblem.getDensityField()
	densityField.writeMED("densityField")
	densityField.writeVTK("densityField")
	print("Density in first cell ", densityField[0])

	densityField=myProblem.getOutputField("Density")
	densityField.writeMED("densityField2")
	densityField.writeVTK("densityField2")
	print("Density in first cell ", densityField[0])

	momentumField=myProblem.getMomentumField()
	momentumField.writeMED("momentumField")
	momentumField.writeVTK("momentumField")
	print("Momentum X in first cell ", momentumField[0,0])
	print("Momentum Y in first cell ", momentumField[0,1])

	momentumField=myProblem.getOutputField("Momentum")
	momentumField.writeMED("momentumField2")
	momentumField.writeVTK("momentumField2")
	print("Momentum X in first cell ", momentumField[0,0])
	print("Momentum Y in first cell ", momentumField[0,1])

	velocityXField=myProblem.getVelocityXField()
	velocityXField.writeMED("velocityXField")
	velocityXField.writeVTK("velocityXField")
	print("Velocity X in first cell ", velocityXField[0])

	velocityXField=myProblem.getOutputField("VelocityX")
	velocityXField.writeMED("velocityXField2")
	velocityXField.writeVTK("velocityXField2")
	print("Velocity X in first cell ", velocityXField[0])

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    IsothermalSinglePhase_2DPoiseuilleFlow()
