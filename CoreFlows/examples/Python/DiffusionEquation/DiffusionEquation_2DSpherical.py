#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab

def DiffusionEquation_2DSpherical(FECalculation):

    # Prepare for the mesh
	inputfile="../resources/BoxWithMeshWithTriangularCells";
	fieldName="Temperature";
	spaceDim=2
	
    # set the limit field for each boundary
	Temperature=623;
	cp_ur=300
	rho_ur=10000
	lambda_ur=5

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,rho_ur,cp_ur,lambda_ur);

     #Set heat exchanges
	fluidTemp=573.;#fluid mean temperature
	heatTransfertCoeff=1000.;#fluid/solid exchange coefficient
	phi=1e5;#heat power ddensity
	myProblem.setFluidTemperature(fluidTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	myProblem.setHeatSource(phi);

	print("type(solverlab.CELLS) = ", type(solverlab.CELLS) )
	M=solverlab.Mesh(inputfile+".med")
	numberOfComponents=1
	time_iteration=0
	#temperature_field_cells=solverlab.Field(fieldName, solverlab.CELLS, M, numberOfComponents, time)
	#myProblem.setInitialField(temperature_field_cells)
	myProblem.setInitialField(inputfile,fieldName,time_iteration, solverlab.CELLS)

    #Initial field load
	print("Loading unstructured mesh and initial data" )
	if( FECalculation):
		myProblem.setInitialField(inputfile,fieldName,time_iteration, solverlab.NODES)
	else:
		myProblem.setInitialField(inputfile,fieldName,time_iteration, solverlab.CELLS)

    # the boundary conditions :
	if( FECalculation):
		boundaryNodeGroupNames=myProblem.getMesh().getNameOfNodeGroups()
		print(len(boundaryNodeGroupNames), " Boundary Node Group detected : ", boundaryNodeGroupNames)
	else:
		boundaryFaceGroupNames=myProblem.getMesh().getNameOfFaceGroups()
		print(len(boundaryFaceGroupNames), " Boundary Face Group detected : ", boundaryFaceGroupNames)

	myProblem.setNeumannBoundaryCondition("GAUCHE");
	myProblem.setNeumannBoundaryCondition("DROITE");
	myProblem.setNeumannBoundaryCondition("HAUT");
	myProblem.setNeumannBoundaryCondition("BAS");

    # set the numerical method
	myProblem.setTimeScheme( solverlab.Explicit);
	# myProblem.setLinearSolver(GMRES,ILU,True);

    # name of result file
	if( FECalculation):
		fileName = "2DSpherical_FE";
	else:
		fileName = "2DSpherical_FV";

    # computation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 0.95;
	maxTime = 100000000;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

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
    if len(sys.argv) >1 :
        FECalculation = bool(int(sys.argv[1]))
        DiffusionEquation_2DSpherical(FECalculation)
    else :
        raise ValueError("DiffusionEquation_2DSpherical : missing one argument")
