#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf
import cdmath as cm

def DriftModel_1DPorosityJump():

    # Prepare for the mesh
    print("Building mesh " );
    xinf = 0 ;
    xsup=4.2;
    nx=100;
    M=cm.Mesh(xinf,xsup,nx)
    spaceDim = M.getSpaceDimension()

    # set the limit field for each boundary
    initialConc=0;
    initialVelocityX=1;
    initialTemperature=600;
    initialPressure=155e5;

    # physical parameters
    porosityField=cm.Field("Porosity",cm.CELLS,M,1);
    for i in xrange(M.getNumberOfCells()):
        x=M.getCell(i).x();
        if x > (xsup-xinf)/3. and x<(xsup-xinf)*2./3:
            porosityField[i]=0.5;
        else:
            porosityField[i]=1;
        pass
    porosityField.writeVTK("PorosityField");

    myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
    nVar =  myProblem.getNumberOfVariables();

    # Prepare for the initial condition
    VV_Constant =[0]*nVar;

	# constant vector
    VV_Constant[0] = initialConc;
    VV_Constant[1] = initialPressure ;
    VV_Constant[2] = initialVelocityX;
    VV_Constant[3] = initialTemperature ;


    #Initial field creation
    print("Building initial data " ); 
    #myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Inlet","Outlet");
    myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Neumann","Neumann");

    # set the boundary conditions
    #myProblem.setInletBoundaryCondition("Inlet",initialTemperature,initialConc,initialVelocityX)
    #myProblem.setOutletBoundaryCondition("Outlet", initialPressure,[xsup]);
    myProblem.setNeumannBoundaryCondition("Neumann");

    # set physical parameters
    myProblem.setPorosityField(porosityField);

    # set the numerical method
    myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    myProblem.setWellBalancedCorrection(True);  
    myProblem.setNonLinearFormulation(cf.VFFC) 
    
    # name of result file
    fileName = "1DPorosityJumpUpwindWB";

    # simulation parameters 
    MaxNbOfTimeStep = 3 ;
    freqSave = 1;
    cfl = 0.95;
    maxTime = 500;
    precision = 1e-5;

    myProblem.setCFL(cfl);
    myProblem.setPrecision(precision);
    myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
    myProblem.setTimeMax(maxTime);
    myProblem.setFreqSave(freqSave);
    myProblem.setFileName(fileName);
    myProblem.setNewtonSolver(precision,20);
    myProblem.saveConservativeField(True);
 
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
    DriftModel_1DPorosityJump()
