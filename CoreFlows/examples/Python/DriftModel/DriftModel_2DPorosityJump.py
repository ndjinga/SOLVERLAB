#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf
import cdmath as cm
import math 

def DriftModel_2DPorosityJump():
    spaceDim = 2;

    # Prepare for the mesh
    print("Building mesh " );
    xinf = 0 ;
    xsup=1.5;
    yinf=0.0;
    ysup=0.8;
    nx=20;
    ny=20; 
    discontinuity=(xinf+xsup)/2
    M=cm.Mesh(xinf,xsup,nx,yinf,ysup,ny)
    # set the limit field for each boundary
    eps=1e-6;
    M.setGroupAtPlan(xsup,0,eps,"Outlet")
    M.setGroupAtPlan(xinf,0,eps,"Inlet")
    M.setGroupAtPlan(ysup,1,eps,"Wall")
    M.setGroupAtPlan(yinf,1,eps,"Wall")
    dx=(xsup-xinf)/nx
    ndis=(discontinuity-xinf)/dx    
    
    
    # set the limit field for each boundary
    inletConc=0;    
    inletVelocityX=1. ;
    inletVelocityY=0 ;
    inletTemperature=563 ;
    outletPressure=155e5 ;


    # physical parameters
    porosityField=cm.Field("Porosity",cm.CELLS,M,1);
    for i in xrange(M.getNumberOfCells()):
        x=M.getCell(i).x();
        y=M.getCell(i).y();
        if x > (xsup-xinf)/3. and x<(xsup-xinf)*2./3 and  y> (ysup-yinf)/4. and y<(ysup-yinf)*3./4:
            porosityField[i]=0.5;
        else:
            porosityField[i]=1;
        pass
    porosityField.writeVTK("PorosityField");
    
    myProblem = cf.DriftModel(cf.around155bars600K,spaceDim);
    nVar =myProblem.getNumberOfVariables();

        # Prepare for the initial condition
    VV_Constant =cm.Vector(nVar)
    
    # constant vectorsetOutletBoundaryCondition        
    VV_Constant[0] = inletConc;
    VV_Constant[1] = outletPressure ;
    VV_Constant[2] = inletVelocityX;
    VV_Constant[3] = inletVelocityY;
    VV_Constant[4] = inletTemperature ;
    
    # set physical parameters
    myProblem.setPorosityField(porosityField);
    
    print("Building initial data");
    myProblem.setInitialFieldConstant(M,VV_Constant)

    
    # the boundary conditions
    myProblem.setOutletBoundaryCondition("Outlet", outletPressure);
    myProblem.setInletBoundaryCondition("Inlet", inletTemperature,inletConc, inletVelocityX, inletVelocityY);
    myProblem.setWallBoundaryCondition("Wall", inletTemperature, 0,0);
    myProblem.setNonLinearFormulation(cf.VFFC) 


    # set the numerical method
    myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
    myProblem.setWellBalancedCorrection(True);  

    # name file save
    fileName = "2DPorosityJump";

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
    DriftModel_2DPorosityJump()
