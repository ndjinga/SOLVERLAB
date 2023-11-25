#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab as svl

def IsothermalSinglePhase_2DChannel(my_mesh, mesh_name, scheme, isCompressible ):
    spaceDim = 2;

    myProblem = svl.IsothermalSinglePhase(svl.Liquid,svl.around155bars600K,spaceDim,isCompressible);
    nVar =myProblem.getNumberOfVariables();

    # physical constants
    outletPressure     = 155e5
    inletVelocityX     = 1
    inletVelocityY     = 0

    print("Setting initial data" );
    initial_field=svl.Field("Initial field", svl.CELLS, my_mesh,3)
    for i in range( 0 , my_mesh.getNumberOfCells() ):
        Ci=my_mesh.getCell(i)
        x=Ci.x()
        y=Ci.y()
        initial_field[i,0] =  outletPressure
        initial_field[i,1] =  inletVelocityX  #x component of the velocity
        initial_field[i,2] =  inletVelocityY  #y component of the velocity
        
    myProblem.setInitialField(initial_field)

    # the boundary conditions
    myProblem.setWallBoundaryCondition("Wall");
    myProblem.setNeumannBoundaryCondition("Inlet");
    myProblem.setNeumannBoundaryCondition("Outlet");

    # set the numerical method
    myProblem.setNumericalScheme(getattr(svl, scheme), svl.Implicit )
    myProblem.setLinearSolver(svl.GMRES,svl.ILU);
    
    # name file save
    if( isCompressible ) :
        fileName = "2D"+mesh_name+"_Compressible_"+scheme+"_CFL1";
    else:
        fileName = "2D"+mesh_name+"_Incompressible_"+scheme+"_CFL1";
        
    # parameters calculation
    MaxNbOfTimeStep = 3 ;
    freqSave = 1;
    cfl = 1;
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
    
    msg = "Available fields : "
    for i in range(numberOfFields):
        msg += message[i]+" "
    print(msg)
    
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
    if len(sys.argv) >4:
        # Load the mesh
        print( "Loading mesh file " +  sys.argv[1]);
        my_mesh = svl.Mesh(sys.argv[1])
        mesh_name = sys.argv[2]
        scheme = sys.argv[3]
        isCompressible = bool(sys.argv[4])
        IsothermalSinglePhase_2DChannel(my_mesh, mesh_name, scheme, isCompressible)
    else :   
        raise ValueError("Provide mesh file name + mesh name + scheme + is compressible !!!")

