#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab
import medcoupling
import MEDLoader

#===============================================================================================================================
# Name        : Simulation of a 2D heat equation 
# Description : Test solving the diffusion of the temperature T in a solid
#               \rho cp dT/dt-\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T_{solid}) 
#               Neumann or Dirichlet boundary conditions passed through boundary fields (not yet active)
#               Finite elements or finite volumes
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
#================================================================================================================================


def DiffusionEquation_2DSpherical(FECalculation, fileName):

    """ Description : Test solving the diffusion of the temperature T in a solid (default is Uranium). 
        Equation : Thermal diffusion equation  \rho cp dT/dt-\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T_{solid})
                Heat capacity cp, density \rho, and conductivity \lambda of the solid MUST be defined
                The solid may be refrigerated by a fluid with temperature T_{solid} transfer coefficient \lambda_{sf} using functions setFluidTemperature and setHeatTransfertCoeff
                The solid may receive some extra heat power due to nuclear fissions using function setHeatSource.
    """
    #Space dimension of the problem
    spaceDim=2
    
    # Mandatory physical values
    solid_specific_heat=300# specific heat capacity, default value 300
    solid_density=10000# density, default value 10000
    solid_conductivity=5# conductivity, default value 5

    myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,solid_density,solid_specific_heat,solid_conductivity);

    # Definition of field support parameter
    if( FECalculation):
        supportOfField=solverlab.NODES
        typeOfField=medcoupling.ON_NODES
    else:
        supportOfField=solverlab.CELLS    
        typeOfField=medcoupling.ON_CELLS
    
    # Set the mesh and initial data
    use_field_in_memory=False
    if(not use_field_in_memory):#On lit le champ dans un fichier
        initial_data_inputfile="../resources/meshSquare";
        initial_data_fieldName="Solid temperature";
        print("Loading unstructured mesh and initial data", " in file ", initial_data_inputfile )
        initial_data_time_iteration=0# default value is 0
        initial_data_time_sub_iteration=0# default value is 0
        initial_data_time_meshLevel=0# default value is 0
        myProblem.setInitialField(initial_data_inputfile, initial_data_fieldName, initial_data_time_iteration, initial_data_time_sub_iteration, initial_data_time_meshLevel, supportOfField)
    else:#On utilise un champ en mémoire
        my_field_in_memory=MEDLoader.ReadField(typeOfField,"../resources/meshSquare.med","Mesh_1",0,"Solid temperature",0,0)
        myProblem.setInitialField( my_field_in_memory )
        
    #### Optional physical values (default value is zero) : fluid temperature field, heat transfert coefficient, heat power field 
    # Loading and setting fluid temperature field
    fluid_temperature_inputfile="../resources/meshSquare";
    fluid_temperature_fieldName="Fluid temperature";
    fluid_temperature_time_iteration=0# default value is 0
    fluid_temperature_time_sub_iteration=0# default value is 0
    fluid_temperature_meshLevel=0# default value is 0
    print("Loading field :", fluid_temperature_fieldName, " in file ", fluid_temperature_inputfile)
    fluidTemperatureField=solverlab.Field(fluid_temperature_inputfile, supportOfField, fluid_temperature_fieldName, fluid_temperature_time_iteration, fluid_temperature_time_sub_iteration, fluid_temperature_meshLevel)
    myProblem.setFluidTemperatureField(fluidTemperatureField)
    # Setting heat transfert coefficient
    heatTransfertCoeff=1000.;#fluid/solid exchange coefficient, default value is 0
    myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
    # Loading heat power field
    heat_power_inputfile="../resources/meshSquare";
    heat_power_fieldName="Heat power";
    heat_power_time_iteration=0# default value is 0
    heat_power_time_sub_iteration=0# default value is 0
    heat_power_meshLevel=0# default value is 0
    print("Loading field :", heat_power_fieldName, " in file ", heat_power_inputfile)
    heatPowerField=solverlab.Field(heat_power_inputfile, supportOfField, heat_power_fieldName, heat_power_time_iteration, heat_power_time_sub_iteration, heat_power_meshLevel)
    myProblem.setHeatPowerField(heatPowerField)

    # the boundary conditions :
    if( FECalculation):
        boundaryGroupNames=myProblem.getMesh().getNameOfNodeGroups()
        print(len(boundaryGroupNames), " Boundary Node Group detected : ", boundaryGroupNames)
    else:
        boundaryGroupNames=myProblem.getMesh().getNameOfFaceGroups()
        print(len(boundaryGroupNames), " Boundary Face Group detected : ", boundaryGroupNames)

    # for each boundary we load the boundary field (replace by a loop over the boundaries)
    boundary1_type=solverlab.NeumannDiffusion# or solverlab.DirichletDiffusion
    boundary1_inputfile="../resources/meshSquare";
    boundary1_fieldName="Left temperature";
    boundary1_time_iteration=0# default value is 0
    boundary1_time_sub_iteration=0# default value is 0
    boundary1_meshLevel=0# default value is 0
    print("Boundary ", boundaryGroupNames[3], ", loading field :", boundary1_fieldName, " in file ", boundary1_inputfile)
    boundary1Field=solverlab.Field(boundary1_inputfile, supportOfField, boundary1_fieldName, boundary1_time_iteration, boundary1_time_sub_iteration, boundary1_meshLevel)
    boundary2_type=solverlab.DirichletDiffusion# or solverlab.NeumannDiffusion
    boundary2_inputfile="../resources/meshSquare";
    boundary2_fieldName="Right temperature";
    boundary2_time_iteration=0# default value is 0
    boundary2_time_sub_iteration=0# default value is 0
    boundary2_meshLevel=0# default value is 0
    print("Boundary ", boundaryGroupNames[2], ", loading field :", boundary2_fieldName, " in file ", boundary2_inputfile)
    boundary2Field=solverlab.Field(boundary2_inputfile, supportOfField, boundary2_fieldName, boundary2_time_iteration, boundary2_time_sub_iteration, boundary2_meshLevel)
    boundary3_type=solverlab.NeumannDiffusion# or solverlab.DirichletDiffusion
    boundary3_inputfile="../resources/meshSquare";
    boundary3_fieldName="Top temperature";
    boundary3_time_iteration=0# default value is 0
    boundary3_time_sub_iteration=0# default value is 0
    boundary3_meshLevel=0# default value is 0
    print("Boundary ", boundaryGroupNames[4], ", loading field :", boundary3_fieldName, " in file ", boundary3_inputfile)
    boundary3Field=solverlab.Field(boundary3_inputfile, supportOfField, boundary3_fieldName, boundary3_time_iteration, boundary3_time_sub_iteration, boundary3_meshLevel)
    boundary4_type=solverlab.DirichletDiffusion# or solverlab.NeumannDiffusion
    boundary4_inputfile="../resources/meshSquare";
    boundary4_fieldName="Bottom temperature";
    boundary4_time_iteration=0# default value is 0
    boundary4_time_sub_iteration=0# default value is 0
    boundary4_meshLevel=0# default value is 0
    print("Boundary ", boundaryGroupNames[1], ", loading field :", boundary4_fieldName, " in file ", boundary4_inputfile)
    boundary4Field=solverlab.Field(boundary4_inputfile, supportOfField, boundary4_fieldName, boundary4_time_iteration, boundary4_time_sub_iteration, boundary4_meshLevel)

    # for each boundary we need to know if we want a Neumann or a Dirichlet boundary condition
    if boundary1_type==solverlab.NeumannDiffusion :
        myProblem.setNeumannBoundaryCondition("Left", boundary1Field)
    elif boundary1_type==solverlab.DirichletDiffusion :
        myProblem.setDirichletBoundaryCondition("Left", boundary1Field)
    if boundary2_type==solverlab.NeumannDiffusion :
        myProblem.setNeumannBoundaryCondition("Right", boundary2Field)
    elif boundary2_type==solverlab.DirichletDiffusion :
        myProblem.setDirichletBoundaryCondition("Right", boundary2Field)
    if boundary3_type==solverlab.NeumannDiffusion :
        myProblem.setNeumannBoundaryCondition("Top", boundary3Field)
    elif boundary3_type==solverlab.DirichletDiffusion :
        myProblem.setDirichletBoundaryCondition("Top", boundary3Field);
    if boundary4_type==solverlab.NeumannDiffusion :
        myProblem.setNeumannBoundaryCondition("Bottom", boundary4Field)
    elif boundary4_type==solverlab.DirichletDiffusion :
        myProblem.setDirichletBoundaryCondition("Bottom", boundary4Field);

    # set the numerical method
    myProblem.setTimeScheme( solverlab.Explicit)# Otherwise solverlab.Implicit
    max_nb_its_lin_solver = 50
    myProblem.setLinearSolver(solverlab.GMRES, solverlab.ILU, max_nb_its_lin_solver );

    # computation parameters
    MaxNbOfTimeStep = 3 ;# default value is 10
    freqSave = 1;# default value is 1
    cfl = 0.95;# default value is 1
    maxTime = 100000000;# default value is 10
    precision = 1e-6;# default value is 1e-6
    result_directory="."# default value = current directory

    myProblem.setCFL(cfl);
    myProblem.setPrecision(precision);
    myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
    myProblem.setTimeMax(maxTime);
    myProblem.setFreqSave(freqSave);
    myProblem.setFileName(fileName);
    myProblem.setResultDirectory(result_directory)
    myProblem.setSaveFileFormat(solverlab.MED)#default value is solverlab.VTK

    # evolution
    myProblem.initialize();
    print("Running python test "+ fileName );

    ok = myProblem.run();
    if (ok):
        print( "Python simulation " + fileName + " is successful !" );
        pass
    else:
        print( "Python simulation " + fileName + "  failed ! " );
        pass

    print( "------------ End of simulation !!! -----------" );

    myProblem.terminate();
    return ok

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        FECalculation = bool(int(sys.argv[1]))
        # name of result file
        if( FECalculation):
            fileName = "2DSpherical_FE";# default value is ""
        else:
            fileName = "2DSpherical_FV";# default value is ""
        DiffusionEquation_2DSpherical(FECalculation, fileName)
    else :
        raise ValueError("DiffusionEquation_2DSpherical : missing one argument")
