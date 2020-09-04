#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF décalée du système des ondes 2D :
#                \partial_t p + c^2 \div q = 0 
#                \partial_t q +    \grad p = 0 
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Recherche d'un régime stationnaire
#               Utilisation du schéma décalé implicite sur un maillage cartésien
#               Initialisation par un régime stationnaire
#               Conditions aux limites périodiques
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================

from math import sin, cos, pi
import time, json
import cdmath
import PV_routines
import VTK_routines

test_desc={}

rho0=1000.#reference density
c0=1500.#reference sound speed
p0=rho0*c0*c0#reference pressure
precision=1e-5

def initial_conditions_wave_system_staggered(my_mesh):
    test_desc["Initial_data"]="Constant pressure, divergence free velocity"
    
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    if(not my_mesh.isStructured()):
        raise ValueError("WaveSystemStaggered: the mesh should be structured");

    NxNyNz=my_mesh.getCellGridStructure()
    DxDyDz=my_mesh.getDXYZ()
    dx=DxDyDz[0]
    if(dim>=2):
        dy=DxDyDz[1]
        if(dim==3):
            dz=DxDyDz[2]

    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        Ci=my_mesh.getCell(i)
        x = Ci.x()

        pressure_field[i] = p0
        
        #We take only the normal component of the velocity on a cartesian grid
        #We save the x component from the back face, the y component from the left face and the z component from the bottom face
        #Warning : boundary values should be the same for left and right as well as top and down (front and back in 3D) boundaries
        if(dim==2):
            y = Ci.y()
            velocity_field[i,0] =  sin(2*pi*(x-0.5*dx))*cos(2*pi*y) # value on the left face
            velocity_field[i,1] = -sin(2*pi*(y-0.5*dy))*cos(2*pi*x) # value on the bottom face
            velocity_field[i,2] =  0
        elif(dim==3):
            y = Ci.y()
            z = Ci.z()
            velocity_field[i,0] =    sin(2*pi*(x-0.5*dx))*cos(2*pi*y)*cos(2*pi*z) # value on the back face
            velocity_field[i,1] =    sin(2*pi*(y-0.5*dy))*cos(2*pi*x)*cos(2*pi*z) # value on the left face
            velocity_field[i,2] = -2*sin(2*pi*(z-0.5*dz))*cos(2*pi*x)*cos(2*pi*y) # value on the bottom face
        else :
                raise ValueError("initial_conditions_wave_system_staggered: the 1D mesh not yet implemented");
        
    return pressure_field, velocity_field
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1

    if(not my_mesh.isStructured()):
        raise ValueError("WaveSystemStaggered: the mesh should be structured");

    NxNyNz=my_mesh.getCellGridStructure()
    DxDyDz=my_mesh.getDXYZ()
    
    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    idMoinsJacCL=cdmath.Matrix(nbComp)
    
    if( dim == 1) :    
        nx=NxNyNz[0]
        dx=DxDyDz[0]
            
        if( scaling==0 ):
            for k in range(nbCells):
                implMat.addValue(k,1*nbCells +  k      , -c0*c0*dt/dx)
                implMat.addValue(k,1*nbCells + (k+1)%nx,  c0*c0*dt/dx)
    
                implMat.addValue(  1*nbCells +  k      ,k,  dt/dx)
                implMat.addValue(  1*nbCells + (k+1)%nx,k, -dt/dx)
        else : # scaling >0    
            for k in range(nbCells):
                implMat.addValue(k,1*nbCells +  k      , -c0*dt/dx)
                implMat.addValue(k,1*nbCells + (k+1)%nx,  c0*dt/dx)
    
                implMat.addValue(  1*nbCells +  k      ,k,  c0*dt/dx)
                implMat.addValue(  1*nbCells + (k+1)%nx,k, -c0*dt/dx)
    
    elif( dim == 2) :# k = j*nx+i
        nx=NxNyNz[0]
        ny=NxNyNz[1]
        dx=DxDyDz[0]
        dy=DxDyDz[1]
                
        if( scaling==0 ):
            for k in range(nbCells):
                i = k % nx
                j = k //nx
    
                implMat.addValue(k,1*nbCells + j*nx +  i      ,   -c0*c0*dt/dx)
                implMat.addValue(k,1*nbCells + j*nx + (i+1)%nx,    c0*c0*dt/dx)
    
                implMat.addValue(k,2*nbCells +   j       *nx + i, -c0*c0*dt/dy)
                implMat.addValue(k,2*nbCells + ((j+1)%ny)*nx + i,  c0*c0*dt/dy)
    
                implMat.addValue(  1*nbCells + j*nx +  i      ,  k,  dt/dx)
                implMat.addValue(  1*nbCells + j*nx + (i+1)%nx,  k, -dt/dx)
    
                implMat.addValue(  2*nbCells +   j       *nx + i,k,  dt/dy)
                implMat.addValue(  2*nbCells + ((j+1)%ny)*nx + i,k, -dt/dy)
    
        else :# scaling >0
            for k in range(nbCells):
                i = k % nx
                j = k //nx
    
                implMat.addValue(k,1*nbCells + j*nx +  i      ,   -c0*dt/dx)
                implMat.addValue(k,1*nbCells + j*nx + (i+1)%nx,    c0*dt/dx)
    
                implMat.addValue(k,2*nbCells +   j       *nx + i, -c0*dt/dy)
                implMat.addValue(k,2*nbCells + ((j+1)%ny)*nx + i,  c0*dt/dy)
    
                implMat.addValue(  1*nbCells + j*nx +  i      ,  k,  c0*dt/dx)
                implMat.addValue(  1*nbCells + j*nx + (i+1)%nx,  k, -c0*dt/dx)
    
                implMat.addValue(  2*nbCells +   j       *nx + i,k,  c0*dt/dy)
                implMat.addValue(  2*nbCells + ((j+1)%ny)*nx + i,k, -c0*dt/dy)
    
    elif( dim == 3) :# k = l*nx*ny+j*nx+i
        nx=NxNyNz[0]
        ny=NxNyNz[1]
        nz=NxNyNz[2]
        dx=DxDyDz[0]
        dy=DxDyDz[1]
        dz=DxDyDz[2]
                
        if( scaling==0 ):
            for k in range(nbCells):
                i =  k % nx
                j = (k //nx)%ny 
                l =  k //(nx*ny)
                
                implMat.addValue(k,1*nbCells + l*nx*ny + j*nx +  i      ,  -c0*c0*dt/dx)
                implMat.addValue(k,1*nbCells + l*nx*ny + j*nx + (i+1)%nx,   c0*c0*dt/dx)
    
                implMat.addValue(k,2*nbCells + l*nx*ny +   j       *nx + i, -c0*c0*dt/dy)
                implMat.addValue(k,2*nbCells + l*nx*ny + ((j+1)%ny)*nx + i,  c0*c0*dt/dy)
    
                implMat.addValue(k,3*nbCells +   l*nx*ny        + j*nx + i, -c0*c0*dt/dz)
                implMat.addValue(k,3*nbCells + ((l+1)%nz)*nx*ny + j*nx + i,  c0*c0*dt/dz)
    
                implMat.addValue(  1*nbCells + l*nx*ny + j*nx +  i      ,  k,  dt/dx)
                implMat.addValue(  1*nbCells + l*nx*ny + j*nx + (i+1)%nx,  k, -dt/dx)
    
                implMat.addValue(  2*nbCells + l*nx*ny +   j       *nx + i,k,  dt/dy)
                implMat.addValue(  2*nbCells + l*nx*ny + ((j+1)%ny)*nx + i,k, -dt/dy)
    
                implMat.addValue(  3*nbCells +   l*nx*ny        + j*nx + i,k,  dt/dz)
                implMat.addValue(  3*nbCells + ((l+1)%nz)*nx*ny + j*nx + i,k, -dt/dz)

        else:# scaling >0
            for k in range(nbCells):
                i =  k % nx
                j = (k //nx)%ny 
                l =  k //(nx*ny)
                
                implMat.addValue(k,1*nbCells + l*nx*ny + j*nx +  i      ,  -c0*dt/dx)
                implMat.addValue(k,1*nbCells + l*nx*ny + j*nx + (i+1)%nx,   c0*dt/dx)
    
                implMat.addValue(k,2*nbCells + l*nx*ny +   j       *nx + i, -c0*dt/dy)
                implMat.addValue(k,2*nbCells + l*nx*ny + ((j+1)%ny)*nx + i,  c0*dt/dy)
    
                implMat.addValue(k,3*nbCells +   l*nx*ny        + j*nx + i, -c0*dt/dz)
                implMat.addValue(k,3*nbCells + ((l+1)%nz)*nx*ny + j*nx + i,  c0*dt/dz)
    
                implMat.addValue(  1*nbCells + l*nx*ny + j*nx +  i      ,  k,  c0*dt/dx)
                implMat.addValue(  1*nbCells + l*nx*ny + j*nx + (i+1)%nx,  k, -c0*dt/dx)
    
                implMat.addValue(  2*nbCells + l*nx*ny +   j       *nx + i,k,  c0*dt/dy)
                implMat.addValue(  2*nbCells + l*nx*ny + ((j+1)%ny)*nx + i,k, -c0*dt/dy)
    
                implMat.addValue(  3*nbCells +   l*nx*ny        + j*nx + i,k,  c0*dt/dz)
                implMat.addValue(  3*nbCells + ((l+1)%nz)*nx*ny + j*nx + i,k, -c0*dt/dz)

    #Add the identity matrix on the diagonal
    #divMat.diagonalShift(1)#only after  filling all coefficients
    for j in range(nbCells*nbComp):
        implMat.addValue(j,j,1)#/(c0*c0)

    return implMat

def WaveSystemStaggered(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,scaling,test_bc):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    iterGMRESMax=50
    
    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial data …")
    pressure_field, velocity_field = initial_conditions_wave_system_staggered(my_mesh)
    initial_pressure, initial_velocity = initial_conditions_wave_system_staggered(my_mesh)

    for k in range(nbCells):
        Un[k + 0*nbCells] =      initial_pressure[k]
        Un[k + 1*nbCells] = rho0*initial_velocity[k,0] # value on the left face
        Un[k + 2*nbCells] = rho0*initial_velocity[k,1] # value on the bottom face
        if(dim==3):
            Un[k + 3*nbCells] = rho0*initial_velocity[k,2]
    if( scaling>0):
        Vn = Un.deepCopy()
        for k in range(nbCells):
            Vn[k] = Vn[k]/c0
            
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity");
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_initial")
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_initial")

    total_pressure_initial=pressure_field.integral()#For conservation test later
    total_velocity_initial=velocity_field.integral()#For conservation test later
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0
    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling)

    if( scaling==0):
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")
    else:
        LS=cdmath.LinearSolver(divMat,Vn,iterGMRESMax, precision, "GMRES","ILU")
    LS.setComputeConditionNumber()

    test_desc["Linear_solver_algorithm"]=LS.getNameOfMethod()
    test_desc["Linear_solver_preconditioner"]=LS.getNameOfPc()
    test_desc["Linear_solver_precision"]=LS.getTolerance()
    test_desc["Linear_solver_maximum_iterations"]=LS.getNumberMaxOfIter()
    test_desc["Numerical_parameter_space_step"]=dx_min
    test_desc["Numerical_parameter_time_step"]=dt
    test_desc["Linear_solver_with_scaling"]=scaling

    test_desc['Linear_system_max_actual_iterations_number']=0
    test_desc["Linear_system_max_actual_error"]=0
    test_desc["Linear_system_max_actual_condition number"]=0

    print("Starting computation of the linear wave system with a pseudo staggered scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        dUn=Un.deepCopy()
        if( scaling==0):
            LS.setSndMember(Un)
            Un=LS.solve()
        else:#( scaling > 0)
            LS.setSndMember(Vn)
            Vn=LS.solve()
            Un = Vn.deepCopy()
            for k in range(nbCells):
                Un[k] = c0*Vn[k]
            
        if(not LS.getStatus()):
            print "Linear system did not converge ", iterGMRES, " GMRES iterations"
            raise ValueError("Pas de convergence du système linéaire");
        dUn-=Un
        
        test_desc["Linear_system_max_actual_iterations_number"]=max(LS.getNumberOfIter(),test_desc["Linear_system_max_actual_iterations_number"])
        test_desc["Linear_system_max_actual_error"]=max(LS.getResidu(),test_desc["Linear_system_max_actual_error"])
        test_desc["Linear_system_max_actual_condition number"]=max(LS.getConditionNumber(),test_desc["Linear_system_max_actual_condition number"])

        max_dp=0 ;        max_dq=0
        for k in range(nbCells):
            max_dp = max(max_dp,abs(dUn[k]))
            for i in range(dim):
                max_dq=max(max_dq,abs(dUn[k+(1+i)*nbCells]))
                
        isStationary= max_dp/p0<precision and max_dq/rho0<precision

        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            print "Variation temporelle relative : pressure ", max_dp/p0 ,", velocity ", max_dq/rho0
            print "Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations"

            delta_press=0
            delta_v=cdmath.Vector(dim)
            for k in range(nbCells):
                pressure_field[k]=Un[k]
                velocity_field[k,0]=Un[k+1*nbCells]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k+2*nbCells]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k+3*nbCells]/rho0
                if (abs(initial_pressure[k]-pressure_field[k])>delta_press):
                    delta_press=abs(initial_pressure[k]-pressure_field[k])
                if (abs(initial_velocity[k,0]-velocity_field[k,0])>delta_v[0]):
                    delta_v[0]=abs(initial_velocity[k,0]-velocity_field[k,0])
                if (abs(initial_velocity[k,1]-velocity_field[k,1])>delta_v[1]):
                    delta_v[1]=abs(initial_velocity[k,1]-velocity_field[k,1])
                if(dim==3):
                    if (abs(initial_velocity[k,2]-velocity_field[k,2])>delta_v[2]):
                        delta_v[2]=abs(initial_velocity[k,2]-velocity_field[k,2])
                
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity",False);

            print "Ecart au stationnaire exact : error_p= ",delta_press/p0," error_||u||= ",delta_v.maxVector()[0]
            print
    print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
    print "Variation temporelle relative : pressure ", max_dp/p0 ,", velocity ", max_dq/rho0 
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        print "Mass loss: ", (total_pressure_initial-pressure_field.integral()).norm()/p0, " precision required= ", precision
        print "Momentum loss: ", (total_velocity_initial-velocity_field.integral()).norm()/velocity_field.normL1().norm(), " precision required= ", precision
        assert (total_pressure_initial-pressure_field.integral()).norm()/p0<precision
        if(test_bc=="Periodic"):
            assert (total_velocity_initial-velocity_field.integral()).norm()<2*precision
        print "------------------------------------------------------------------------------------"

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        if(dim==2):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        elif(dim==3):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,0,0],[1,1,1], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,0,0],[1,1,1], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat")
        
        return delta_press/p0, delta_v.maxVector()[0], nbCells, time, it, velocity_field.getNormEuclidean().max(), diag_data_press, diag_data_vel,test_desc["Linear_system_max_actual_condition number"]
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh,meshName,resolution,scaling, meshType, testColor,cfl):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    test_name="Resolution of the Wave system in dimension " +str( my_mesh.getMeshDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    test_name_comment="New scheme for low Mach flows"
    test_model="Wave system"
    test_method="Staggered"
    test_initial_data="Constant pressure, divergence free velocity"
    test_bc="Periodic"
    print test_name
    print "Numerical method : ", test_method
    print "Initial data : ", test_initial_data
    print "Boundary conditions : ",test_bc
    print "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells"
    if( scaling>0):
        print "Use of scaling strategy for a better preconditioning"

    # Problem data
    tmax = 1000.
    ntmax = 35000
    output_freq = 1000

    error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel, cond_number = WaveSystemStaggered(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,scaling,test_bc)
    end = time.time()

    test_desc["Global_name"]=test_name
    test_desc["Global_comment"]=test_name_comment
    test_desc["PDE_model"]=test_model
    test_desc["PDE_is_stationary"]=False
    test_desc["PDE_search_for_stationary_solution"]=True
    test_desc["Numerical_method_name"]=test_method
    test_desc["Numerical_method_space_discretization"]="Finite volumes"
    test_desc["Numerical_method_time_discretization"]="Implicit"
    test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
    test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
    test_desc["Mesh_is_unstructured"]=True
    test_desc["Mesh_cell_type"]=my_mesh.getElementTypesNames()
    test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfCells()
    test_desc["Mesh_max_number_of_neighbours"]=10
    test_desc["Geometry"]="Square"
    test_desc["Boundary_conditions"]=test_bc
    test_desc["Initial_data"]=test_initial_data
    test_desc["Part_of_mesh_convergence_analysis"]=True
    test_desc["Numerical_parameter_cfl"]=cfl
    test_desc["Simulation_parameter_maximum_time_step"]=ntmax
    test_desc["Simulation_parameter_maximum_time"]=tmax
    test_desc["Simulation_output_frequency"]=output_freq
    test_desc["Simulation_final_time_after_run"]=t_final
    test_desc["Simulation_final_number_of_time_steps_after_run"]=ndt_final
    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Absolute_error"]=max(error_p*p0,error_u)
    test_desc["Relative_error"]=max(error_p,error_u)

    with open('test_WaveSystem'+str(my_mesh.getMeshDimension())+'DStaggered_'+meshName+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)
    
    return error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel, end - start, cond_number

def solve_file( filename,meshName, resolution,scaling, meshType, testColor,cfl):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution,scaling, meshType, testColor,cfl)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0,1,20,0,1,20)
    cfl=0.5
    solve(M1,"SquareWithSquares",100,2,"Regular squares","Green",cfl)
