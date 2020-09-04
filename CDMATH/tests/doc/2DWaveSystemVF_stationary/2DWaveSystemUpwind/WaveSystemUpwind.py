#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import time, json
import cdmath
import PV_routines
import VTK_routines

test_desc={}

rho0=1000#reference density
c0=1500#reference sound speed
p0=rho0*c0*c0#reference pressure
precision=1e-5

def initial_conditions_wave_system(my_mesh):
    test_desc["Initial_data"]="Constant pressure, divergence free velocity"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()
        z = my_mesh.getCell(i).z()

        pressure_field[i] = p0
        if(dim==2):
            velocity_field[i,0] =  sin(pi*x)*cos(pi*y)
            velocity_field[i,1] = -sin(pi*y)*cos(pi*x)
            velocity_field[i,2] = 0
        if(dim==3):
            velocity_field[i,0] =    sin(pi*x)*cos(pi*y)*cos(pi*z)
            velocity_field[i,1] =    sin(pi*y)*cos(pi*x)*cos(pi*z)
            velocity_field[i,2] = -2*sin(pi*z)*cos(pi*x)*cos(pi*y)
        
    return pressure_field, velocity_field

def jacobianMatrices(normal, coeff,scaling):
    test_desc["Numerical_method_name"]="Upwind"
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)
    absA=cdmath.Matrix(dim+1,dim+1)

    absA[0,0]=c0*coeff
    for i in range(dim):
        for j in range(dim):
            absA[i+1,j+1]=c0*normal[i]*normal[j]*coeff
        if( scaling==0):
            A[0,i+1]=c0*c0*normal[i]*coeff
            A[i+1,0]=      normal[i]*coeff
        elif( scaling==1):
            A[0,i+1]=      normal[i]*coeff
            A[i+1,0]=      normal[i]*coeff
        else:
            A[0,i+1]=   c0*normal[i]*coeff
            A[i+1,0]=   c0*normal[i]*coeff
       
    return (A-absA)/2
    
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling,test_bc):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1
    normal=cdmath.Vector(dim)

    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    idMoinsJacCL=cdmath.Matrix(nbComp)
    
    for j in range(nbCells):#On parcourt les cellules
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

            Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure(),scaling);

            cellAutre =-1
            if ( not Fk.isBorder()) :
                # hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if (Fk.getCellsId()[0] == j) :
                    # hypothese verifiée 
                    cellAutre = Fk.getCellsId()[1];
                elif(Fk.getCellsId()[1] == j) :
                    # hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                else :
                    raise ValueError("computeFluxes: problem with mesh, unknown cel number")
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                if( test_bc=="Periodic" and Fk.getGroupName() != "Wall" and Fk.getGroupName() != "Paroi" and Fk.getGroupName() != "Neumann"):#Periodic boundary condition unless Wall/Neumann specified explicitly
                    test_desc["Boundary_conditions"]="Periodic"
                    indexFP = my_mesh.getIndexFacePeriodic(indexFace)
                    Fp = my_mesh.getFace(indexFP)
                    cellAutre = Fp.getCellsId()[0]
                    
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif( test_bc=="Wall" or Fk.getGroupName() == "Wall" or Fk.getGroupName() == "Paroi"):#Wall boundary condition
                    test_desc["Boundary_conditions"]="Wall"
                    v=cdmath.Vector(dim+1)
                    for i in range(dim) :
                        v[i+1]=normal[i]
                    idMoinsJacCL=v.tensProduct(v)*2
                    
                    implMat.addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL)
                    
                elif(test_bc!="Neumann" and Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print Fk.getGroupName()
                    raise ValueError("computeFluxes: Unknown boundary condition name");
                
    return implMat

def WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,scaling,test_bc):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    isImplicit=False
    
    #iteration vectors
    Un=cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_wave_system(my_mesh)
    initial_pressure, initial_velocity = initial_conditions_wave_system(my_mesh)

    for k in range(nbCells):
        Un[k*(dim+1)+0] =     initial_pressure[k]
        Un[k*(dim+1)+1] =rho0*initial_velocity[k,0]
        Un[k*(dim+1)+2] =rho0*initial_velocity[k,1]
        if(dim==3):
            Un[k*(dim+1)+3] =rho0*initial_velocity[k,2]
    if( scaling==1):
        Vn = Un.deepCopy()
        for k in range(nbCells):
            Vn[k*(dim+1)+0] = Vn[k*(dim+1)+0]/(c0*c0)
    elif( scaling==2):
        Vn = Un.deepCopy()
        for k in range(nbCells):
            Vn[k*(dim+1)+0] = Vn[k*(dim+1)+0]/c0
            
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity");
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_initial")
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_initial")

    total_pressure_initial=pressure_field.integral()#For conservation test later
    total_velocity_initial=velocity_field.integral()#For conservation test later
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling,test_bc)
    if(isImplicit):
        #Adding the identity matrix on the diagonal
        if( scaling==0 or  scaling==2):
            divMat.diagonalShift(1)#only after  filling all coefficients
            #for j in range(nbCells*(dim+1)):
            #    divMat.addValue(j,j,1)
        else:
            for j in range(nbCells):
                divMat.addValue(j*(dim+1),j*(dim+1),1/(c0*c0))#/(c0*c0)
                for i in range(dim):
                    divMat.addValue(j*(dim+1)+1+i,j*(dim+1)+1+i,1)
        
        iterGMRESMax=50
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

    print("Starting computation of the linear wave system with an Upwind scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            if( scaling==0):
                LS.setSndMember(Un)
            else:
                LS.setSndMember(Vn)
            if( scaling<2):
                Un=LS.solve();
                if( scaling==1):
                    Vn = Un.deepCopy()
                    for k in range(nbCells):
                        Vn[k*(dim+1)+0] = Vn[k*(dim+1)+0]/(c0*c0)
            else:#( scaling==2)
                Vn=LS.solve()
                Un = Vn.deepCopy()
                for k in range(nbCells):
                    Un[k*(dim+1)+0] = c0*Vn[k*(dim+1)+0]
            if(not LS.getStatus()):
                print "Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations"
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

            test_desc["Linear_system_max_actual_iterations_number"]=max(LS.getNumberOfIter(),test_desc["Linear_system_max_actual_iterations_number"])
            test_desc["Linear_system_max_actual_error"]=max(LS.getResidu(),test_desc["Linear_system_max_actual_error"])
            test_desc["Linear_system_max_actual_condition number"]=max(LS.getConditionNumber(),test_desc["Linear_system_max_actual_condition number"])

        else:
            dUn=divMat*Un
            Un-=dUn
        
        maxVector=dUn.maxVector(dim+1)
        isStationary= maxVector[0]/p0<precision and maxVector[1]/rho0<precision and maxVector[2]/rho0<precision;
        if(dim==3):
            isStationary=isStationary and maxVector[3]/rho0<precision
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
            if(isImplicit):
                print "Linear system converged in ", iterGMRES, " GMRES iterations"

            delta_press=0
            delta_v=cdmath.Vector(dim)
            for k in range(nbCells):
                pressure_field[k]=Un[k*(dim+1)+0]
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0
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
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity",False);

            print "Ecart au stationnaire exact : error_p= ",delta_press/p0," error_||u||= ",delta_v.maxVector()[0]
            print
    print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
    print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
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
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        if(dim==2):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        elif(dim==3):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,0,0],[1,1,1], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,0,0],[1,1,1], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat")
        
        return delta_press/p0, delta_v.maxVector()[0], nbCells, time, it, velocity_field.getNormEuclidean().max(), diag_data_press, diag_data_vel
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh, meshName, resolution, scaling, meshType, testColor, cfl, test_bc):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    test_name="Resolution of the Wave system in dimension " +str( my_mesh.getSpaceDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    test_name_comment="Classical characteristic based scheme"
    test_model="Wave system"
    test_method="Upwind"
    test_initial_data="Constant pressure, divergence free velocity"

    print test_name
    print "Numerical method : ", test_method
    print "Initial data : ", test_initial_data
    print "Boundary conditions : ",test_bc
    print "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells"
    if( scaling):
        print "Use of scaling strategy for better preconditioning"
    
    # Problem data
    tmax = 10000.
    ntmax = 30000
    output_freq = 10000

    error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel = WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution, scaling,test_bc)
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
    test_desc["Mesh_cell_type"]=my_mesh.getElementTypes()
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

    with open('test_WaveSystem'+str(my_mesh.getMeshDimension())+'DUpwind_'+meshName+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel, end - start

def solve_file( filename,meshName, resolution, scaling, meshType, testColor,cfl,test_bc):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution, scaling, meshType, testColor,cfl,test_bc)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0.,1.,20,0.,1.,20,0)
    cfl=0.5
    solve(M1,"SquareRegularTriangles",100,"Regular triangles","Green",cfl,"Periodic")

    M2=cdmath.Mesh(0.,1.,10,0.,1.,10,0.,1.,10,6)
    cfl=1./3
    solve(M2,"CubeRegularTetrahedra",100,"Regular tetrahedra","Green",cfl,"Wall")
