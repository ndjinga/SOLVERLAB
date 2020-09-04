#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système des ondes 2D avec terme source régulier
#                \partial_t p + c^2 \div q = c^2*(sign(x-1./2)*(y-1./2)*abs(y-1./2) + (x-1./2)*abs(x-1./2)*sign(y-1./2))
#                \partial_t q +    \grad p =-q
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Recherche d'un régime stationnaire
#               Utilisation du schéma centré implicite sur un maillage général
#               Initialisation par un champ stationnaire
#               Conditions aux limites périodiques
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sin, cos, pi
import time, json
import cdmath
import PV_routines
import VTK_routines
from numpy import sign

test_desc={}

rho0=1000.#reference density
c0=1500.#reference sound speed
p0=rho0*c0*c0#reference pressure
precision=1e-5

def source_term_and_stat_solution_wave_system(my_mesh):
    test_desc["Source_term"]="True"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    source_vector = cdmath.Vector(nbCells*(dim+1))

    stat_pressure_field = cdmath.Field("Pressure", cdmath.CELLS, my_mesh, 1)
    stat_momentum_field = cdmath.Field("Momentum", cdmath.CELLS, my_mesh, 3)

    for k in range(nbCells):
        x = my_mesh.getCell(k).x()
        y = my_mesh.getCell(k).y()
        z = my_mesh.getCell(k).z()

        if(dim==1):
            source_vector[k*(dim+1)+0] = -2*c0*c0*sign(x-1./2)
            source_vector[k*(dim+1)+1] = 0

            stat_pressure_field[k]   =    abs(x-1./2)*(x-1./2)
            stat_momentum_field[k,0] = -2*abs(x-1./2)
            stat_momentum_field[k,1] = 0
            stat_momentum_field[k,2] = 0
        elif(dim==2):
            source_vector[k*(dim+1)+0] = -2*c0*c0*(sign(x-1./2)*(y-1./2)*abs(y-1./2) + (x-1./2)*abs(x-1./2)*sign(y-1./2))
            source_vector[k*(dim+1)+1] = 0
            source_vector[k*(dim+1)+2] = 0

            stat_pressure_field[k]   =    (x-1./2)*abs(x-1./2)*(y-1./2)*abs(y-1./2)
            stat_momentum_field[k,0] = -2*         abs(x-1./2)*(y-1./2)*abs(y-1./2)
            stat_momentum_field[k,1] = -2*(x-1./2)*abs(x-1./2)         *abs(y-1./2)
            stat_momentum_field[k,2] = 0
        elif(dim==3):
            source_vector[k*(dim+1)+0] = -2*c0*c0*(sign(x-1./2)*(y-1./2)*abs(y-1./2)*(z-1./2)*abs(z-1./2)+(x-1./2)*abs(x-1./2)*sign(y-1./2)*(z-1./2)*abs(z-1./2)+(x-1./2)*abs(x-1./2)*(y-1./2)*abs(y-1./2)*sign(z-1./2))
            source_vector[k*(dim+1)+1] = 0
            source_vector[k*(dim+1)+2] = 0
            source_vector[k*(dim+1)+3] = 0
        
            stat_pressure_field[k]   =       (x-1./2)*abs(x-1./2)*(y-1./2)*abs(y-1./2)*(z-1./2)*abs(z-1./2)
            stat_momentum_field[k,0] = -2*abs(x-1./2)*(y-1./2)*abs(y-1./2)*(z-1./2)*abs(z-1./2)
            stat_momentum_field[k,1] = -2*(x-1./2)*abs(x-1./2)*abs(y-1./2)*(z-1./2)*abs(z-1./2)
            stat_momentum_field[k,2] = -2*(x-1./2)*abs(x-1./2)*(y-1./2)*abs(y-1./2)*abs(z-1./2)

    return source_vector, stat_pressure_field, stat_momentum_field

def jacobianMatrices(normal, coeff, scaling):
    test_desc["Numerical_method_name"]="Centered scheme"
    
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)

    for i in range(dim):
        if( scaling==0):
            A[0,i+1]=c0*c0*normal[i]*coeff
            A[i+1,0]=      normal[i]*coeff
        elif( scaling==1):
            A[0,i+1]=      normal[i]*coeff
            A[i+1,0]=      normal[i]*coeff
        else:
            A[0,i+1]=   c0*normal[i]*coeff
            A[i+1,0]=   c0*normal[i]*coeff
    
    return A/2
    
    
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
                    raise ValueError("computeFluxes: problem with mesh, unknown cell number")
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                if( test_bc=="Periodic" and Fk.getGroupName() != "Wall" and Fk.getGroupName() != "Paroi" and Fk.getGroupName() != "Neumann"):#Periodic boundary condition unless Wall/Neumann specified explicitly
                    test_desc["Boundary_conditions"]="Periodic"
                    
                    indexFP = my_mesh.getIndexFacePeriodic(indexFace, my_mesh.getName()== "squareWithBrickWall", my_mesh.getName()== "squareWithHexagons")
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
    iterGMRESMax=50
    
    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(meshName.find("square")==-1 and meshName.find("Square")==-1 and meshName.find("cube")==-1 and meshName.find("Cube")==-1):
        print "Mesh name : ", meshName
        raise ValueError("Mesh name should contain substring square or cube to use wave system with stiff source term")
    else:
        S, stat_pressure, stat_momentum=source_term_and_stat_solution_wave_system(my_mesh)
        pressure_field = stat_pressure.deepCopy()
        momentum_field = stat_momentum.deepCopy()

    for k in range(nbCells):
        Un[k*(dim+1)+0] =     pressure_field[k]
        Un[k*(dim+1)+1] =     momentum_field[k,0]
        Un[k*(dim+1)+2] =     momentum_field[k,1]
        if(dim==3):
            Un[k*(dim+1)+3] =momentum_field[k,2]

    if( scaling==1):
        Vn = Un.deepCopy()
        for k in range(nbCells):
            Vn[k*(dim+1)+0] = Vn[k*(dim+1)+0]/(c0*c0)
            S[ k*(dim+1)+0] = S[ k*(dim+1)+0]/(c0*c0)
    elif( scaling==2):
        Vn = Un.deepCopy()
        for k in range(nbCells):
            Vn[k*(dim+1)+0] = Vn[k*(dim+1)+0]/c0
            S[ k*(dim+1)+0] = S[ k*(dim+1)+0]/c0

    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure");
    momentum_field.setTime(time,it);
    momentum_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum");
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_initial")
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum"+'_0.vtu',"Momentum",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum_initial")

    total_pressure_initial=pressure_field.integral()#For conservation test later
    total_momentum_initial=momentum_field.integral()#For conservation test later
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0
    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling,test_bc)
    #Adding the momentum friction term
    for j in range(nbCells):
        for i in range(dim):
            divMat.addValue(j*(dim+1)+1+i,j*(dim+1)+1+i,dt)

    #Add the identity matrix on the diagonal
    if( scaling==0 or  scaling==2):
        divMat.diagonalShift(1)#only after  filling all coefficients
    else:
        for j in range(nbCells):
            divMat.addValue(j*(dim+1),j*(dim+1),1/(c0*c0))#/(c0*c0)
            for i in range(dim):
                divMat.addValue(j*(dim+1)+1+i,j*(dim+1)+1+i,1)

    if( scaling==0):
        LS=cdmath.LinearSolver(divMat,Un+S*dt,iterGMRESMax, precision, "GMRES","ILU")
    else:
        LS=cdmath.LinearSolver(divMat,Vn+S*dt,iterGMRESMax, precision, "GMRES","ILU")
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

    print("Starting computation of the linear wave system with a centered scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        dUn=Un.deepCopy()
        if( scaling==0):
            LS.setSndMember(Un+S*dt)
        else:
            LS.setSndMember(Vn+S*dt)
        if( scaling<2):
            Un=LS.solve()
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
            print "Linear system did not converge ", iterGMRES, " GMRES iterations"
            raise ValueError("Pas de convergence du système linéaire");
        dUn-=Un
        
        test_desc["Linear_system_max_actual_iterations_number"]=max(LS.getNumberOfIter(),test_desc["Linear_system_max_actual_iterations_number"])
        test_desc["Linear_system_max_actual_error"]=max(LS.getResidu(),test_desc["Linear_system_max_actual_error"])
        test_desc["Linear_system_max_actual_condition number"]=max(LS.getConditionNumber(),test_desc["Linear_system_max_actual_condition number"])

        maxVector=dUn.maxVector(dim+1)

        isStationary = maxVector[0]  < precision and maxVector[1]/rho0<precision and maxVector[2]/rho0<precision;

        if(dim==3):
            isStationary=isStationary and maxVector[3]/rho0<precision
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            print "Variation temporelle relative : pressure ", maxVector[0]    ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
            print "Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations"

            delta_press=0
            delta_q=cdmath.Vector(dim)
            for k in range(nbCells):
                pressure_field[k]  =Un[k*(dim+1)+0]
                momentum_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    momentum_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        momentum_field[k,2]=Un[k*(dim+1)+3]/rho0
                if (abs(stat_pressure[k]-pressure_field[k])>delta_press):
                    delta_press=abs(stat_pressure[k]-pressure_field[k])
                if (abs(stat_momentum[k,0]-momentum_field[k,0])>delta_q[0]):
                    delta_q[0]=abs(stat_momentum[k,0]-momentum_field[k,0])
                if (abs(stat_momentum[k,1]-momentum_field[k,1])>delta_q[1]):
                    delta_q[1]=abs(stat_momentum[k,1]-momentum_field[k,1])
                if(dim==3):
                    if (abs(stat_momentum[k,2]-momentum_field[k,2])>delta_q[2]):
                        delta_q[2]=abs(stat_momentum[k,2]-momentum_field[k,2])
                
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure",False);
            momentum_field.setTime(time,it);
            momentum_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum",False);

            print "Ecart au stationnaire exact : error_p= ",delta_press   ," error_||q||= ",delta_q.maxVector()[0]
            print
    print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
    print "Variation temporelle relative : pressure ", maxVector[0]    ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        print "------------------------------------------------------------------------------------"

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat");
        momentum_field.setTime(time,0);
        momentum_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum_Stat");

        #Postprocessing : Extraction of the diagonal data
        if(dim==2):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(momentum_field,[0,1,0],[1,0,0], resolution)    
        elif(dim==3):
            diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,0,0],[1,1,1], resolution)    
            diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(momentum_field,[0,0,0],[1,1,1], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum_Stat"+'_0.vtu',"Momentum",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_momentum_Stat")
        
        return delta_press   , delta_q.maxVector()[0], nbCells, time, it, momentum_field.getNormEuclidean().max(), diag_data_press, diag_data_vel,test_desc["Linear_system_max_actual_condition number"]
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh,meshName,resolution,scaling, meshType, testColor,cfl,test_bc="Periodic"):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    test_name="Resolution of the Wave system in dimension " +str( my_mesh.getMeshDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    test_name_comment="No upwinding"
    test_model="Wave system"
    test_method="Centered"
    test_initial_data="Stiff stationary state"
    print test_name
    print "Numerical method : ", test_method
    print "Initial data : ", test_initial_data
    print "Boundary conditions : ",test_bc
    print "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells"
    if( scaling>0):
        print "Use of scaling strategy for better preconditioning"

    # Problem data
    tmax = 1000.
    ntmax = 35000
    output_freq = 10000

    error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel, cond_number = WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,scaling,test_bc)
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
    if(meshName.find("square")>-1 or meshName.find("Square")>-1):
        test_desc["Geometry"]="Square"
    elif(meshName.find("disk")>-1) or meshName.find("Disk")>-1:
        test_desc["Geometry"]="Disk"
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

    with open('test_WaveSystem'+str(my_mesh.getMeshDimension())+'DCentered_'+meshName+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)
    
    return error_p, error_u, nbCells, t_final, ndt_final, max_vel, diag_data_press, diag_data_vel, end - start, cond_number

def solve_file( filename,meshName, resolution,scaling, meshType, testColor,cfl,test_bc):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution,scaling, meshType, testColor,cfl,test_bc)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0.,1.,20,0.,1.,20)
    cfl=0.5
    scaling=0
    solve(M1,"SquareRegularSquares",100,scaling,"Regular_squares","Green",cfl,"Periodic")

    M2=cdmath.Mesh(0.,1.,10,0.,1.,10,0.,1.,10,6)
    cfl=1./3
    scaling=2
    solve(M2,"CubeRegularTetrahedra",100,scaling,"Regular_tetrahedra","Green",cfl,"Wall")
