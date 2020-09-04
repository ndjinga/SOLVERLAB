#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import time, json
import cdmath
import PV_routines
import VTK_routines

test_desc={}

precision=1e-5

def initial_conditions_transport_equation(my_mesh):
    test_desc["Initial_data"]="Shock"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    initial_field = cdmath.Field("unknown", cdmath.CELLS, my_mesh, 1)

    rayon = 0.15
    xcentre = 0.5
    ycentre = 0.5
    zcentre = 0.5

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()
        z = my_mesh.getCell(i).z()

        valX = (x - xcentre) * (x - xcentre)
        
        if(dim==1):
            val=valX
        if(dim==2):
            valY = (y - ycentre) * (y - ycentre)
            val =  sqrt(valX + valY)
        if(dim==3):
            valY = (y - ycentre) * (y - ycentre)
            valZ = (z - zcentre) * (z - zcentre)
            val =  sqrt(valX + valY + valZ)

        if val < rayon:
            initial_field[i] = 1
            pass
        else:
            initial_field[i] = 0
            pass
        pass
        
    return initial_field

def source_term_and_stat_solution_transport_equation(my_mesh,velocity):
    test_desc["Source_term"]="True"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    source_vector = cdmath.Vector(nbCells)

    stat_field = cdmath.Field("Stationary field", cdmath.CELLS, my_mesh, 1)

    for k in range(nbCells):
        x = my_mesh.getCell(k).x()
        y = my_mesh.getCell(k).y()
        z = my_mesh.getCell(k).z()

        if(dim==1):
            source_vector[k] = 2*pi*velocity[0]*cos(2*pi*x)
            stat_field[k]    =                  sin(2*pi*x)
        elif(dim==2):
            source_vector[k*(dim+1)+0] = 2*pi*(velocity[0]*cos(2*pi*x)*sin(2*pi*y)+velocity[1]*sin(2*pi*x)*cos(2*pi*y))
            stat_field[k]   =                              sin(2*pi*x)*sin(2*pi*y)
        elif(dim==3):
            source_vector[k*(dim+1)+0] = 2*pi*(velocity[0]*cos(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)+velocity[1]*sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z)+velocity[2]*sin(2*pi*x)*sin(2*pi*y)*cos(2*pi*z))
            stat_field[k]   =                              sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

    return source_vector, stat_field

def upwinding_coeff(normal, coeff, velocity):
    test_desc["Numerical_method_name"]="Upwind"
    dim=normal.size()
    
    if(velocity*normal<0.):
        return velocity*normal*coeff
    else:
        return 0.
    
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,test_bc,velocity):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=1
    normal=cdmath.Vector(dim)

    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    for j in range(nbCells):#On parcourt les cellules
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

            Am=upwinding_coeff( normal,dt*Fk.getMeasure()/Cj.getMeasure(),velocity);

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
                if( test_bc=="Periodic" and Fk.getGroupName() != "Neumann"):#Periodic boundary condition unless Wall/Neumann specified explicitly
                    test_desc["Boundary_conditions"]="Periodic"
                    indexFP = my_mesh.getIndexFacePeriodic(indexFace, my_mesh.getName()== "squareWithBrickWall", my_mesh.getName()== "squareWithHexagons")
                    Fp = my_mesh.getFace(indexFP)
                    cellAutre = Fp.getCellsId()[0]
                    
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif(test_bc!="Neumann" and Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print( Fk.getGroupName() )
                    raise ValueError("computeFluxes: Unknown boundary condition name");
                
    return implMat

def TransportEquationVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,test_bc,velocity,isImplicit,with_source=False):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
    #iteration vectors
    Un =cdmath.Vector(nbCells)
    dUn=cdmath.Vector(nbCells)
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(with_source):
        unknown_field = cdmath.Field("Unknown", cdmath.CELLS, my_mesh, 1)
        for k in range(nbCells):# fields initialisation
            unknown_field[k]   = 0
        S, stat_field=source_term_and_stat_solution_transport_equation(my_mesh)
    else:#The initial datum is a stationary field
        stat_field = initial_conditions_transport_equation(my_mesh)
        for k in range(nbCells):
            Un[k] =     stat_field[k]
        unknown_field = initial_conditions_transport_equation(my_mesh)
        S = cdmath.Vector(nbCells)#source term is zero
            
    #sauvegarde de la donnée initiale
    stat_field.setTime(time,it);
    stat_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName);
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("TransportEquation"+str(dim)+"DUpwind"+meshName+'_0.vtu',"unknown",'CELLS',"TransportEquation"+str(dim)+"DUpwind"+meshName+"_initial")

    total_mass_initial=unknown_field.integral()#For conservation test later
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / velocity.norm()

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,test_bc,velocity)
    if(isImplicit):
        #Adding the identity matrix on the diagonal
        divMat.diagonalShift(1)#only after  filling all coefficients
        #for j in range(nbCells*(dim+1)):
        #    divMat.addValue(j,j,1)
        
        iterGMRESMax=50
        LS=cdmath.LinearSolver(divMat,Un+S*dt,iterGMRESMax, precision, "GMRES","ILU")
        LS.setComputeConditionNumber()
        
        test_desc["Linear_solver_algorithm"]=LS.getNameOfMethod()
        test_desc["Linear_solver_preconditioner"]=LS.getNameOfPc()
        test_desc["Linear_solver_precision"]=LS.getTolerance()
        test_desc["Linear_solver_maximum_iterations"]=LS.getNumberMaxOfIter()
        test_desc["Numerical_parameter_space_step"]=dx_min
        test_desc["Numerical_parameter_time_step"]=dt
    
        test_desc['Linear_system_max_actual_iterations_number']=0
        test_desc["Linear_system_max_actual_error"]=0
        test_desc["Linear_system_max_actual_condition number"]=0

    print("Starting computation of the linear transport equation with an Upwind scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            LS.setSndMember(Un+S*dt)
            Un=LS.solve();
            if(not LS.getStatus()):
                print( "Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations" )
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

            test_desc["Linear_system_max_actual_iterations_number"]=max(LS.getNumberOfIter(),test_desc["Linear_system_max_actual_iterations_number"])
            test_desc["Linear_system_max_actual_error"]=max(LS.getResidu(),test_desc["Linear_system_max_actual_error"])
            test_desc["Linear_system_max_actual_condition number"]=max(LS.getConditionNumber(),test_desc["Linear_system_max_actual_condition number"])

        else:
            dUn=divMat*Un-S*dt
            Un-=dUn
        
        maxVector=dUn.maxVector(1)
        isStationary= maxVector[0]<precision 

        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            print("Variation temporelle relative : ", maxVector[0] )
            if(isImplicit):
                print("Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations")

            delta=0
            for k in range(nbCells):
                unknown_field[k]  =Un[k]
                if (abs(stat_field[k]-unknown_field[k])>delta):
                    delta=abs(stat_field[k]-unknown_field[k])
                
            unknown_field.setTime(time,it);
            unknown_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName,False);

            print("Ecart au stationnaire exact = ",delta)
            print()
    print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) =
    print( "Variation temporelle relative : ", maxVector[0] )

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint" )
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time )
        if(not with_source and test_bc=="Periodic"):
            print( "Mass loss: ", (total_masse_initial-unknown_field.integral()).norm(), " precision required= ", precision )
            assert (total_mass_initial-unknown_field.integral()).norm()<precision
            if():
        print("------------------------------------------------------------------------------------")

        unknown_field.setTime(time,0);
        unknown_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat");

        #Postprocessing : Extraction of the diagonal data
        if(dim==2):
            diag_data_u=VTK_routines.Extract_field_data_over_line_to_numpyArray(unknown_field,[0,1,0],[1,0,0], resolution)    
        elif(dim==3):
            diag_data_u=VTK_routines.Extract_field_data_over_line_to_numpyArray(unknown_field,[0,0,0],[1,1,1], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat"+'_0.vtu',"unknown",'CELLS',"TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat")
        
        return delta, nbCells, time, it, unknown_field.getNormEuclidean().max(), diag_data_u
    else:
        print("Temps maximum Tmax= ", tmax, " atteint" )
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh, meshName, resolution, meshType, testColor, cfl, test_bc, with_source=False):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    test_name="Resolution of the Transport Equation in dimension " +str( my_mesh.getSpaceDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    test_name_comment="Classical characteristic based scheme"
    test_model="Transport_Equation"
    test_method="Upwind"
    if(with_source):
        test_initial_data="Initial unknown is zero"
    else:
        test_initial_data="Spherical shock"

    print( test_name )
    print("Numerical method : ", test_method )
    isImplicit=True
    if( isImplicit>0):
        print("Time scheme : Implicit" )
        test_desc["Numerical_method_time_discretization"]="Implicit"
    else:
        print("Time scheme : Explicit" )
        test_desc["Numerical_method_time_discretization"]="Explicit"
    print("Initial data : ", test_initial_data )
    print("Boundary conditions : ",test_bc )
    print("Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells" )
    
    # Problem data
    tmax = 10000.
    ntmax = 30000
    output_freq = 10000
    
    velocity=cdmath.Vector(dim,1)

    error, nbCells, t_final, ndt_final, max_unknown, diag_data_u = TransportEquationVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution, test_bc,velocity,isImplicit,with_source)
    end = time.time()

    test_desc["Global_name"]=test_name
    test_desc["Global_comment"]=test_name_comment
    test_desc["PDE_model"]=test_model
    test_desc["PDE_is_stationary"]=False
    test_desc["PDE_search_for_stationary_solution"]=True
    test_desc["Numerical_method_name"]=test_method
    test_desc["Numerical_method_space_discretization"]="Finite volumes"
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
    test_desc["Absolute_error"]=error

    with open('test_TransportEquation'+str(my_mesh.getMeshDimension())+'DUpwind_'+meshName+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return error, nbCells, t_final, ndt_final, max_unknown, diag_data_u, end - start

def solve_file( filename,meshName, resolution,meshType, testColor,cfl,test_bc,with_source=False):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution, meshType, testColor,cfl,test_bc,with_source)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0.,1.,20,0.,1.,20,0)
    cfl=0.5
    solve(M1,"SquareRegularTriangles",100,"Regular triangles","Green",cfl,"Periodic",True)

    M2=cdmath.Mesh(0.,1.,10,0.,1.,10,0.,1.,10,6)
    cfl=1./3
    solve(M2,"CubeRegularTetrahedra",100,"Regular tetrahedra","Green",cfl,"Wall")
