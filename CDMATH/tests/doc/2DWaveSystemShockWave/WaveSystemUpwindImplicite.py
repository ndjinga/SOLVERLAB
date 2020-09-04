#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import time, json
import cdmath
import PV_routines
import VTK_routines

p0=155e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_square_shock(my_mesh):
    print "Spherical wave in a square : initial data"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    rayon = 0.15
    xcentre = 0.5
    ycentre = 0.5
    zcentre = 0.5


    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()

        velocity_field[i,0] = 0
        velocity_field[i,1] = 0
        velocity_field[i,2] = 0

        valX = (x - xcentre) * (x - xcentre)
        valY = (y - ycentre) * (y - ycentre)

        if(dim==2):
            val =  sqrt(valX + valY)
        if(dim==3):
            z = my_mesh.getCell(i).z()
            valZ = (z - zcentre) * (z - zcentre)
            val =  sqrt(valX + valY + valZ)

        if val < rayon:
            pressure_field[i] = p0
            pass
        else:
            pressure_field[i] = p0/2
            pass
        pass

    return pressure_field, velocity_field


def jacobianMatrices(normal, coeff,scaling):
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
                    raise ValueError("computeFluxes: problem with mesh, unknown cell number")
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                indexFP = my_mesh.getIndexFacePeriodic(indexFace, my_mesh.getName()== "squareWithBrickWall", my_mesh.getName()== "squareWithHexagons")
                Fp = my_mesh.getFace(indexFP)
                cellAutre = Fp.getCellsId()[0]
                
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                
    return implMat

def WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,scaling,test_bc,with_source=False):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    isImplicit=True
    
    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(meshName.find("square")>-1 or meshName.find("Square")>-1 or meshName.find("cube")>-1 or meshName.find("Cube")>-1):
        pressure_field, velocity_field = initial_conditions_square_shock(my_mesh)
        stat_pressure, stat_velocity   = initial_conditions_square_shock(my_mesh)
    elif(meshName.find("disk")>-1 or meshName.find("Disk")>-1):
        pressure_field, velocity_field = initial_conditions_disk_vortex(my_mesh)
        stat_pressure, stat_velocity   = initial_conditions_disk_vortex(my_mesh)
    else:
        print "Mesh name : ", meshName
        raise ValueError("Mesh name should contain substring square, cube or disk")

    for k in range(nbCells):
        Un[k*(dim+1)+0] =     pressure_field[k]
        Un[k*(dim+1)+1] =rho0*velocity_field[k,0]
        Un[k*(dim+1)+2] =rho0*velocity_field[k,1]
        if(dim==3):
            Un[k*(dim+1)+3] =rho0*velocity_field[k,2]

    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity");
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_initial")
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_initial")

    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling,test_bc)
    if( isImplicit):
        #Adding the identity matrix on the diagonal
        divMat.diagonalShift(1)#only after  filling all coefficients
        #for j in range(nbCells*(dim+1)):
        #    divMat.addValue(j,j,1)
        
        iterGMRESMax=50
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")

        LS.setComputeConditionNumber()
        
    print("Starting computation of the linear wave system with an Upwind scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            LS.setSndMember(Un)
            Un=LS.solve();
            if(not LS.getStatus()):
                print "Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations"
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

        else:
            dUn=divMat*Un-S*dt
            Un-=dUn
        
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            if(isImplicit):
                print "Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations"

            for k in range(nbCells):
                pressure_field[k]  =Un[k*(dim+1)+0]
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity",False);

            print
    print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
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
        
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh, meshName, resolution, scaling, meshType, testColor, cfl, test_bc, with_source=False):
    print "Resolution of the Wave system in dimension " +str( my_mesh.getSpaceDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    print "Initial data : ", "Riemann problem"
    print "Boundary conditions : ", "Neumann"
    print "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells"
    
    # Problem data
    tmax = 10000.
    ntmax = 500
    output_freq = 1

    WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution, scaling,test_bc,with_source)

    return 

def solve_file( filename,meshName, resolution, scaling, meshType, testColor,cfl,test_bc,with_source=False):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution, scaling, meshType, testColor,cfl,test_bc,with_source)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0.,1.,50,0.,1.,50)
    cfl=0.5
    scaling=0
    solve(M1,"SquareRegularSquares",100,scaling,"RegularSquares","Green",cfl,"Periodic",False)
