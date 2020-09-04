#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Test de préservation d'un état stationnaire
#               Utilisation du schéma MAC sur un maillage cartésien
#               Initialisation par un vortex stationnaire
#               Conditions aux limites périodiques
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sin, cos, pi, sqrt
import cdmath
import PV_routines
import VTK_routines
import sys

p0=155e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars, 600K
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_square_vortex(my_mesh):
    print( "Initial data : Square vortex (Constant pressure, divergence free velocity)" )
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    pressure_field = cdmath.Field("Pressure", cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity", cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()
        z = my_mesh.getCell(i).z()

        pressure_field[i] = p0
        if(dim==1):
            velocity_field[i,0] = 1
            velocity_field[i,1] = 0
            velocity_field[i,2] = 0
        elif(dim==2):
            velocity_field[i,0] =  sin(pi*x)*cos(pi*y)
            velocity_field[i,1] = -sin(pi*y)*cos(pi*x)
            velocity_field[i,2] = 0
        elif(dim==3):
            velocity_field[i,0] =    sin(pi*x)*cos(pi*y)*cos(pi*z)
            velocity_field[i,1] =    sin(pi*y)*cos(pi*x)*cos(pi*z)
            velocity_field[i,2] = -2*sin(pi*z)*cos(pi*x)*cos(pi*y)
        
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

    return implMat

def WaveSystemStaggered(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    scaling=0
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    iterGMRESMax=50
    
    #iteration vectors
    Un=cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_square_vortex(my_mesh)

    for k in range(nbCells):
        Un[k + 0*nbCells] =      pressure_field[k]
        Un[k + 1*nbCells] = rho0*velocity_field[k,0] # value on the left face
        Un[k + 2*nbCells] = rho0*velocity_field[k,1] # value on the bottom face
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
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt,scaling)

    #Add the identity matrix on the diagonal
    for j in range(nbCells*(dim+1)):
        divMat.addValue(j,j,1)

    LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","LU")

    print("Starting computation of the linear wave system with an pseudo staggered scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        dUn=Un.deepCopy()
        LS.setSndMember(Un)
        Un=LS.solve();
        cvgceLS=LS.getStatus();
        iterGMRES=LS.getNumberOfIter();
        if(not cvgceLS):
            print( "Linear system did not converge ", iterGMRES, " GMRES iterations")
            raise ValueError("Pas de convergence du système linéaire");
        dUn-=Un
        
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
            print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
            print( "Variation temporelle relative : pressure ", max_dp/p0 ,", velocity ", max_dq/rho0)
            print( "Linear system converged in ", iterGMRES, " GMRES iterations")

            for k in range(nbCells):
                pressure_field[k]=Un[k]
                velocity_field[k,0]=Un[k+1*nbCells]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k+2*nbCells]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k+3*nbCells]/rho0
                
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity",False);

    print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
    print( "Variation temporelle relative : pressure ", max_dp/p0 ,", velocity ", max_dq/rho0 )
    print()

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print( "------------------------------------------------------------------------------------")

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat");

        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat")
        
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint" )


def solve(my_mesh,meshName,resolution):
    print( "Resolution of the Wave system in dimension ", my_mesh.getSpaceDimension())
    print( "Numerical method : staggered scheme")
    print( "Initial data : stationary solution (constant pressure, divergence free velocity)")
    print( "Periodic boundary conditions")
    print( "Mesh name : ",meshName , my_mesh.getNumberOfCells(), " cells")
    
    # Problem data
    tmax = 1000.
    ntmax = 100
    cfl = 1./my_mesh.getSpaceDimension()
    output_freq = 100

    WaveSystemStaggered(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution)

def solve_file( filename,meshName, resolution):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution)
    

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        my_mesh = cdmath.Mesh(sys.argv[1])
        solve(my_mesh,my_mesh.getName(),100)
    else :
        nx=50
        my_mesh = cdmath.Mesh(0,1,nx,0,1,nx)
        solve(my_mesh,my_mesh.getName(),100)
