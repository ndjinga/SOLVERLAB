#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Propagation d'une onde de choc sphérique
#               Utilisation du schéma centré (implicite) sur un maillage général
#               Initialisation par une surpression sphérique
#               Conditions aux limites parois
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sqrt
from numpy import sign
import cdmath
import PV_routines
import VTK_routines
import sys

p0=155e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars, 600K
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_shock(my_mesh, isCircle):
    print( "Initial data : Spherical wave" )
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    rayon = 0.15
    if(not isCircle):
        xcentre = 0.5
        ycentre = 0.5
        zcentre = 0.5
    else:
        xcentre = 0.
        ycentre = 0.
        zcentre = 0.

    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        velocity_field[i,0] = 0
        velocity_field[i,1] = 0
        velocity_field[i,2] = 0

        x = my_mesh.getCell(i).x()
        valX = (x - xcentre) * (x - xcentre)

        if(dim==1):
            val =  sqrt(valX)
        if(dim==2):
            y = my_mesh.getCell(i).y()
            valY = (y - ycentre) * (y - ycentre)
            val =  sqrt(valX + valY)
        if(dim==3):
            y = my_mesh.getCell(i).y()
            z = my_mesh.getCell(i).z()
            valY = (y - ycentre) * (y - ycentre)
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

def jacobianMatrices(normal, coeff, signun):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)

    for i in range(dim):
        A[i+1,0]=normal[i]*coeff
        A[0,i+1]=c0*c0*normal[i]*coeff
    
    return A*(1./2)
    
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1
    normal=cdmath.Vector(dim)

    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    idMoinsJacCL=cdmath.Matrix(nbComp)

    v0=cdmath.Vector(dim)
    for i in range(dim) :
        v0[i] = 1.

    for j in range(nbCells):#On parcourt les cellules
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

            signun=sign(normal*v0)
            Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure(),signun);

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
                if( Fk.getGroupName() != "Periodic" and Fk.getGroupName() != "Neumann"):#Wall boundary condition unless Periodic/Neumann specified explicitly
                    v=cdmath.Vector(dim+1)
                    for i in range(dim) :
                        v[i+1]=normal[i]
                    idMoinsJacCL=v.tensProduct(v)*2
                    
                    implMat.addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL)
                    
                elif( Fk.getGroupName() == "Periodic"):#Periodic boundary condition
                    indexFP=my_mesh.getIndexFacePeriodic(indexFace)
                    Fp = my_mesh.getFace(indexFP)
                    cellAutre = Fp.getCellsId()[0]
                    
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif(Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print( Fk.getGroupName() )
                    raise ValueError("computeFluxes: Unknown boundary condition name");
                
    return implMat

def WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    iterGMRESMax=50
    
    #iteration vectors
    Un=cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(dim==1 or filename.find("square")>-1 or filename.find("Square")>-1 or filename.find("cube")>-1 or filename.find("Cube")>-1):
        pressure_field, velocity_field = initial_conditions_shock(my_mesh,False)
    elif(filename.find("disk")>-1 or filename.find("Disk")>-1 or filename.find("Hexagon")>-1):
        pressure_field, velocity_field = initial_conditions_shock(my_mesh,True)
    else:
        print( "Mesh name : ", filename )
        raise ValueError("Mesh name should contain substring square, cube, Hexagon or disk")

    for k in range(nbCells):
        Un[k*(dim+1)+0] =      pressure_field[k]
        Un[k*(dim+1)+1] = rho0*velocity_field[k,0]
        if(dim>=2):
            Un[k + 2*nbCells] = rho0*velocity_field[k,1] # value on the bottom face
            if(dim==3):
                Un[k + 3*nbCells] = rho0*velocity_field[k,2]
            
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity");
    #Postprocessing : save 2D picture
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_initial")
    PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity_initial")
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt)

    # Add the identity matrix on the diagonal
    for j in range(nbCells*(dim+1)):
        divMat.addValue(j,j,1.)
    LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","LU")

    print("Starting computation of the linear wave system with a centered scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        dUn=Un.deepCopy()
        LS.setSndMember(Un)
        Un=LS.solve();
        cvgceLS=LS.getStatus();
        iterGMRES=LS.getNumberOfIter();
        if(not cvgceLS):
            print( "Linear system did not converge ", iterGMRES, " GMRES iterations" )
            raise ValueError("Pas de convergence du système linéaire");
        dUn-=Un
        
        maxVector=dUn.maxVector(dim+1)
        isStationary= maxVector[0]/p0<precision and maxVector[1]/rho0<precision and maxVector[2]/rho0<precision;
        if(dim==3):
            isStationary=isStationary and maxVector[3]/rho0<precision
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
            print( "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0 )
            print( "Linear system converged in ", iterGMRES, " GMRES iterations" )

            for k in range(nbCells):
                pressure_field[k]=Un[k*(dim+1)+0]
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0
                
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity",False);

    print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
    print( "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0 )
    print()

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint" )
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time )
        print( "------------------------------------------------------------------------------------" )

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity_Stat");

        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DCentered"+meshName+"_velocity_Stat")
        
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint" )


def solve(my_mesh,meshName,resolution):
    print( "Resolution of the Wave system in dimension ", my_mesh.getSpaceDimension() )
    print( "Numerical method : implicit centered")
    print( "Initial data : spherical wave")
    print( "Wall boundary conditions")
    print( "Mesh name : ",meshName , my_mesh.getNumberOfCells(), " cells")
    
    # Problem data
    tmax = 1000.
    ntmax = 100
    cfl = 1./my_mesh.getSpaceDimension()
    output_freq = 1

    WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution)

def solve_file( filename,meshName, resolution):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, filename+str(my_mesh.getNumberOfCells()),resolution)
    

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        filename=sys.argv[1]
        my_mesh = cdmath.Mesh(filename)
        solve(my_mesh,filename,100)
    else :
        raise ValueError("WaveSystemUpwind.py expects a mesh file name")
