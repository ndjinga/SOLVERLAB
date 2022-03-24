#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Système linéaire issu du système des ondes 2D sans terme source
#               Discrétisation VF du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Schéma implicite (centré, upwind ou pseudo-décalé) sur un maillage général
#               Initialisation par un vortex stationnaire
#               Conditions aux limites parois ou Neumann ou périodique
#		        Solveur linéaire exact ou itératif (GMRES, CGNE)
#		        Résolution d'un seul pas de temps
#================================================================================================================================


from math import sin, cos, pi, sqrt
from numpy import sign
import time
import cdmath
import sys

p0=155e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars, 600K
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_square_vortex(my_mesh):
    print( "Initial data : Square vortex (Constant pressure, divergence free velocity)")
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

def jacobianMatrices(normal, coeff, signun, num_scheme):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)

    if num_scheme.lower() in ['centered','centred','centré']:
        for i in range(dim):
            A[i+1,0]=normal[i]*coeff
            A[0,i+1]=c0*c0*normal[i]*coeff
        
        return A*(1./2)
    elif num_scheme.lower() in ['upwind','amont','décentré']:
        absA=cdmath.Matrix(dim+1,dim+1)
    
        absA[0,0]=c0*coeff
        for i in range(dim):
            A[i+1,0]=      normal[i]*coeff
            A[0,i+1]=c0*c0*normal[i]*coeff
            for j in range(dim):
                absA[i+1,j+1]=c0*normal[i]*normal[j]*coeff
        
        return (A - absA)*(1./2)
    elif num_scheme.lower() in ['staggered','stag','décalé','pstag','pseudodécalé']:
        absA=cdmath.Matrix(dim+1,dim+1)
    
        for i in range(dim):
            A[i+1,0]=normal[i]*coeff
            absA[i+1,0]=-signun*A[i+1,0]
            A[0,i+1]=c0*c0*normal[i]*coeff
            absA[0,i+1]=signun*A[0,i+1]
        
        return (A-absA)*(1./2)
    else:
        print('Numerical scheme requested is ', num_scheme)
        raise ValueError("WaveSystemSpectrum : schemes allowed are upwind, centred or pstag")
        
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt, num_scheme ):
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
            Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure(),signun, num_scheme );

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

def WaveSystemVF( cfl, my_mesh, filename, num_scheme, lin_solver):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    dt = 0.
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    iterGMRESMax=200
    
    #iteration vectors
    Un=cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(filename.find("square")>-1 or filename.find("Square")>-1 or filename.find("cube")>-1 or filename.find("Cube")>-1):
        pressure_field, velocity_field = initial_conditions_square_vortex(my_mesh)
    elif(filename.find("disk")>-1 or filename.find("Disk")>-1 or filename.find("Hexagon")>-1):
        pressure_field, velocity_field = initial_conditions_disk_vortex(my_mesh)
    else:
        print( "Mesh name : ", filename)
        raise ValueError("Mesh name should contain substring square, cube, Hexagon or disk")

    for k in range(nbCells):
        Un[k*(dim+1)+0] =      pressure_field[k]
        Un[k*(dim+1)+1] = rho0*velocity_field[k,0]
        Un[k*(dim+1)+2] = rho0*velocity_field[k,1]
        if(dim==3):
            Un[k*(dim+1)+3] = rho0*velocity_field[k,2]
            
    #sauvegarde de la donnée initiale
    pressure_field.setTime(0,0);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"D"+num_scheme+meshName+"_pressure");
    velocity_field.setTime(0,0);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"D"+num_scheme+meshName+"_velocity");
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt, num_scheme)

    # Add the identity matrix on the diagonal
    for j in range(nbCells*(dim+1)):
        divMat.addValue(j,j,1.)
        
    if lin_solver=='gmres' :
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")
    elif lin_solver=='cgne' :
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "CGNE","ICC")
    elif lin_solver=='lu' :
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","LU")
    else :
        raise ValueError("Linear solver "+lin_solver+" not accepted. Use gmres, cgne or lu")

    # Solving the linear system
    dUn=Un.deepCopy()
    LS.setSndMember(Un)
    start = time.time()
    Un=LS.solve();
    end = time.time()
    cvgceLS=LS.getStatus();
    iterGMRES=LS.getNumberOfIter();
    print( "Linear system converged in ", iterGMRES, " iterations with solver : ", LS.getNameOfMethod()," and préconditioner : ", LS.getNameOfPc() )
    if(not cvgceLS):
        print( "Linear system did not converge ", iterGMRES, " iterations of ", LS.getNameOfMethod())
        raise ValueError("Pas de convergence du système linéaire");
        
    dUn-=Un
        
    maxVector=dUn.maxVector(dim+1)
    print( "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0 )

    for k in range(nbCells):
        pressure_field[k]=Un[k*(dim+1)+0]
        velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
        if(dim>1):
            velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
            if(dim>2):
                velocity_field[k,2]=Un[k*(dim+1)+3]/rho0
        
    pressure_field.setTime(dt,1);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"D"+num_scheme+meshName+"_pressure",False);
    velocity_field.setTime(dt,1);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"D"+num_scheme+meshName+"_velocity",False);

    return iterGMRES, end - start, nbCells
    
def solve(my_mesh,meshName,num_scheme, lin_solver, cfl):
    print( "Resolution of the Wave system in dimension ", my_mesh.getSpaceDimension() )
    print( "Numerical method : implicit ", num_scheme, " CFL ", cfl)
    print( "Linear solver ", lin_solver)
    print( "Mesh name : ",meshName , my_mesh.getNumberOfCells(), " cells" )
    
    return WaveSystemVF( cfl, my_mesh, meshName, num_scheme, lin_solver)

def solve_file( filename,meshName, num_scheme, lin_solver,cfl):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, filename+str(my_mesh.getNumberOfCells()),num_scheme, lin_solver, cfl)
    

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        filename=sys.argv[1]
        my_mesh = cdmath.Mesh(filename)
        if len(sys.argv) >2 :
            num_scheme = sys.argv[2]
            if len(sys.argv) >3 :
                lin_solver = sys.argv[3]
            else:
               raise ValueError("WaveSystemCentered.py expects a linear solver name")
        else:
           raise ValueError("WaveSystemCentered.py expects a numerical method name")
        cfl = 1000000.
        solve(my_mesh,filename, num_scheme, lin_solver, cfl)
    else :
        raise ValueError("WaveSystemCentered.py expects a mesh file name")
