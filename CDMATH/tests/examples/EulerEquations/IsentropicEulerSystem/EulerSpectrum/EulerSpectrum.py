#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Calcul VF du spectre d'Euler isotherme
#                \partial_t rho + \div q = 0
#                \partial_t q   + \div q\otimes q/rho  +  \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Utilisation du schéma pstag, centré ou amont sur un maillage général
#              
#================================================================================================================================


from math import sqrt
from numpy import sign
from matplotlib import pyplot as plt
import cdmath
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

def jacobianMatrices(normal, coeff, signun, num_scheme, rho_l,q_l,rho_r,q_r):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)
    absA=cdmath.Matrix(dim+1,dim+1)

    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r )
        raise ValueError("Negative density")

    u_l = cdmath.Vector(dim); 
    u_r = cdmath.Vector(dim); 
    u = cdmath.Vector(dim);
    
    for i in range(dim):
        u_l[i]=q_l[i]/rho_l
        u_r[i]=q_r[i]/rho_r
        u[i] = (u_l[i]*sqrt(rho_l)+u_r[i]*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));

    un=u*normal;

    #Compute jacobian matrix
    for i in range(dim):
        A[0,i+1] =        normal[i]
        A[i+1,0] = (c0*c0*normal[i] - un*u[i])
        A[i+1,i+1] = un
        for j in range(dim):
            A[i+1,1+j]+=normal[j]*u[i]

    if num_scheme.lower() in ['centered','centred','centré']:
        
        return A*(coeff/2), un
        
    elif num_scheme.lower() in ['upwind','amont','décentré']:
        absA=cdmath.Matrix(dim+1,dim+1)
    
        tangent=cdmath.Vector(dim);
        if dim ==2 :
            tangent[0]= normal[1];
            tangent[1]=-normal[0];
        elif dim >2:
            raise ValueError("Implementation not ready for dimension 3")
            
        #Il faudrait swiger l'opérateur de produit vectoriel : utiliser .tensProd
        #subMatrix=(abs(un+c0)*((u-c0*normal)^normal)-abs(un-c0)*((u-c0*normal)^normal))/(2*c0)+abs(un)*(tangent^tangent)
        
        absA[0,0]=(abs(un-c0)*(un+c0)+abs(un+c0)*(c0-un))/(2*c0);
        for i in range(dim):
            absA[i+1,0]=(abs(un-c0)*(un+c0)*(u[i]-c0*normal[i])-abs(un+c0)*(un-c0)*(u[i]+c0*normal[i]))/(2*c0)-abs(un)*(u*tangent)*tangent[i];
            absA[0,i+1]=(abs(un+c0)-abs(un-c0))/(2*c0)*normal[i];
            #Il faudrait swiger l'opérateur de produit vectoriel
            #for j in range(dim):
            #    absA[i+1,j+1]=subMatrix[i,j]

        if dim ==1 :
            absA[1,1]=(abs(un+c0)*((u[0]+c0*normal[0])*normal[0])-abs(un-c0)*((u[0]-c0*normal[0])*normal[0]))/(2*c0)
        elif dim ==2 :
            absA[1,1]=(abs(un+c0)*((u[0]+c0*normal[0])*normal[0])-abs(un-c0)*((u[0]-c0*normal[0])*normal[0]))/(2*c0)+abs(un)*(tangent[0]*tangent[0]);
            absA[1,2]=(abs(un+c0)*((u[0]+c0*normal[0])*normal[1])-abs(un-c0)*((u[0]-c0*normal[0])*normal[1]))/(2*c0)+abs(un)*(tangent[0]*tangent[1]);
            absA[2,1]=(abs(un+c0)*((u[1]+c0*normal[1])*normal[0])-abs(un-c0)*((u[1]-c0*normal[1])*normal[0]))/(2*c0)+abs(un)*(tangent[1]*tangent[0]);
            absA[2,2]+=(abs(un+c0)*((u[1]+c0*normal[1])*normal[1])-abs(un-c0)*((u[1]-c0*normal[1])*normal[1]))/(2*c0)+abs(un)*(tangent[1]*tangent[1]);

        return (A - absA)*(coeff/2), un
    elif num_scheme.lower() in ['staggered','stag','décalé','pstag','pseudodécalé']:
        absA=cdmath.Matrix(dim+1,dim+1)
    
        for i in range(dim):
            absA[i+1,0]=-signun*A[i+1,0]
            absA[0,i+1]=signun*A[0,i+1]
        
        return (A-absA)*(coeff/2), un
    else:
        print('Numerical scheme requested is ', num_scheme)
        raise ValueError("IsentropicEulerSystemSpectrum : schemes allowed are upwind, centred or pstag")
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax, num_scheme, Un):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1
    normal=cdmath.Vector(dim)
    maxAbsEigVa = 0

    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    #Pour les conditions limites
    q_l=cdmath.Vector(dim);
    q_r=cdmath.Vector(dim);
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

            cellAutre =-1#ghost cell
            if ( not Fk.isBorder()) :
                # hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if (Fk.getCellsId()[0] == j) :
                    # hypothese verifiée 
                    cellAutre = Fk.getCellsId()[1];
                elif(Fk.getCellsId()[1] == j) :
                    # hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                else :
                    raise ValueError("computeDivergenceMatrix: problem with mesh, unknown cel number")
                    
                q_l[0]=Un[j*nbComp+1]
                q_l[1]=Un[j*nbComp+2]
                q_r[0]=Un[cellAutre*nbComp+1]
                q_r[1]=Un[cellAutre*nbComp+2]
                Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(), signun, num_scheme,Un[j*nbComp], q_l,Un[cellAutre*nbComp],q_r);
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                if( Fk.getGroupName() != "Periodic" and Fk.getGroupName() != "Neumann"):#Wall boundary condition unless Periodic/Neumann specified explicitly
                    q_l[0]=Un[j*nbComp+1]
                    q_l[1]=Un[j*nbComp+2]
                    q_r=q_l-normal*2*(q_l*normal)
                    Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(), signun, num_scheme,Un[j*nbComp], q_l,Un[j*nbComp],q_r);

                    v=cdmath.Vector(dim+1)
                    for i in range(dim) :
                        v[i+1]=normal[i]
                    idMoinsJacCL=v.tensProduct(v)*2
                    
                    implMat.addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL)
                    
                elif( Fk.getGroupName() == "Periodic"):#Periodic boundary condition
                    indexFP=my_mesh.getIndexFacePeriodic(indexFace)
                    Fp = my_mesh.getFace(indexFP)
                    cellAutre = Fp.getCellsId()[0]
                    
                    q_l[0]=Un[j*nbComp+1]
                    q_l[1]=Un[j*nbComp+2]
                    q_r[0]=Un[cellAutre*nbComp+1]
                    q_r[1]=Un[cellAutre*nbComp+2]
                    Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(), signun, num_scheme, Un[j*nbComp],q_l,Un[cellAutre*nbComp],q_r);
            
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif(Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print( Fk.getGroupName() )
                    raise ValueError("computeDivergenceMatrix: Unknown boundary condition name");
                
            maxAbsEigVa = max(maxAbsEigVa,abs(un+c0),abs(un-c0));

    return implMat, maxAbsEigVa

def IsentropicEulerSystemSpectrum( cfl, my_mesh, filename, num_scheme):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    # Initial conditions #
    print("Construction of the initial condition …")
    if(dim==1 or filename.find("square")>-1 or filename.find("Square")>-1 or filename.find("cube")>-1 or filename.find("Cube")>-1):
        pressure_field, velocity_field = initial_conditions_shock(my_mesh,False)
    elif(filename.find("disk")>-1 or filename.find("Disk")>-1 or filename.find("Hexagon")>-1 or filename.find("HEXAON")>-1):
        pressure_field, velocity_field = initial_conditions_shock(my_mesh,True)
    else:
        print( "Mesh name : ", filename )
        raise ValueError("Mesh name should contain substring square, cube, disk or hexagon")

    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    for k in range(nbCells):
        Un[k*(dim+1)+0] =                pressure_field[k]/(c0*c0)
        Un[k*(dim+1)+1] =Un[k*(dim+1)+0]*velocity_field[k,0]
        if(dim>=2):
            Un[k*(dim+1)+2] = Un[k*(dim+1)+0]*velocity_field[k,1]
            if(dim==3):
                Un[k*(dim+1)+3] = Un[k*(dim+1)+0]*velocity_field[k,2]
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
    dx_min=my_mesh.minRatioVolSurf()

    divMat, maxAbsEigVa=computeDivergenceMatrix(my_mesh,nbVoisinsMax, num_scheme, Un)

    dt = cfl * dx_min /  maxAbsEigVa

    # Add the identity matrix on the diagonal
    divMat.diagonalShift(1/dt)#only after  filling all coefficients
    divMat.viewNonZeroStructure( 0, "FiniteVolumesMatrixOn"+meshName+"_IsentropicEulerSystem"+num_scheme)
    divMat.saveToFile( "FiniteVolumesMatrixOn"+meshName+"_IsentropicEulerSystem"+num_scheme,  True)
    
    #Plot the spectrum of the linear system matrix
    X,Y=divMat.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_IsentropicEulerSystem"+num_scheme)
    plt.xlim((min(X)-50)*1.1, (max(X)+50)*1.1)
    plt.ylim((min(Y)-10)*1.1, (max(Y)+10)*1.1)
    plt.scatter(X,Y, label=num_scheme+' scheme')#
    plt.title('Eigenvalues of the '+num_scheme+ ' finite volume method \n for the Isentropic Euler system')
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(block=False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_IsentropicEulerSystem.png")

    plt.clf()
    #Plot the spectrum of the symmetrisation of the linear system matrix
    divMatSquare = divMat*divMat.transpose()
    X2,Y2=divMatSquare.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_IsentropicEulerSystem"+num_scheme+"_symmetrised")
    plt.xlim((min(X2)-50)*1.1, (max(X2)+50)*1.1)
    plt.ylim((min(Y2)-10)*1.1, (max(Y2)+10)*1.1)
    plt.scatter(X2,Y2, label=num_scheme+' scheme'+" symmetrised")#
    plt.title('Eigenvalues of tAA for the '+num_scheme+ ' finite volume method \n for the Isentropic Euler system')
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(block=False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_IsentropicEulerSystem_symmetrised.png")

    assert min(map(abs,Y2))<precision
    sigma_min = min(X2)
    sigma_max = max(X2)
    
    if sigma_min >0.:#avoid division by zero
        print("Minimum singular value   A = ", sqrt(sigma_min), ", maximum singular value   A = ", format(sqrt(sigma_max),'.1E'), ", conditionnement   A = ", format(sqrt(sigma_max/sigma_min),'.1E') )
        print("Minimum singular value tAA = ",      sigma_min,  ", maximum singular value tAA = ", format(     sigma_max ,'.1E'), ", conditionnement tAA = ", format(     sigma_max/sigma_min ,'.1E') )
    else:
        print("Minimum singular value   A = ", 0, ", maximum singular value   A = ", format(sqrt(sigma_max),'.1E'), ", conditionnement   A = infinity" )
        print("Minimum singular value tAA = ", 0, ", maximum singular value tAA = ", format(     sigma_max ,'.1E'), ", conditionnement tAA = infinity" )

def solveSpectrum(my_mesh,meshName, num_scheme):
    print( "Spectrum of the IsentropicEuler system in dimension ", my_mesh.getSpaceDimension() )
    if   num_scheme.lower() in ['upwind','amont','décentré']:
        print( "Numerical method : ", "Upwind finite volume" )
    elif num_scheme.lower() in ['centered','centred','centré']:    
        print( "Numerical method : ", "Centered finite volume" )
    elif num_scheme.lower() in ['staggered','stag','décalé','pstag','pseudodécalé']:    
        print( "Numerical method : ", "pseudo-staggered finite volume" )
    else:
        print('Numerical scheme requested is ', num_scheme)
        raise ValueError("IsentropicEulerSystemSpectrum : schemes allowed are upwind, centred or pstag")
    print( "Wall boundary conditions" )
    print( "Mesh name : ",meshName )
    print(  my_mesh.getNumberOfCells(), " cells" )
    
    # Problem data
    cfl = 100000./my_mesh.getSpaceDimension()

    IsentropicEulerSystemSpectrum( cfl, my_mesh, meshName, num_scheme)

def solve_file_spectrum( filename,meshName, num_scheme):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, filename+str(my_mesh.getNumberOfCells()), num_scheme)
    

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        filename=sys.argv[1]
        if len(sys.argv) >2 :
            num_scheme = sys.argv[2]
        else:
           num_scheme = 'Upwind'
        my_mesh = cdmath.Mesh(filename)
        solveSpectrum(my_mesh,filename, num_scheme)
    else :
        raise ValueError("IsentropicEulerSystemSpectrum.py expects a mesh file name")
