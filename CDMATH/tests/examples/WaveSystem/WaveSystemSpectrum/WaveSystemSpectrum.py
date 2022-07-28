#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Calcul VF du spectre du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Utilisation du schéma colocalisé centré ou amont sur un maillage général
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
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax, num_scheme):
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
            Am=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(),signun, num_scheme);

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

def WaveSystemSpectrum( cfl, my_mesh, filename, num_scheme):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax, num_scheme)

    # Add the identity matrix on the diagonal
    divMat.diagonalShift(1/dt)#only after  filling all coefficients
    divMat.viewNonZeroStructure( 0, "FiniteVolumesMatrixOn"+meshName+"_WaveSystem"+num_scheme)
    divMat.saveToFile( "FiniteVolumesMatrixOn"+meshName+"_WaveSystem"+num_scheme,  True)
    
    #Plot the spectrum of the linear system matrix
    X,Y=divMat.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_WaveSystem"+num_scheme)
    plt.xlim((min(X)-50)*1.1, (max(X)+50)*1.1)
    plt.ylim((min(Y)-10)*1.1, (max(Y)+10)*1.1)
    plt.title('Eigenvalues of the '+num_scheme+ ' finite volume method \n for the wave system')
    plt.scatter(X,Y, label=num_scheme+' scheme')#
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(block=False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_WaveSystem.png")

    divMatSquare = divMat*divMat.transpose()
    X2,Y2=divMatSquare.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_WaveSystem"+num_scheme+"_symmetrised")
    plt.xlim((min(X2)-50)*1.1, (max(X2)+50)*1.1)
    plt.ylim((min(Y2)-10)*1.1, (max(Y2)+10)*1.1)
    plt.title('Eigenvalues of tAA for the '+num_scheme+ ' finite volume method \n for the wave system')
    plt.scatter(X2,Y2, label=num_scheme+' scheme'+" symmetrised")#
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(block=False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_WaveSystem_symmetrised.png")

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
    print( "Spectrum of the Wave system in dimension ", my_mesh.getSpaceDimension() )
    if   num_scheme.lower() in ['upwind','amont','décentré']:
        print( "Numerical method : ", "Upwind" )
    elif num_scheme.lower() in ['centered','centred','centré']:    
        print( "Numerical method : ", "Centered" )
    elif num_scheme.lower() in ['staggered','stag','décalé','pstag','pseudodécalé']:    
        print( "Numerical method : ", "pseudo-staggered" )
    else:
        print('Numerical scheme requested is ', num_scheme)
        raise ValueError("WaveSystemSpectrum : schemes allowed are upwind, centred or pstag")
    print( "Wall boundary conditions" )
    print( "Mesh name : ",meshName , my_mesh.getNumberOfCells(), " cells" )
    
    # Problem data
    cfl = 100000./my_mesh.getSpaceDimension()

    WaveSystemSpectrum( cfl, my_mesh, meshName, num_scheme)

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
        raise ValueError("WaveSystemSpectrum.py expects a mesh file name")
