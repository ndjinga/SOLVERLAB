#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Calcul VF du spectre de l'opérateur de transport 2D \ddt u + \vec a\cdot\nabla u = 0 
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
# Description : Utilisation de la méthode des volumes finis avec champs discrétisés aux cellules d'un maillage quelconque
#		Schéma amont ou centré
#================================================================================================================================

from math import sin, cos, pi, sqrt
import sys
from matplotlib import pyplot as plt
import cdmath

precision=1e-5

def upwinding_coeff(normal, coeff, velocity, is_upwind):
    dim=normal.size()
    
    if(not is_upwind):
        return (velocity*normal)/2
    else:
        if(velocity*normal<0.):
            return velocity*normal*coeff
        else:
            return 0.
    
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,test_bc,velocity, is_upwind):
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

            Am=upwinding_coeff( normal,Fk.getMeasure()/Cj.getMeasure(),velocity, is_upwind);

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
                if( test_bc=="Periodic" and Fk.getGroupName() != "Neumann"):#Periodic boundary condition unless Wall/Neumann specified explicitly
                    indexFP = my_mesh.getIndexFacePeriodic(indexFace, my_mesh.getName()== "squareWithBrickWall", my_mesh.getName()== "squareWithHexagons")
                    Fp = my_mesh.getFace(indexFP)
                    cellAutre = Fp.getCellsId()[0]
                    
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif(test_bc!="Neumann" and Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print( Fk.getGroupName() )
                    raise ValueError("computeFluxes: Unknown boundary condition name");
                
    return implMat


def solveSpectrum(my_mesh, meshName, meshType, cfl, test_bc, is_upwind):
    print( "Spectrum of the Transport Equation in dimension ", my_mesh.getMeshDimension() )
    if( is_upwind ):
        print( "Numerical method : ", "Upwind" )
    else:
        print( "Numerical method : ", "Centered" )
    print( "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells" )
    
    
    dim=my_mesh.getMeshDimension()
    velocity=cdmath.Vector(dim)
    for i in range(dim) :
        velocity[i] = 1

    nbCells = my_mesh.getNumberOfCells()
    
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / velocity.norm()

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,test_bc,velocity, is_upwind)
    #Adding the identity matrix on the diagonal
    divMat.diagonalShift(1/dt)#only after  filling all coefficients
    
    if is_upwind:
        num_scheme='Upwind'
    else:
        num_scheme='Centred'
    divMat.viewNonZeroStructure(0, "FiniteVolumesMatrixOn"+meshName+"_TransportEquation"+num_scheme)
    divMat.saveToFile( "FiniteVolumesMatrixOn"+meshName+"_TransportEquation"+num_scheme, True)
    X,Y=divMat.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_TransportEquation"+num_scheme)

    plt.xlim((min(X)-50)*1.1, (max(X)+50)*1.1)
    plt.ylim((min(Y)-10)*1.1, (max(Y)+10)*1.1)
    plt.title('Eigenvalues of the '+num_scheme+ ' finite volume method \n for the transport equation')
    #Plot the spectrum of the linear system matrix
    plt.scatter(X,Y, label=num_scheme+' scheme')
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_TransportEquation.png")
    
    divMatSquare = divMat*divMat.transpose()
    X2,Y2=divMatSquare.plotEigenvalues("FiniteVolumesEigenvaluesOn"+meshName+"_TransportEquation"+num_scheme+"_symmetrised")
    plt.scatter(X2,Y2, label=num_scheme+' scheme'+" symmetrised")#
    plt.xlabel("Real part")
    plt.ylabel("Imaginary part")
    plt.axvline(x=0, color='black')
    plt.axhline(y=0, color='black')
    plt.legend()
    plt.show(False)
    plt.savefig("FiniteVolumesEigenvaluesOn"+meshName+"_"+num_scheme+"Scheme"+"_TransportEquation_symmetrised.png")
    
    assert abs(min(Y2))<precision
    sigma_min = min(X2)
    sigma_max = max(X2)
    print("Minimum singular value   A = ", sqrt(sigma_min), ", maximum singular value   A = ", format(sqrt(sigma_max),'.1E'), ", conditionnement   A = ", format(sqrt(sigma_max/sigma_min),'.1E') )
    print("Minimum singular value tAA = ",      sigma_min,  ", maximum singular value tAA = ", format(     sigma_max ,'.1E'), ", conditionnement tAA = ", format(     sigma_max/sigma_min ,'.1E') )

def solve_file( filename,meshName, meshType, cfl, test_bc, is_upwind):
    my_mesh = cdmath.Mesh(filename+".med")

    return solveSpectrum(my_mesh, meshName+str(my_mesh.getNumberOfCells()), meshType, cfl, test_bc, is_upwind)
    

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        is_upwind = sys.argv[1].lower() in ['false', '0', 'f', 'n', 'no']
    else :
        is_upwind = True

    M1=cdmath.Mesh(0.,1.,12,0.,1.,12)
    cfl=1000000
    solveSpectrum(M1,"SquareRegularTriangles","Regular triangles",cfl,"Periodic", is_upwind)
	
