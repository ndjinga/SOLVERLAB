#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Propagation d'une onde de choc droite
#               Utilisation du schéma upwind explicite ou implicite sur un maillage général
#               Initialisation par une discontinuité verticale
#               Conditions aux limites de Neumann
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

def initial_conditions_RiemannProblem(my_mesh):
    print( "Initial data : Riemann problem" )
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    xcentre = 0

    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()

        velocity_field[i,0] = 0
        velocity_field[i,1] = 0
        velocity_field[i,2] = 0

        if x < xcentre:
            pressure_field[i] = p0
            pass
        else:
            pressure_field[i] = p0/2
            pass
        pass

    return pressure_field, velocity_field

def jacobianMatrices(normal,coeff):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)
    absA=cdmath.Matrix(dim+1,dim+1)

    absA[0,0]=c0*coeff
    for i in range(dim):
        A[i+1,0]=      normal[i]*coeff
        A[0,i+1]=c0*c0*normal[i]*coeff
        for j in range(dim):
            absA[i+1,j+1]=c0*normal[i]*normal[j]*coeff
    
    return (A - absA)*(1./2)
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt):
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

            Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure());

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
                
    return implMat

def WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution, isImplicit):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
    SumFluxes = cdmath.Field("Fluxes", cdmath.CELLS, my_mesh, dim+1)

    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_RiemannProblem(my_mesh)

    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    for k in range(nbCells):
        Un[k*(dim+1)+0] =     pressure_field[k]
        Un[k*(dim+1)+1] =rho0*velocity_field[k,0]
        Un[k*(dim+1)+2] =rho0*velocity_field[k,1]
        if(dim==3):
            Un[k*(dim+1)+3] =rho0*velocity_field[k,2]

    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity");

    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0
    
    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt)
    if( isImplicit):
        #Adding the identity matrix on the diagonal
        divMat.diagonalShift(1)#only after  filling all coefficients
        
        iterGMRESMax=50
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")

        LS.setComputeConditionNumber()
        
    print("Starting computation of the linear wave system with an explicit UPWIND scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            LS.setSndMember(Un)
            Un=LS.solve();
            if(not LS.getStatus()):
                print( "Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations" )
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

        else:
            dUn=divMat*Un
            Un-=dUn
        
        time=time+dt;
        it=it+1;
 
         #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )

            for k in range(nbCells):
                pressure_field[k]  =Un[k*(dim+1)+0]
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity",False);

    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    print()

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time)

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
        diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat")
        
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint")


def solve(my_mesh,filename,resolution, isImplicit):
    print( "Resolution of the Wave system in dimension ", my_mesh.getSpaceDimension() )
    print( "Numerical method : upwind")
    print( "Initial data : single straight discontinuity (Riemann problem)")
    print( "Neumann boundary conditions")
    print( "Mesh name : ",filename , my_mesh.getNumberOfCells(), " cells")

    # Problem data
    tmax = 1.
    ntmax = 50
    cfl = 1./my_mesh.getSpaceDimension()
    output_freq = 1

    WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution, isImplicit)

def solve_file( filename,resolution, isImplicit):
    my_mesh = cdmath.Mesh(filename+".med")
    solve(my_mesh, filename,resolution, isImplicit)
    
if __name__ == """__main__""":
    if len(sys.argv) >2 :
        filename=sys.argv[1]
        isImplicit=bool(int(sys.argv[2]))
        my_mesh = cdmath.Mesh(filename)
        solve(my_mesh,filename,100, isImplicit)
    else :
        raise ValueError("WaveSystemUpwind.py expects a mesh file name and a boolean (isImplicit)")
