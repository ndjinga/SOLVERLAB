#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import time, json
import cdmath
import PV_routines
import VTK_routines

precision=1e-5

def initial_conditions_transport_equation(my_mesh):
    print( "Initial_data","Shock")
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    initial_field = cdmath.Field("unknown", cdmath.CELLS, my_mesh, 1)

    rayon = 0.15
    xcentre = 0.5
    ycentre = 0.5
    zcentre = 0.5

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()

        valX = (x - xcentre) * (x - xcentre)
        
        if(dim==1):
            val=sqrt(valX)
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
            initial_field[i] = 1
            pass
        else:
            initial_field[i] = 0
            pass
        pass
        
    return initial_field

def upwinding_coeff(normal, coeff, velocity):
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

def TransportEquationVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution,test_bc,velocity,isImplicit):
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
    stat_field = initial_conditions_transport_equation(my_mesh)
    for k in range(nbCells):
        Un[k] =     stat_field[k]
    unknown_field = initial_conditions_transport_equation(my_mesh)
    S = cdmath.Vector(nbCells)#source term is zero
            
    total_mass_initial=unknown_field.integral()#For conservation test later
    
    #sauvegarde de la donnée initiale
    unknown_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName);

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
        #LS.setComputeConditionNumber()#Computes the condition number of the preconditioned operator
        
    print("Starting computation of the linear transport equation with an Upwind scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            LS.setSndMember(Un+S*dt)
            Un=LS.solve();
            if(not LS.getStatus()):
                print("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

        else:
            dUn=divMat*Un-S*dt
            Un-=dUn
        
        maxVector=dUn.maxVector(1)
        isStationary= maxVector[0]<precision 

        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
            print( "Variation temporelle relative : ", maxVector[0] )
            if(isImplicit):
                print( "Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations" )

            for k in range(nbCells):
                unknown_field[k]  =Un[k]

            unknown_field.setTime(time,it);
            unknown_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName,False);

    print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
    print( "Variation temporelle relative : ", maxVector[0] )

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint" )
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        if(test_bc=="Periodic"):
            print( "Mass loss: ", (total_mass_initial-unknown_field.integral()).norm(), " precision required= ", precision )
            assert (total_mass_initial-unknown_field.integral()).norm()<precision
        print( "------------------------------------------------------------------------------------")

        unknown_field.setTime(time,0);
        unknown_field.writeVTK("TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat");

        #Postprocessing : Extraction of the diagonal data
        if(dim==2):
            diag_data_u=VTK_routines.Extract_field_data_over_line_to_numpyArray(unknown_field,[0,1,0],[1,0,0], resolution)    
        elif(dim==3):
            diag_data_u=VTK_routines.Extract_field_data_over_line_to_numpyArray(unknown_field,[0,0,0],[1,1,1], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat"+'_0.vtu',"unknown",'CELLS',"TransportEquation"+str(dim)+"DUpwind"+meshName+"_Stat")
        
        return nbCells, time, it, unknown_field.getNormEuclidean().max(), diag_data_u        
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint" )
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh, meshName, resolution, meshType, cfl, test_bc):
    print( "Resolution of the Transport Equation in dimension ", my_mesh.getMeshDimension() )
    print( "Numerical method : ", "Upwind" )
    print( "Initial data : ", "Spherical shock" )
    print( "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells" )
    
    # Problem data
    tmax = 10000.
    ntmax = 30000
    output_freq = 1000

    isImplicit=True
    
    dim=my_mesh.getMeshDimension()
    velocity=cdmath.Vector(dim)
    for i in range(dim) :
        velocity[i] = 1

    nbCells, t_final, ndt_final, max_unknown, diag_data_u = TransportEquationVF(ntmax, tmax, cfl, my_mesh, output_freq, meshName, resolution, test_bc, velocity,isImplicit)

    return nbCells, t_final, ndt_final, max_unknown, diag_data_u

def solve_file( filename,meshName, resolution,meshType, cfl, test_bc):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, meshName+str(my_mesh.getNumberOfCells()),resolution, meshType, cfl, test_bc)
    

if __name__ == """__main__""":
    M1=cdmath.Mesh(0.,1.,20,0.,1.,20,0)
    cfl=0.5
    solve(M1,"SquareRegularTriangles",100,"Regular triangles",cfl,"Periodic")

