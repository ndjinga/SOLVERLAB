#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système Euler 2D barotrope sans terme source
#                \partial_t rho + \div q = 0
#                \partial_t q   + \div q\otimes q/rho  +  \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Propagation d'une onde de choc sphérique
#               Utilisation du schéma de Roe upwind explicite ou implicite sur un maillage général
#               Initialisation par une surpression sphérique
#               Conditions aux limites parois
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sqrt
import cdmath
import PV_routines
import VTK_routines
import sys

p0=155e5#reference pressure in a pressurised nuclear vessel (used for initial data)
c0=700.#reference sound speed for water at 155 bars, 600K (used for eox p= rho c0**2
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
            pressure_field[i] = p0
            pass
        pass

    return pressure_field, velocity_field

    
def jacobianMatrices(normale,coeff,rho_l,q_l,rho_r,q_r):
    RoeMat   = cdmath.Matrix(3,3);
    AbsRoeMa = cdmath.Matrix(3,3);

    tangent=cdmath.Vector(2);
    tangent[0]= normale[1];
    tangent[1]=-normale[0];

    u_l = cdmath.Vector(2); u_l[0]*=q_l[0]/rho_l; u_l[1]*=q_l[1]/rho_l;
    u_r = cdmath.Vector(2); u_r[0]*=q_r[0]/rho_r; u_r[1]*=q_r[1]/rho_r;
    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r )
        raise ValueError("Negative density")
    u=cdmath.Vector(2);
    u[0] = (u_l[0]*sqrt(rho_l)+u_r[0]*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   
    u[1] = (u_l[1]*sqrt(rho_l)+u_r[1]*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   
    un=u*normale;

    RoeMat[0,0]   = 0
    RoeMat[0,1]   =     normale[0]
    RoeMat[0,2]   =     normale[1]
    RoeMat[1,0]   = c0*c0*normale[0] - un*u[0]
    RoeMat[2,0]   = c0*c0*normale[1] - un*u[1]
    RoeMat[1,1]   = un + normale[0]*u[0]
    RoeMat[1,2]   =      normale[1]*u[0]
    RoeMat[2,2]   = un + normale[1]*u[1]
    RoeMat[2,1]   =      normale[0]*u[1]

    AbsRoeMa[0,0]=(abs(un-c0)*(un+c0)+abs(un+c0)*(c0-un))/(2*c0);
    AbsRoeMa[0,1]= (abs(un+c0)-abs(un-c0))/(2*c0)*normale[0];
    AbsRoeMa[0,2]= (abs(un+c0)-abs(un-c0))/(2*c0)*normale[1];
    AbsRoeMa[1,0]=(abs(un-c0)*(un+c0)*(u[0]-c0*normale[0])-abs(un+c0)*(un-c0)*(u[0]+c0*normale[0]))/(2*c0)-abs(un)*(u*tangent)*tangent[0];
    AbsRoeMa[2,0]=(abs(un-c0)*(un+c0)*(u[1]-c0*normale[1])-abs(un+c0)*(un-c0)*(u[1]+c0*normale[1]))/(2*c0)-abs(un)*(u*tangent)*tangent[1];
    #subMatrix=(abs(un+c0)*((u-c0*normale)^normale)-abs(un-c0)*((u-c0*normale)^normale))/(2*c0)+abs(un)*(tangent^tangent);
    AbsRoeMa[1,1]=(abs(un+c0)*((u[0]-c0*normale[0])*normale[0])-abs(un-c0)*((u[0]-c0*normale[0])*normale[0]))/(2*c0)+abs(un)*(tangent[0]*tangent[0]);#subMatrix[0,0];
    AbsRoeMa[1,2]=(abs(un+c0)*((u[0]-c0*normale[0])*normale[1])-abs(un-c0)*((u[0]-c0*normale[0])*normale[1]))/(2*c0)+abs(un)*(tangent[0]*tangent[1]);#subMatrix[0,1];
    AbsRoeMa[2,1]=(abs(un+c0)*((u[1]-c0*normale[1])*normale[0])-abs(un-c0)*((u[1]-c0*normale[1])*normale[0]))/(2*c0)+abs(un)*(tangent[1]*tangent[0]);
    AbsRoeMa[2,2]=(abs(un+c0)*((u[1]-c0*normale[1])*normale[1])-abs(un-c0)*((u[1]-c0*normale[1])*normale[1]))/(2*c0)+abs(un)*(tangent[1]*tangent[1]);

    return (RoeMat-AbsRoeMa)*coeff*0.5,un;

def computeDivergenceMatrix(my_mesh,implMat,Un):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1
    normal=cdmath.Vector(dim)
    maxAbsEigVa = 0
    q_l=cdmath.Vector(2);
    q_r=cdmath.Vector(2);

    idMoinsJacCL=cdmath.Matrix(nbComp)
    
    for j in range(nbCells):#On parcourt les cellules
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

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
                    
                q_l[0]=Un[j*nbComp+1]
                q_l[1]=Un[j*nbComp+2]
                q_r[0]=Un[cellAutre*nbComp+1]
                q_r[1]=Un[cellAutre*nbComp+2]
                Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(),Un[j*nbComp],q_l,Un[cellAutre*nbComp],q_r);
            
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                if( Fk.getGroupName() != "Periodic" and Fk.getGroupName() != "Neumann"):#Wall boundary condition unless Periodic/Neumann specified explicitly
                    q_l[0]=Un[j*nbComp+1]
                    q_l[1]=Un[j*nbComp+2]
                    q_r=q_l-normal*2*(q_l*normal)
                    Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(),Un[j*nbComp],q_l,Un[j*nbComp],q_r);
            
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
                    Am, un=jacobianMatrices( normal,Fk.getMeasure()/Cj.getMeasure(),Un[j*nbComp],q_l,Un[cellAutre*nbComp],q_r);
            
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
                elif( Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print( Fk.getGroupName() )
                    raise ValueError("computeFluxes: Unknown boundary condition name");
            
            maxAbsEigVa = max(maxAbsEigVa,abs(un+c0),abs(un-c0));
    
    return maxAbsEigVa

def EulerSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution, isImplicit):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    
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
        print(Un[k*(dim+1)+0],Un[k*(dim+1)+1],Un[k*(dim+1)+2])
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity");

    dx_min=my_mesh.minRatioVolSurf()

    divMat=cdmath.SparseMatrixPetsc(nbCells*(1+dim),nbCells*(1+dim),(nbVoisinsMax+1)*(1+dim))
    if( isImplicit):
        iterGMRESMax=50
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")

    print("Starting computation of the linear wave system with an explicit UPWIND scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        divMat.zeroEntries()#sets the matrix coefficients to zero
        vp_max=computeDivergenceMatrix(my_mesh,divMat,Un)#To update at each time step
        dt = cfl * dx_min / vp_max#To update at each time step
            
        if(isImplicit):
            #Adding the identity matrix on the diagonal
            divMat.diagonalShift(1)#only after  filling all coefficients
            dUn=Un.deepCopy()
            LS.setSndMember(Un)
            Un=LS.solve();
            if(not LS.getStatus()):
                print( "Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un

        else:
            dUn=divMat*Un
#            print(Un)
            Un-=dUn
        
        time=time+dt;
        it=it+1;
 
         #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

            for k in range(nbCells):
                pressure_field[k]  =Un[k*(dim+1)+0]*c0*c0
                velocity_field[k,0]=Un[k*(dim+1)+1]/Un[k*(dim+1)+0]
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/Un[k*(dim+1)+0]
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/Un[k*(dim+1)+0]
#                print(Un[k*(dim+1)+0],Un[k*(dim+1)+1],Un[k*(dim+1)+2])

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity",False);

    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    print()

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time)

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
        diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"EulerSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("EulerSystem"+str(dim)+"DUpwind"+"_isImplicit"+str(isImplicit)+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"EulerSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat")
        
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint")


def solve(my_mesh,filename,resolution, isImplicit):
    print( "Resolution of the Euler system in dimension ", my_mesh.getSpaceDimension())
    if( my_mesh.getSpaceDimension()!=2):
        raise ValueError("Only dimension 2 simulations allowed")
    print( "Numerical method : upwind")
    print( "Initial data : spherical wave")
    print( "Wall boundary conditions")
    print( "Mesh name : ",filename , my_mesh.getNumberOfCells(), " cells")

    # Problem data
    tmax = 1.
    ntmax = 1
    cfl = 1./my_mesh.getSpaceDimension()
    output_freq = 1

    EulerSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution, isImplicit)

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
        raise ValueError("EulerSystemUpwind.py expects a mesh file name and a boolean (isImplicit)")
