#!/usr/bin/env python
# -*-coding:utf-8 -*

import CDMATH
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import sys

p0=155.e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_Riemann_problem(a,b,nx):
    print("Initial data Riemann problem")

    dx = (b - a) / nx #space step
    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)

    u_initial = [ 0 ]*nx
    p_initial = [ (xi<(a+b)/2)*p0 + (xi>=(a+b)/2)*p0/2  for xi in x]

    return p_initial, u_initial

def jacobianMatrices(coeff,scaling):
    dim=1
    A=CDMATH.Matrix(dim+1,dim+1)
    absA=CDMATH.Matrix(dim+1,dim+1)

    absA[0,0]=c0*coeff
    for i in range(dim):
        for j in range(dim):
            absA[i+1,j+1]=c0*coeff
        if( scaling==0):
            A[0,i+1]=c0*c0*coeff
            A[i+1,0]=      coeff
        elif( scaling==1):
            A[0,i+1]=      coeff
            A[i+1,0]=      coeff
        else:
            A[0,i+1]=   c0*coeff
            A[i+1,0]=   c0*coeff
       
    return A,absA
    
    
def computeUpwindDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling):
    nbCells = nx
    dx=(b-a)/nx
    dim=1
    nbComp=dim+1
    normal=CDMATH.Vector(dim)

    implMat=CDMATH.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    A,absA= jacobianMatrices(dt/dx,scaling)
    for j in range(nbCells):#On parcourt les cellules
        if ( j==0) :
            implMat.addValue(j*nbComp,(j+1)*nbComp,(A-absA)/2)
            implMat.addValue(j*nbComp,    j*nbComp,(A-absA)/2*(-1.))
        elif ( j==nbCells-1) :
            implMat.addValue(j*nbComp,    j*nbComp,(A+absA)/2)
            implMat.addValue(j*nbComp,(j-1)*nbComp,(A+absA)/2*(-1.))
        else :
            implMat.addValue(j*nbComp,(j+1)*nbComp,(A-absA)/2)
            implMat.addValue(j*nbComp,    j*nbComp,(A-absA)/2*(-1.))

            implMat.addValue(j*nbComp,    j*nbComp,(A+absA)/2)
            implMat.addValue(j*nbComp,(j-1)*nbComp,(A+absA)/2*(-1.))
                
    return implMat

def staggeredMatrices(coeff,scaling):
    dim=1
    S1=CDMATH.Matrix(dim+1,dim+1)
    S2=CDMATH.Matrix(dim+1,dim+1)

    for i in range(dim):
        if( scaling==0):
            S1[0,i+1]=c0*c0*coeff
            S2[i+1,0]=      coeff
        else:
            S1[0,i+1]=   c0*coeff
            S2[i+1,0]=   c0*coeff
       
    return S1,S2
    
    
def computeStaggeredDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling):
    nbCells = nx
    dx=(b-a)/nx
    dim=1
    nbComp=dim+1
    normal=CDMATH.Vector(dim)

    implMat=CDMATH.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    S1,S2 = staggeredMatrices(dt/dx,scaling)
    for k in range(nbCells):#On parcourt les cellules
        if ( k==0) :
            implMat.addValue(k*nbComp, (k+1)*nbComp, S1) 
            implMat.addValue(k*nbComp,  k   *nbComp, S1*(-1.))
        elif ( k==nbCells-1) :
            implMat.addValue(k*nbComp,  k   *nbComp, S2)
            implMat.addValue(k*nbComp, (k-1)*nbComp, S2*(-1.))
        else :
            implMat.addValue(k*nbComp, (k+1)*nbComp, S1)
            implMat.addValue(k*nbComp,  k   *nbComp, S1*(-1.))
            
            implMat.addValue(k*nbComp,  k   *nbComp, S2)
            implMat.addValue(k*nbComp, (k-1)*nbComp, S2*(-1.))
        
    return implMat

def WaveSystemVF(ntmax, tmax, cfl, a,b,nx, output_freq, meshName,scaling):
    dim=1
    nbCells = nx
    
    dt = 0.
    time = 0.
    it=0. # ; ?
    isStationary=False
    
    dx=(b-a)/nx
    dt = cfl * dx / c0

    nbVoisinsMax=2
    isImplicit=False
    
    #iteration vectors
    Un =CDMATH.Vector(nbCells*(dim+1))
    dUn=CDMATH.Vector(nbCells*(dim+1))
    
    Un_staggered =CDMATH.Vector(nbCells*(dim+1))
    dUn_staggered=CDMATH.Vector(nbCells*(dim+1))

    Un_implicite =CDMATH.Vector(nbCells*(dim+1))
    dUn_implicite=CDMATH.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_Riemann_problem(a,b,nx)
    pressure_field_staggered, velocity_field_staggered = initial_conditions_Riemann_problem(a,b,nx)
    pressure_field_implicite, velocity_field_implicite = initial_conditions_Riemann_problem(a,b,nx)
    max_initial_p=max(pressure_field)
    min_initial_p=min(pressure_field)
    max_initial_v=max(velocity_field)
    min_initial_v=min(velocity_field)
    
    for k in range(nbCells):
        Un[k*(dim+1)+0] =     pressure_field[k]
        Un[k*(dim+1)+1] =rho0*velocity_field[k]
        Un_staggered[k*(dim+1)+0] =     pressure_field_staggered[k]
        Un_staggered[k*(dim+1)+1] =rho0*velocity_field_staggered[k]
        Un_implicite[k*(dim+1)+0] =     pressure_field_implicite[k]
        Un_implicite[k*(dim+1)+1] =rho0*velocity_field_implicite[k]

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Finite volumes schemes for the 1D Wave System", artist = "CEA Saclay", comment="Shock propagation")
    writer=FFMpegWriter(fps=10, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DWaveSystem_FV"+".mp4", ntmax):
        #sauvegarde de la donnée initiale
        plt.xlabel('x (m)')
        plt.ylabel('Pressure -Pa)')
        plt.xlim(a,b)
        plt.ylim( min_initial_p - 0.1*(max_initial_p-min_initial_p), max_initial_p +  0.1*(max_initial_p-min_initial_p) )
        plt.title("Riemann problem for Wave system on " + str(nx) + " cells")
        line1, = plt.plot([a+0.5*dx + i*dx for i in range(nx)], pressure_field, label='Explicit upwind scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line2, = plt.plot([a+0.5*dx + i*dx for i in range(nx)], pressure_field_staggered, label='Staggered scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line3, = plt.plot([a+0.5*dx + i*dx for i in range(nx)], pressure_field_implicite, label='Implicit upwind scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        plt.legend()
        writer.grab_frame()
        plt.savefig("WaveSystem"+str(dim)+"DUpwindExplicit"+meshName+"_pressure"+"_0"+".png")
        np.savetxt( "WaveSystem"+str(dim)+"DUpwindExplicit"+meshName+"_pressure"+"_0"+".txt", pressure_field, delimiter="\n")
        np.savetxt( "WaveSystem"+str(dim)+"DUpwindExplicit"+meshName+"_velocity"+"_0"+".txt", velocity_field, delimiter="\n")
        plt.savefig("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+"_0"+".png")
        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+"_0"+".txt", pressure_field, delimiter="\n")
        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity"+"_0"+".txt", velocity_field, delimiter="\n")
        plt.savefig("WaveSystem"+str(dim)+"DUpwindImplicit"+meshName+"_pressure"+"_0"+".png")
        np.savetxt( "WaveSystem"+str(dim)+"DUpwindImplicit"+meshName+"_pressure"+"_0"+".txt", pressure_field, delimiter="\n")
        np.savetxt( "WaveSystem"+str(dim)+"DUpwindImplicit"+meshName+"_velocity"+"_0"+".txt", velocity_field, delimiter="\n")
        
        divMat=computeUpwindDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling)
        divMat_staggered=computeStaggeredDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling)
        divMat_implicit=computeUpwindDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling)
            
        iterGMRESMax=50
    
        divMat_implicit.diagonalShift(1)#only after  filling all coefficients
        if( scaling==0):
            LS=CDMATH.LinearSolver(divMat_implicit,Un,iterGMRESMax, precision, "GMRES","ILU")
        else:
            LS=CDMATH.LinearSolver(divMat_implicit,Un,iterGMRESMax, precision, "GMRES","ILU") #Vn 
    
        divMat_staggered.diagonalShift(1)#only after  filling all coefficients
        LS_staggered=CDMATH.LinearSolver(divMat_staggered,Un_staggered,iterGMRESMax, precision, "GMRES","ILU")
    
        print("Starting computation of the linear wave system with upwind and staggered scheme …")
        
        # Starting time loop
        while (it<ntmax and time <= tmax and not isStationary):
            dUn=divMat*Un
            Un-=dUn
    
            dUn_implicite=Un_implicite.deepCopy()
            LS.setSndMember(Un_implicite)
            Un_implicite=LS.solve();
            if(not LS.getStatus()):
                print ("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            dUn_implicite-=Un_implicite
    
            dUn_staggered=Un_staggered.deepCopy()
            LS_staggered.setSndMember(Un_staggered)
            Un_staggered=LS_staggered.solve();
            if(not LS_staggered.getStatus()):
                print ("Linear system did not converge for staggered scheme ", LS.getNumberOfIter(), " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            dUn_staggered-=Un_staggered
    
            for k in range(nbCells):
                pressure_field[k] = Un[k*(dim+1)+0]
                velocity_field[k] = Un[k*(dim+1)+1] / rho0
                pressure_field_staggered[k] = Un_staggered[k*(dim+1)+0]
                velocity_field_staggered[k] = Un_staggered[k*(dim+1)+1] / rho0
                pressure_field_implicite[k] = Un_implicite[k*(dim+1)+0]
                velocity_field_implicite[k] = Un_implicite[k*(dim+1)+1] / rho0
    
            line1.set_ydata(pressure_field)
            line2.set_ydata(pressure_field_staggered)
            line3.set_ydata(pressure_field_implicite)
            writer.grab_frame()
    
            time=time+dt;
            it=it+1;
        
            #Sauvegardes
            if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                if(isImplicit):
                    print("Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations")
    
                np.savetxt("WaveSystem"+str(dim)+"DUpwind"   +meshName+"_pressure"+str(it)+".txt", pressure_field          , delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DUpwind"   +meshName+"_velocity"+str(it)+".txt", velocity_field          , delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+str(it)+".txt", pressure_field_staggered, delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity"+str(it)+".txt", velocity_field_staggered, delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DUpwindImplicit"   +meshName+"_pressure"+str(it)+".txt", pressure_field          , delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DUpwindImplicit"   +meshName+"_velocity"+str(it)+".txt", velocity_field          , delimiter="\n")
                plt.savefig("WaveSystem"+str(dim)+"DUpwindStaggered"+meshName+"_pressure"+str(it)+".png")
    
                print
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

    if(it>=ntmax):
        print("Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
        return
    elif(isStationary):
        print("Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print("------------------------------------------------------------------------------------")

        np.savetxt("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat.txt", pressure_field, delimiter="\n")
        np.savetxt("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat.txt", velocity_field, delimiter="\n")
        plt.savefig("WaveSystem"+str(dim)+"DUpwindStaggered"+meshName+"_pressure_Stat.png")

        return
    else:
        print("Temps maximum Tmax= ", tmax, " atteint")
        return


def solve( a,b,nx, meshName, scaling, meshType, cfl):

    print("Resolution of the Wave system in dimension 1 on "+str(nx)+ " cells")
    print("Initial data : ", "Riemann problem")
    print("Boundary conditions : ", "Neumann")
    print("Mesh name : ",meshName , ", ", nx, " cells")
    
    # Problem data
    tmax = 10000.
    ntmax = 50
    output_freq = 1

    WaveSystemVF(ntmax, tmax, cfl, a,b,nx, output_freq, meshName, scaling)
    
    return
    

if __name__ == """__main__""":
    a=0.
    b=1.
    nx=100
    cfl=0.99
    scaling=0
    solve( a,b,nx,"SquareRegularSquares",scaling,"RegularSquares",cfl)

