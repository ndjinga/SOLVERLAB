#!/usr/bin/env python3
# -*-coding:utf-8 -*

import cdmath
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
    print( "Initial data Riemann problem" )

    dx = (b - a) / nx #space step
    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)

    u_initial = [ 0 ]*nx
    p_initial = [ (xi<(a+b)/2)*p0 + (xi>=(a+b)/2)*p0/2  for xi in x]

    return p_initial, u_initial

def staggeredMatrices(coeff,scaling):
    dim=1
    S1=cdmath.Matrix(dim+1,dim+1)
    S2=cdmath.Matrix(dim+1,dim+1)

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
    normal=cdmath.Vector(dim)

    implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

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
    it=0;
    isStationary=False
    
    dx=(b-a)/nx
    dt = cfl * dx / c0

    nbVoisinsMax=2
    
    #iteration vectors
    Un_staggered =cdmath.Vector(nbCells*(dim+1))
    dUn_staggered=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field_staggered, velocity_field_staggered = initial_conditions_Riemann_problem(a,b,nx)
    max_initial_p=max(pressure_field_staggered)
    min_initial_p=min(pressure_field_staggered)
    max_initial_v=max(velocity_field_staggered)
    min_initial_v=min(velocity_field_staggered)
    
    for k in range(nbCells):
        Un_staggered[k*(dim+1)+0] =     pressure_field_staggered[k]
        Un_staggered[k*(dim+1)+1] =rho0*velocity_field_staggered[k]

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Finite volumes schemes for the 2D Wave System", artist = "CEA Saclay", comment="Shock propagation")
    writer=FFMpegWriter(fps=10, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "2DWaveSystem_Staggered"+".mp4", ntmax):
        #sauvegarde de la donnée initiale
        plt.xlabel('x (m)')
        plt.ylabel('Pressure -Pa)')
        plt.xlim(a,b)
        plt.ylim( min_initial_p - 0.1*(max_initial_p-min_initial_p), max_initial_p +  0.1*(max_initial_p-min_initial_p) )
        plt.title("Riemann problem for Wave system on " + str(nx) + " cells")
        line2, = plt.plot([a+0.5*dx + i*dx for i in range(nx)], pressure_field_staggered, label='Staggered scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        plt.legend()
        writer.grab_frame()
        plt.savefig("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+"_0"+".png")
        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+"_0"+".txt", pressure_field_staggered, delimiter="\n")
        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity"+"_0"+".txt", velocity_field_staggered, delimiter="\n")
        
        divMat_staggered=computeStaggeredDivergenceMatrix(a,b,nx,nbVoisinsMax,dt,scaling)
            
        iterGMRESMax=50
    
        divMat_staggered.diagonalShift(1)#only after  filling all coefficients
        LS_staggered=cdmath.LinearSolver(divMat_staggered,Un_staggered,iterGMRESMax, precision, "GMRES","ILU")
    
        print("Starting computation of the linear wave system with staggered scheme …")
        
        # Starting time loop
        while (it<ntmax and time <= tmax and not isStationary):
            dUn_staggered=Un_staggered.deepCopy()
            LS_staggered.setSndMember(Un_staggered)
            Un_staggered=LS_staggered.solve();
            if(not LS_staggered.getStatus()):
                print( "Linear system did not converge for staggered scheme ", LS.getNumberOfIter(), " GMRES iterations" )
                raise ValueError("Pas de convergence du système linéaire");
            dUn_staggered-=Un_staggered
    
            for k in range(nbCells):
                pressure_field_staggered[k] = Un_staggered[k*(dim+1)+0]
                velocity_field_staggered[k] = Un_staggered[k*(dim+1)+1] / rho0
    
            line2.set_ydata(pressure_field_staggered)
            writer.grab_frame()
    
            time=time+dt;
            it=it+1;
        
            #Sauvegardes
            if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
                print( "-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) )
                print( "Linear system converged in ", LS_staggered.getNumberOfIter(), " GMRES iterations" )
    
                np.savetxt("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+str(it)+".txt", pressure_field_staggered, delimiter="\n")
                np.savetxt("WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity"+str(it)+".txt", velocity_field_staggered, delimiter="\n")
                plt.savefig("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure"+str(it)+".png")
    
                print()
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

    if(it>=ntmax):
        print( "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint" )
        return
    elif(isStationary):
        print( "Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print( "------------------------------------------------------------------------------------")

        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat.txt", pressure_field_staggered, delimiter="\n")
        np.savetxt( "WaveSystem"+str(dim)+"DStaggered"+meshName+"_velocity_Stat.txt", velocity_field_staggered, delimiter="\n")
        plt.savefig("WaveSystem"+str(dim)+"DStaggered"+meshName+"_pressure_Stat.png")

        return
    else:
        print( "Temps maximum Tmax= ", tmax, " atteint")
        return


def solve( a,b,nx, meshName, scaling, meshType, cfl):

    print( "Resolution of the Wave system in dimension 1 on "+str(nx)+ " cells, staggered scheme")
    print( "Initial data : ", "Riemann problem")
    print( "Boundary conditions : ", "Neumann")
    print( "Mesh name : ",meshName , ", ", nx, " cells")
    
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
