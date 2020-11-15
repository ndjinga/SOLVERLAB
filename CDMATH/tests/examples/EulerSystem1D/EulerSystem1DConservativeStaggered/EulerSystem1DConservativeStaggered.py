#!/usr/bin/env python3
# -*-coding:utf-8 -*

# Isothermal Euler system
# d rho/d t + d q/d x =0
# d q/d t + d (q^2/rho+p)/d x = 0, where q = rho*u and p = c^2*rho
# UU = (rho,q) : conservative variable
# Scheme : Conservative stagerred scheme (Ait Ameur et al)
# Author :  Katia Ait Ameur
# Date : November 2020

import cdmath
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from math import sqrt

c0=330.#reference sound speed for water at 15 bars
precision=1e-5
rho0=1

def initial_conditions_Riemann_problem(a,b,nx):
    print("Initial data Riemann problem")
    dx = (b - a) / nx #space step
    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    rho_initial = [ (xi<(a+b)/2)*rho0        + (xi>=(a+b)/2)*rho0*2       for xi in x]
    q_initial   = [ (xi<(a+b)/2)*rho0*(-300) + (xi>=(a+b)/2)*rho0*2*(-300)  for xi in x]

    return rho_initial, q_initial

def matrix_coef(rho_l,q_l,rho_r,q_r):
    m_coef = -0.5*(q_l-q_r)/rho_r
    return m_coef

def jacobianMatricesm(coeff,rho_l,q_l,rho_r,q_r):
    RoeMat   = cdmath.Matrix(2,2);
    Dmac     = cdmath.Matrix(2,2);   
    
    u_l=q_l/rho_l  
    u_r=q_r/rho_r
    Dmac_coef = matrix_coef(rho_l,q_l,rho_r,q_r)
    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   
    RoeMat[0,0]   = 0
    RoeMat[0,1]   = 1
    RoeMat[1,0]   = c0*c0 - u*u
    RoeMat[1,1]   = 2*u
    
    Dmac[0,0]= abs(u)-u;
    Dmac[0,1]= 1;
    Dmac[1,0]= -c0*c0-u*u;
    Dmac[1,1]= abs(u)+u; 

    return (RoeMat-Dmac)*coeff*0.5
 
def jacobianMatricesp(coeff,rho_l,q_l,rho_r,q_r):
    RoeMat   = cdmath.Matrix(2,2);
    Dmac = cdmath.Matrix(2,2);   
    
    u_l=q_l/rho_l 
    u_r=q_r/rho_r
    Dmac_coef = matrix_coef(rho_l,q_l,rho_r,q_r)
    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   
    RoeMat[0,0]   = 0
    RoeMat[0,1]   = 1
    RoeMat[1,0]   = c0*c0 - u*u
    RoeMat[1,1]   = 2*u
    
    Dmac[0,0]= abs(u)-u;
    Dmac[0,1]= 1;
    Dmac[1,0]= -c0*c0-u*u;
    Dmac[1,1]= abs(u)+u; 

    return (RoeMat+Dmac)*coeff*0.5

def EulerSystemStaggered(ntmax, tmax, cfl, a,b,nx, output_freq, meshName):
    dim=1
    nbComp=dim+1
    nbCells = nx
    dt = 0.
    dx=(b-a)/nx
    dt = cfl * dx / c0
    nbVoisinsMax=2

    time = 0.
    it=0;
    iterGMRESMax=50
    newton_max = 50
    isStationary=False

    #iteration vectors
    Un =cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    dUk=cdmath.Vector(nbCells*(dim+1))
    Rhs=cdmath.Vector(nbCells*(dim+1))

    # Initial conditions #
    print("Construction of the initial condition …")
    rho_field, q_field = initial_conditions_Riemann_problem(a,b,nx)
    v_field   = [   q_field[i]/rho_field[i]  for i in range(nx)]
    p_field   = [ rho_field[i]*(c0*c0)       for i in range(nx)]

    max_initial_rho=max(rho_field)
    min_initial_rho=min(rho_field)
    max_initial_q=max(q_field)
    min_initial_q=min(q_field)
    max_initial_p=max(p_field)
    min_initial_p=min(p_field)
    max_initial_v=max(v_field)
    min_initial_v=min(v_field)

    for k in range(nbCells):
        Un[k*nbComp+0] = rho_field[k]
        Un[k*nbComp+1] =   q_field[k]

    print("Starting computation of the non linear Euler system with staggered scheme …")
    divMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    # Picture settings
    fig, ([axDensity, axMomentum],[axVelocity, axPressure]) = plt.subplots(2, 2,sharex=True, figsize=(10,10))
    fig.suptitle('Conservative staggered scheme')
    lineDensity, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field, label='Density') #new picture for video # Returns a tuple of line objects, thus the comma
    axDensity.set(xlabel='x (m)', ylabel='Density')
    axDensity.set_xlim(a,b)
    axDensity.set_ylim(min_initial_rho - 0.1*(max_initial_rho-min_initial_rho), max_initial_rho +  0.1*(max_initial_rho-min_initial_rho) )
    axDensity.legend()
    lineMomentum, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field, label='Momentum')
    axMomentum.set(xlabel='x (m)', ylabel='Momentum')
    axMomentum.set_xlim(a,b)
    axMomentum.set_ylim(min_initial_q - 0.1*(max_initial_q-min_initial_q), max_initial_q +  0.1*(max_initial_q-min_initial_q) )
    axMomentum.legend()
    lineVelocity, = axVelocity.plot([a+0.5*dx + i*dx for i in range(nx)], v_field, label='Velocity')
    axVelocity.set(xlabel='x (m)', ylabel='Velocity')
    axVelocity.set_xlim(a,b)
    axVelocity.set_ylim(min_initial_v - 0.4*abs(min_initial_v), max_initial_v +  0.05*abs(max_initial_v) )
    axVelocity.legend()
    linePressure, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field, label='Pressure')
    axPressure.set(xlabel='x (m)', ylabel='Pressure')
    axPressure.set_xlim(a,b)
    axPressure.set_ylim(min_initial_p - 0.05*abs(min_initial_p), max_initial_p +  0.05*abs(max_initial_p) )
    axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axPressure.legend()
 
    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Conservative staggered scheme for the 1D isothermal Euler System", artist = "CEA Saclay", comment="Shock tube")
    writer=FFMpegWriter(fps=10, metadata=metadata, codec='h264')
    with writer.saving(fig, "1DEuler_System_ConservativeStaggered"+".mp4", ntmax):
        writer.grab_frame()
        plt.savefig("EulerSystem"+str(dim)+"D_ConservativeStaggered"+meshName+"_0"+".png")

        # Starting time loop
        while (it<ntmax and time <= tmax and not isStationary):
            dUn = Un.deepCopy()
            Uk  = Un.deepCopy()
            residu = 1.
            k=0
            while (k<newton_max and residu > precision ):
                #DEBUT BOUCLE NEWTON
                divMat.zeroEntries()#sets the matrix coefficients to zero
                for j in range(nbCells):# 
                    if ( j==0) : 
                        rho_r = Uk[j*nbComp+0]
                        q_r   = Uk[j*nbComp+1]
                        rho_l = rho_r # Conditions de Neumann
                        q_l   =   q_r
                        Am= jacobianMatricesm(dt/dx,rho_l,q_l,rho_r,q_r)
                        divMat.addValue(j*nbComp,(j+1)*nbComp,Am)
                        divMat.addValue(j*nbComp,    j*nbComp,Am*(-1.))
                        dUi=cdmath.Vector(2)
                        dUi[0] = Uk[(j+1)*nbComp+0]-Uk[j*nbComp+0]
                        dUi[1] = Uk[(j+1)*nbComp+1]-Uk[j*nbComp+1]
                        temp=cdmath.Vector(2)
                        temp = Am*dUi
                        #print("Bloc 0 matrix  :   ", Am)
                        Rhs[j*nbComp+0] = -temp[0]-(Uk[j*nbComp+0]-Un[j*nbComp+0]) 
                        Rhs[j*nbComp+1] = -temp[1]-(Uk[j*nbComp+1]-Un[j*nbComp+1]) 
                    elif ( j==nbCells-1) :
                        rho_l = Uk[j*nbComp+0]
                        q_l   = Uk[j*nbComp+1]
                        rho_r = rho_l # Conditions de Neumann
                        q_r   =   q_l
                        Ap= jacobianMatricesp(dt/dx,rho_l,q_l,rho_r,q_r)
                        divMat.addValue(j*nbComp,    j*nbComp,Ap)
                        divMat.addValue(j*nbComp,(j-1)*nbComp,Ap*(-1.))
                        dUi=cdmath.Vector(2)
                        dUi[0] = Uk[(j-1)*nbComp+0]-Uk[j*nbComp+0]
                        dUi[1] = Uk[(j-1)*nbComp+1]-Uk[j*nbComp+1]
                        temp=cdmath.Vector(2)
                        temp = Ap*dUi
                        Rhs[j*nbComp+0] = temp[0]-(Uk[j*nbComp+0]-Un[j*nbComp+0]) 
                        Rhs[j*nbComp+1] = temp[1]-(Uk[j*nbComp+1]-Un[j*nbComp+1]) 
                    else :
                        rho_l = Uk[(j-1)*nbComp+0]
                        q_l   = Uk[(j-1)*nbComp+1]
                        rho_r = Uk[j*nbComp+0]
                        q_r   = Uk[j*nbComp+1]
                        Ap = jacobianMatricesp(dt/dx,rho_l,q_l,rho_r,q_r)
                        ###############################################################
                        rho_l = Uk[j*nbComp+0]
                        q_l   = Uk[j*nbComp+1]
                        rho_r = Uk[(j+1)*nbComp+0]
                        q_r   = Uk[(j+1)*nbComp+1]
                        Am = jacobianMatricesm(dt/dx,rho_l,q_l,rho_r,q_r)
                        divMat.addValue(j*nbComp,(j+1)*nbComp,Am)
                        divMat.addValue(j*nbComp,    j*nbComp,Am*(-1.))
                        divMat.addValue(j*nbComp,    j*nbComp,Ap)
                        divMat.addValue(j*nbComp,(j-1)*nbComp,Ap*(-1.))
                        dUi1=cdmath.Vector(2)
                        dUi2=cdmath.Vector(2)
                        dUi1[0] = Uk[(j+1)*nbComp+0]-Uk[j*nbComp+0]
                        dUi1[1] = Uk[(j+1)*nbComp+1]-Uk[j*nbComp+1]
                        dUi2[0] = Uk[(j-1)*nbComp+0]-Uk[j*nbComp+0]
                        dUi2[1] = Uk[(j-1)*nbComp+1]-Uk[j*nbComp+1] 
                        temp1 = cdmath.Vector(2)
                        temp2 = cdmath.Vector(2)
                        temp1 = Am*dUi1
                        temp2 = Ap*dUi2
                        Rhs[j*nbComp+0] = -temp1[0] + temp2[0] -(Uk[j*nbComp+0]-Un[j*nbComp+0]) 
                        Rhs[j*nbComp+1] = -temp1[1] + temp2[1] -(Uk[j*nbComp+1]-Un[j*nbComp+1]) 
                divMat.diagonalShift(1)#only after  filling all coefficients
                LS=cdmath.LinearSolver(divMat,Rhs,iterGMRESMax, precision, "GMRES","LU")
                dUk=LS.solve(); 
                residu = dUk.norm()
                #stop
                Uk+=dUk
                if(not LS.getStatus()):
                    print("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
                    raise ValueError("Pas de convergence du système linéaire");
                k=k+1
            Un = Uk.deepCopy()
            dUn-=Un
            for k in range(nbCells):
                rho_field[k] = Un[k*(dim+1)+0]
                q_field[k]   = Un[k*(dim+1)+1] 
            v_field   = [   q_field[i]/rho_field[i]  for i in range(nx)]
            p_field   = [ rho_field[i]*(c0*c0)       for i in range(nx)]

            lineDensity.set_ydata(rho_field)
            lineMomentum.set_ydata(q_field)
            lineVelocity.set_ydata(v_field)
            linePressure.set_ydata(p_field)
            writer.grab_frame()

            time=time+dt;
            it=it+1;
            #Sauvegardes
            if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                #print("Last linear system converged in ", LS.getNumberOfIter(), " GMRES iterations", ", residu final:   ",residu)

                plt.savefig("EulerSystem"+str(dim)+"D_ConservativeStaggered"+meshName+"_"+str(it)+".png")
                print
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    if(it>=ntmax):
        print("Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
        return
    elif(isStationary):
        print("Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print("------------------------------------------------------------------------------------")
        plt.savefig("EulerSystem"+str(dim)+"Staggered"+meshName+"_Stat.png")
        return
    else:
        print("Temps maximum Tmax= ", tmax, " atteint")
        return

def solve( a,b,nx, meshName, meshType, cfl):
    print("Resolution of the Euler system in dimension 1 on "+str(nx)+ " cells")
    print("Initial data : ", "Riemann problem")
    print("Boundary conditions : ", "Neumann")
    print("Mesh name : ",meshName , ", ", nx, " cells")
    # Problem data
    tmax = 10.
    ntmax = 50
    output_freq = 10
    EulerSystemStaggered(ntmax, tmax, cfl, a,b,nx, output_freq, meshName)
    return

if __name__ == """__main__""":
    a=0.
    b=1.
    nx=100
    cfl=0.99
    solve( a,b,nx,"SquareRegularSquares","RegularSquares",cfl)
