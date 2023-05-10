#!/usr/bin/env python3
# -*-coding:utf-8 -*

# Isothermal Euler system
# d rho/d t + d q/d x =0
# d q/d t + d (q^2/rho+p)/d x = 0, where q = rho*u and p = c^2*rho
# UU = (rho,q) : conservative variable
# Conservative schemes : centered, upwind, Conservative stagerred scheme (Ait Ameur et al)
# Author :  Michael Ndjinga, Katia Ait Ameur
# Date : May 2023

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

def RoeMatrix(rho_l,q_l,rho_r,q_r):
    RoeMat   = cdmath.Matrix(2,2);
    
    u_l=q_l/rho_l  
    u_r=q_r/rho_r

    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   

    RoeMat[0,0]   = 0
    RoeMat[0,1]   = 1
    RoeMat[1,0]   = c0*c0 - u*u
    RoeMat[1,1]   = 2*u
    
    return RoeMat
 
def staggeredDMatrix(rho_l,q_l,rho_r,q_r):
    Dmac     = cdmath.Matrix(2,2);   
    
    u_l=q_l/rho_l  
    u_r=q_r/rho_r

    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   

    Dmac[0,0]= abs(u)-u;
    Dmac[0,1]= 1;
    Dmac[1,0]= -c0*c0-u*u;
    Dmac[1,1]= abs(u)+u; 

    return Dmac
 
def upwindDMatrix(rho_l,q_l,rho_r,q_r):
    Dupwind     = cdmath.Matrix(2,2);   
    
    u_l=q_l/rho_l  
    u_r=q_r/rho_r

    if rho_l<0 or rho_r<0 :
        print( "rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   

    Dupwind[0,0]=(abs(u-c0)*(u+c0)+abs(u+c0)*(c0-u))/(2*c0);
    Dupwind[0,1]= (abs(u+c0)-abs(u-c0))/(2*c0);
    Dupwind[1,0]=(abs(u-c0)*(u*u-c0*c0)+abs(u+c0)*(c0*c0-u*u))/(2*c0);
    Dupwind[1,1]=((u+c0)*abs(u+c0)-(u-c0)*abs(u-c0))/(2*c0);

    return Dupwind
 
def centeredDMatrix(rho_l,q_l,rho_r,q_r):
    Dcentered     = cdmath.Matrix(2,2);   
    
    return Dcentered
 
def EulerSystemSchemeComparison(ntmax, tmax, cfl, a,b,nx, output_freq, meshName):
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
    Un_staggered =cdmath.Vector(nbCells*(dim+1))
    dUn_staggered=cdmath.Vector(nbCells*(dim+1))
    dUk_staggered=cdmath.Vector(nbCells*(dim+1))
    Rhs_staggered=cdmath.Vector(nbCells*(dim+1))

    Un_upwind =cdmath.Vector(nbCells*(dim+1))
    dUn_upwind=cdmath.Vector(nbCells*(dim+1))
    dUk_upwind=cdmath.Vector(nbCells*(dim+1))
    Rhs_upwind=cdmath.Vector(nbCells*(dim+1))

    Un_centered =cdmath.Vector(nbCells*(dim+1))
    dUn_centered=cdmath.Vector(nbCells*(dim+1))
    dUk_centered=cdmath.Vector(nbCells*(dim+1))
    Rhs_centered=cdmath.Vector(nbCells*(dim+1))

    # Initial conditions #
    print("Construction of the initial condition …")
    rho_field_upwind, q_field_upwind = initial_conditions_Riemann_problem(a,b,nx)
    v_field_upwind   = [   q_field_upwind[i]/rho_field_upwind[i]  for i in range(nx)]
    p_field_upwind   = [ rho_field_upwind[i]*(c0*c0)       for i in range(nx)]

    rho_field_staggered = rho_field_upwind 
    q_field_staggered = q_field_upwind
    v_field_staggered   = v_field_upwind
    p_field_staggered   = p_field_upwind

    rho_field_centered = rho_field_upwind
    q_field_centered = q_field_upwind
    v_field_centered   = v_field_upwind
    p_field_centered   = p_field_upwind

    max_initial_rho=max(rho_field_upwind)
    min_initial_rho=min(rho_field_upwind)
    max_initial_q=max(q_field_upwind)
    min_initial_q=min(q_field_upwind)
    max_initial_p=max(p_field_upwind)
    min_initial_p=min(p_field_upwind)
    max_initial_v=max(v_field_upwind)
    min_initial_v=min(v_field_upwind)

    for k in range(nbCells):
        Un_upwind[k*nbComp+0] = rho_field_upwind[k]
        Un_upwind[k*nbComp+1] =   q_field_upwind[k]
        Un_staggered[k*nbComp+0] = rho_field_staggered[k]
        Un_staggered[k*nbComp+1] =   q_field_staggered[k]
        Un_centered[k*nbComp+0] = rho_field_centered[k]
        Un_centered[k*nbComp+1] =   q_field_centered[k]

    #Pour les boucles sur les cellules
    dUi=cdmath.Vector(2)
    temp=cdmath.Vector(2)
    dUi1=cdmath.Vector(2)
    dUi2=cdmath.Vector(2)
    temp1 = cdmath.Vector(2)
    temp2 = cdmath.Vector(2)

    print("Starting computation of the non linear Euler system with various scheme …")
    divMat_upwind=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)
    divMat_staggered=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)
    divMat_centered=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    # Picture settings
    fig, ([axDensity, axMomentum],[axVelocity, axPressure]) = plt.subplots(2, 2,sharex=True, figsize=(10,10))
    fig.suptitle('Comparison of finite volume schemes')
    lineDensity_staggered=, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field_staggered, label='conservative staggered') #new picture for video # Returns a tuple of line objects, thus the comma
    lineDensity_upwind=, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field_upwind, label='upwind') #new picture for video # Returns a tuple of line objects, thus the comma
    lineDensity_centered=, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field_centered, label='centered') #new picture for video # Returns a tuple of line objects, thus the comma
    axDensity.set(xlabel='x (m)', ylabel='Density')
    axDensity.set_xlim(a,b)
    axDensity.set_ylim(min_initial_rho - 0.1*(max_initial_rho-min_initial_rho), max_initial_rho +  0.1*(max_initial_rho-min_initial_rho) )
    axDensity.legend()
    lineMomentum_staggered=, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field_staggered, label='conservative staggered')
    lineMomentum_upwind=, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field_upwind, label='upwind')
    lineMomentum_centered=, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field_centered, label='centered')
    axMomentum.set(xlabel='x (m)', ylabel='Momentum')
    axMomentum.set_xlim(a,b)
    axMomentum.set_ylim(min_initial_q - 0.1*(max_initial_q-min_initial_q), max_initial_q +  0.1*(max_initial_q-min_initial_q) )
    axMomentum.legend()
    lineVelocity_staggered=, = axVelocity.plot([a+0.5*dx + i*dx for i in range(nx)], v_field_staggered, label='conservative staggered')
    lineVelocity_upwind=, = axVelocity.plot([a+0.5*dx + i*dx for i in range(nx)], v_field_upwind, label='upwind')
    lineVelocity_centered=, = axVelocity.plot([a+0.5*dx + i*dx for i in range(nx)], v_field_centered, label='centered')
    axVelocity.set(xlabel='x (m)', ylabel='Velocity')
    axVelocity.set_xlim(a,b)
    axVelocity.set_ylim(min_initial_v - 0.4*abs(min_initial_v), max_initial_v +  0.05*abs(max_initial_v) )
    axVelocity.legend()
    linePressure_staggered=, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field_staggered, label='conservative staggered')
    linePressure_upwind=, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field_upwind, label='upwind')
    linePressure_centered=, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field_centered, label='centered')
    axPressure.set(xlabel='x (m)', ylabel='Pressure')
    axPressure.set_xlim(a,b)
    axPressure.set_ylim(min_initial_p - 0.05*abs(min_initial_p), max_initial_p +  0.05*abs(max_initial_p) )
    axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axPressure.legend()
 
    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Scheme comparison for the 1D isothermal Euler System", artist = "CEA Saclay", comment="Shock tube")
    writer=FFMpegWriter(fps=10, metadata=metadata, codec='h264')
    with writer.saving(fig, "1DEuler_System_Scheme_Comparison"+".mp4", ntmax):
        writer.grab_frame()
        plt.savefig("EulerSystem"+str(dim)+"D_Scheme_Comparison"+meshName+"_0"+".png")

        # Starting time loop
        while (it<ntmax and time <= tmax and not isStationary):
            dUn_upwind = Un_upwind.deepCopy()
            Uk_upwind  = Un_upwind.deepCopy()
            dUn_staggered = Un_staggered.deepCopy()
            Uk_staggered  = Un_staggered.deepCopy()
            # UPWIND SCHEME
            residu = 1.
            k=0
            while (k<newton_max and residu > precision ):
                #DEBUT BOUCLE NEWTON
                divMat_upwind.zeroEntries()#sets the matrix coefficients to zero
                for j in range(nbCells):# 
                    if ( j==0) : 
                        rho_r = Uk_upwind[j*nbComp+0]
                        q_r   = Uk_upwind[j*nbComp+1]
                        rho_l = rho_r # Conditions de Neumann
                        q_l   =   q_r
                        Am_upwind=  (RoeMatrix(rho_l,q_l,rho_r,q_r) - upwindDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_upwind.addValue(j*nbComp,(j+1)*nbComp,Am_upwind)
                        divMat_upwind.addValue(j*nbComp,    j*nbComp,Am_upwind*(-1.))
                        dUi[0] = Uk_upwind[(j+1)*nbComp+0]-Uk_upwind[j*nbComp+0]
                        dUi[1] = Uk_upwind[(j+1)*nbComp+1]-Uk_upwind[j*nbComp+1]
                        temp = Am_upwind*dUi
                        #print("Bloc 0 matrix  :   ", Am_upwind)
                        Rhs_upwind[j*nbComp+0] = -temp[0]-(Uk_upwind[j*nbComp+0]-Un_upwind[j*nbComp+0]) 
                        Rhs_upwind[j*nbComp+1] = -temp[1]-(Uk_upwind[j*nbComp+1]-Un_upwind[j*nbComp+1]) 
                    elif ( j==nbCells-1) :
                        rho_l = Uk_upwind[j*nbComp+0]
                        q_l   = Uk_upwind[j*nbComp+1]
                        rho_r = rho_l # Conditions de Neumann
                        q_r   =   q_l
                        Ap_upwind= (RoeMatrix(rho_l,q_l,rho_r,q_r) + staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_upwind.addValue(j*nbComp,    j*nbComp,Ap_upwind)
                        divMat_upwind.addValue(j*nbComp,(j-1)*nbComp,Ap_upwind*(-1.))
                        dUi[0] = Uk_upwind[(j-1)*nbComp+0]-Uk_upwind[j*nbComp+0]
                        dUi[1] = Uk_upwind[(j-1)*nbComp+1]-Uk_upwind[j*nbComp+1]
                        temp = Ap_upwind*dUi
                        Rhs_upwind[j*nbComp+0] = temp[0]-(Uk_upwind[j*nbComp+0]-Un_upwind[j*nbComp+0]) 
                        Rhs_upwind[j*nbComp+1] = temp[1]-(Uk_upwind[j*nbComp+1]-Un_upwind[j*nbComp+1]) 
                    else :
                        rho_l = Uk_upwind[(j-1)*nbComp+0]
                        q_l   = Uk_upwind[(j-1)*nbComp+1]
                        rho_r = Uk_upwind[j*nbComp+0]
                        q_r   = Uk_upwind[j*nbComp+1]
                        Ap_upwind = (RoeMatrix(rho_l,q_l,rho_r,q_r) + staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        ###############################################################
                        rho_l = Uk_upwind[j*nbComp+0]
                        q_l   = Uk_upwind[j*nbComp+1]
                        rho_r = Uk_upwind[(j+1)*nbComp+0]
                        q_r   = Uk_upwind[(j+1)*nbComp+1]
                        Am_upwind = (RoeMatrix(rho_l,q_l,rho_r,q_r) - staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_upwind.addValue(j*nbComp,(j+1)*nbComp,Am_upwind)
                        divMat_upwind.addValue(j*nbComp,    j*nbComp,Am_upwind*(-1.))
                        divMat_upwind.addValue(j*nbComp,    j*nbComp,Ap_upwind)
                        divMat_upwind.addValue(j*nbComp,(j-1)*nbComp,Ap_upwind*(-1.))
                        dUi1[0] = Uk_upwind[(j+1)*nbComp+0]-Uk_upwind[j*nbComp+0]
                        dUi1[1] = Uk_upwind[(j+1)*nbComp+1]-Uk_upwind[j*nbComp+1]
                        dUi2[0] = Uk_upwind[(j-1)*nbComp+0]-Uk_upwind[j*nbComp+0]
                        dUi2[1] = Uk_upwind[(j-1)*nbComp+1]-Uk_upwind[j*nbComp+1] 
                        temp1 = Am_upwind*dUi1
                        temp2 = Ap_upwind*dUi2
                        Rhs_upwind[j*nbComp+0] = -temp1[0] + temp2[0] -(Uk_upwind[j*nbComp+0]-Un_upwind[j*nbComp+0]) 
                        Rhs_upwind[j*nbComp+1] = -temp1[1] + temp2[1] -(Uk_upwind[j*nbComp+1]-Un_upwind[j*nbComp+1]) 
                divMat_upwind.diagonalShift(1)#only after  filling all coefficients
                LS_upwind=cdmath.LinearSolver(divMat_upwind,Rhs_upwind,iterGMRESMax, precision, "GMRES","LU")
                dUk_upwind=LS_upwind.solve(); 
                residu = dUk_upwind.norm()
                #stop
                Uk_upwind+=dUk_upwind
                if(not LS_upwind.getStatus()):
                    print("Upwind scheme Linear system did not converge ", LS_upwind.getNumberOfIter(), " GMRES iterations")
                    raise ValueError("Upwind scheme  : Pas de convergence du système linéaire");
                k=k+1
            Un_upwind = Uk_upwind.deepCopy()
            dUn_upwind-=Un_upwind
            for k in range(nbCells):
                rho_field_upwind[k] = Un_upwind[k*(dim+1)+0]
                q_field_upwind[k]   = Un_upwind[k*(dim+1)+1] 
            v_field_upwind   = [   q_field_upwind[i]/rho_field_upwind[i]  for i in range(nx)]
            p_field_upwind   = [ rho_field_upwind[i]*(c0*c0)       for i in range(nx)]

            lineDensity_upwind.set_ydata(rho_field_upwind)
            lineMomentum_upwind.set_ydata(q_field_upwind)
            lineVelocity_upwind.set_ydata(v_field_upwind)
            linePressure_upwind.set_ydata(p_field_upwind)

            # STAGGERED SCHEME
            residu = 1.
            k=0
            while (k<newton_max and residu > precision ):
                #DEBUT BOUCLE NEWTON
                divMat_upwind.zeroEntries()#sets the matrix coefficients to zero
                for j in range(nbCells):# 
                    if ( j==0) : 
                        rho_r = Uk_staggered[j*nbComp+0]
                        q_r   = Uk_staggered[j*nbComp+1]
                        rho_l = rho_r # Conditions de Neumann
                        q_l   =   q_r
                        Am_staggered=  (RoeMatrix(rho_l,q_l,rho_r,q_r) - staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_staggered.addValue(j*nbComp,(j+1)*nbComp,Am_staggered)
                        divMat_staggered.addValue(j*nbComp,    j*nbComp,Am_staggered*(-1.))
                        dUi[0] = Uk_staggered[(j+1)*nbComp+0]-Uk_staggered[j*nbComp+0]
                        dUi[1] = Uk_staggered[(j+1)*nbComp+1]-Uk_staggered[j*nbComp+1]
                        temp = Am_staggered*dUi
                        #print("Bloc 0 matrix  :   ", Am_staggered)
                        Rhs_staggered[j*nbComp+0] = -temp[0]-(Uk_staggered[j*nbComp+0]-Un_staggered[j*nbComp+0]) 
                        Rhs_staggered[j*nbComp+1] = -temp[1]-(Uk_staggered[j*nbComp+1]-Un_staggered[j*nbComp+1]) 
                    elif ( j==nbCells-1) :
                        rho_l = Uk_staggered[j*nbComp+0]
                        q_l   = Uk_staggered[j*nbComp+1]
                        rho_r = rho_l # Conditions de Neumann
                        q_r   =   q_l
                        Ap_staggered= (RoeMatrix(rho_l,q_l,rho_r,q_r) + staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_staggered.addValue(j*nbComp,    j*nbComp,Ap_staggered)
                        divMat_staggered.addValue(j*nbComp,(j-1)*nbComp,Ap_staggered*(-1.))
                        dUi[0] = Uk_staggered[(j-1)*nbComp+0]-Uk_staggered[j*nbComp+0]
                        dUi[1] = Uk_staggered[(j-1)*nbComp+1]-Uk_staggered[j*nbComp+1]
                        temp = Ap_staggered*dUi
                        Rhs_staggered[j*nbComp+0] = temp[0]-(Uk_staggered[j*nbComp+0]-Un_staggered[j*nbComp+0]) 
                        Rhs_staggered[j*nbComp+1] = temp[1]-(Uk_staggered[j*nbComp+1]-Un_staggered[j*nbComp+1]) 
                    else :
                        rho_l = Uk_staggered[(j-1)*nbComp+0]
                        q_l   = Uk_staggered[(j-1)*nbComp+1]
                        rho_r = Uk_staggered[j*nbComp+0]
                        q_r   = Uk_staggered[j*nbComp+1]
                        Ap_staggered = (RoeMatrix(rho_l,q_l,rho_r,q_r) + staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        ###############################################################
                        rho_l = Uk_staggered[j*nbComp+0]
                        q_l   = Uk_staggered[j*nbComp+1]
                        rho_r = Uk_staggered[(j+1)*nbComp+0]
                        q_r   = Uk_staggered[(j+1)*nbComp+1]
                        Am_staggered = (RoeMatrix(rho_l,q_l,rho_r,q_r) - staggeredDMatrix(rho_l,q_l,rho_r,q_r))*(0.5*dt/dx)
                        divMat_staggered.addValue(j*nbComp,(j+1)*nbComp,Am_staggered)
                        divMat_staggered.addValue(j*nbComp,    j*nbComp,Am_staggered*(-1.))
                        divMat_staggered.addValue(j*nbComp,    j*nbComp,Ap_staggered)
                        divMat_staggered.addValue(j*nbComp,(j-1)*nbComp,Ap_staggered*(-1.))
                        dUi1[0] = Uk_staggered[(j+1)*nbComp+0]-Uk_staggered[j*nbComp+0]
                        dUi1[1] = Uk_staggered[(j+1)*nbComp+1]-Uk_staggered[j*nbComp+1]
                        dUi2[0] = Uk_staggered[(j-1)*nbComp+0]-Uk_staggered[j*nbComp+0]
                        dUi2[1] = Uk_staggered[(j-1)*nbComp+1]-Uk_staggered[j*nbComp+1] 
                        temp1 = Am_staggered*dUi1
                        temp2 = Ap_staggered*dUi2
                        Rhs_staggered[j*nbComp+0] = -temp1[0] + temp2[0] -(Uk_staggered[j*nbComp+0]-Un_staggered[j*nbComp+0]) 
                        Rhs_staggered[j*nbComp+1] = -temp1[1] + temp2[1] -(Uk_staggered[j*nbComp+1]-Un_staggered[j*nbComp+1]) 
                divMat_staggered.diagonalShift(1)#only after  filling all coefficients
                LS_staggered=cdmath.LinearSolver(divMat_staggered,Rhs_staggered,iterGMRESMax, precision, "GMRES","LU")
                dUk_staggered=LS_staggered.solve(); 
                residu = dUk_staggered.norm()
                #stop
                Uk_staggered+=dUk_staggered
                if(not LS_staggered.getStatus()):
                    print("Staggered scheme Linear system did not converge ", LS_staggered.getNumberOfIter(), " GMRES iterations")
                    raise ValueError("Staggered scheme : Pas de convergence du système linéaire");
                k=k+1
            Un_staggered = Uk_staggered.deepCopy()
            dUn_staggered-=Un_staggered
            for k in range(nbCells):
                rho_field_staggered[k] = Un_staggered[k*(dim+1)+0]
                q_field_staggered[k]   = Un_staggered[k*(dim+1)+1] 
            v_field_staggered   = [   q_field_staggered[i]/rho_field_staggered[i]  for i in range(nx)]
            p_field_staggered   = [ rho_field_staggered[i]*(c0*c0)       for i in range(nx)]

            lineDensity_staggered.set_ydata(rho_field_staggered)
            lineMomentum_staggered.set_ydata(q_field_staggered)
            lineVelocity_staggered.set_ydata(v_field_staggered)
            linePressure_staggered.set_ydata(p_field_staggered)

            writer.grab_frame()

            time=time+dt;
            it=it+1;
            #Sauvegardes
            if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                #print("Upwind : Last linear system converged in ", LS_upwind.getNumberOfIter(), " GMRES iterations", ", residu final:   ",residu)

                plt.savefig("EulerSystem"+str(dim)+"D_Scheme_Comparison"+meshName+"_"+str(it)+".png")
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
    EulerSystemSchemeComparison(ntmax, tmax, cfl, a,b,nx, output_freq, meshName)
    return

if __name__ == """__main__""":
    a=0.
    b=1.
    nx=100
    cfl=0.99
    solve( a,b,nx,"SquareRegularSquares","RegularSquares",cfl)
