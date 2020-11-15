# Isothermal Euler system
# d rho/d t + d q/d x =0
# d q/d t + d (q^2/rho+p)/d x = 0, where q = rho*u and p = c^2*rho
# UU = (rho,q) : conservative variable
# Scheme : Roe scheme (without any entropy correction)
# Comment : the solution displays a non entropic (non physical) shock instead of a rarefaction wave
# Date : November 2020

from cdmath import *
import numpy as np
from math import sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

c0=330.#reference sound speed for water at 15 bars
precision=1e-5

def flux(rho,q):
    fl = Vector(2);
    fl[0] = q;
    fl[1] = q*q/rho + c0*c0*rho;    
    return fl

def compTimeStep(UU,dx,CFL):
    M = UU.getMesh();
    maxAbsEigVa = 0;
    for i in range(M.getNumberOfCells()):
        maxAbsEigVa = max(maxAbsEigVa,abs(UU[i,1]/UU[i,0]+c0),abs(UU[i,1]/UU[i,0]-c0));
        pass
    dt = CFL*dx/maxAbsEigVa;
    return dt

def AbsRoeMatrix(rho_l,q_l,rho_r,q_r):
    AbsRoeMa = Matrix(2,2);
    u_l = q_l/rho_l;
    u_r = q_r/rho_r;
    u = (u_l*sqrt(rho_l)+u_r*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));   
    AbsRoeMa[0,0]=(abs(u-c0)*(u+c0)+abs(u+c0)*(c0-u))/(2*c0);
    AbsRoeMa[0,1]= (abs(u+c0)-abs(u-c0))/(2*c0);
    AbsRoeMa[1,0]=(abs(u-c0)*(u*u-c0*c0)+abs(u+c0)*(c0*c0-u*u))/(2*c0);
    AbsRoeMa[1,1]=((u+c0)*abs(u+c0)-(u-c0)*abs(u-c0))/(2*c0);
    return AbsRoeMa;
    
def main(meshName):
    # mesh definition    
    xmin=.0;
    xmax= 1;
    nx=100#1000;
    ntmax=50;
    dx = (xmax-xmin)/nx;
    M=Mesh(xmin,xmax,nx);
    dim=1
    
    # initial time
    time=0.;
    # maximum time T
    tmax=0.01;
    it=0;
    freqSortie=10;
    # CFL condition
    CFL=0.99;
    # conservative variable (unknown)
    UU=Field("Conservative variables",CELLS,M,2);

    # initial condition for Riemann problem
    print("Construction of the initial condition …")
    rhoL=1
    vL=-300
    rhoR=1.45
    vR=-300
    for i in range(M.getNumberOfCells()):
        x=M.getCell(i).x();
        if x < (xmax-xmin)/2:
            UU[i,0] = rhoL;
            UU[i,1] = rhoL*vL;
        else:
            UU[i,0] = rhoR;
            UU[i,1] = rhoR*vR;
        pass

    rho_field = [ UU[i,0]          for i in range(nx)]
    q_field   = [ UU[i,1]          for i in range(nx)]
    v_field   = [ UU[i,1]/UU[i,0]  for i in range(nx)]
    p_field   = [ UU[i,0]*(c0*c0)  for i in range(nx)]

    max_initial_rho=max(rhoL,rhoR)
    min_initial_rho=min(rhoL,rhoR)
    max_initial_p=max_initial_rho*c0*c0
    min_initial_p=min_initial_rho*c0*c0
    max_initial_v=max(vL,vR)
    min_initial_v=min(vL,vR)
    max_initial_q=max(vL*rhoL,vR*rhoR)
    min_initial_q=min(vL*rhoL,vR*rhoR)

    print( "Numerical solution of the 1D Euler equation")

    # prepare some memory
    UU_limL = Vector(2);
    UU_limR = Vector(2);
    Del_UU_RL = Vector(2);
    UU_new = UU;
    
    # Picture settings
    fig, ([axDensity, axMomentum],[axVelocity, axPressure]) = plt.subplots(2, 2,sharex=True, figsize=(10,10))
    fig.suptitle('Explicit Upwind scheme')
    lineDensity, = axDensity.plot([xmin+0.5*dx + i*dx for i in range(nx)], rho_field, label='Density') #new picture for video # Returns a tuple of line objects, thus the comma
    axDensity.set(xlabel='x (m)', ylabel='Density')
    axDensity.set_xlim(xmin,xmax)
    axDensity.set_ylim(min_initial_rho - 0.1*(max_initial_rho-min_initial_rho), max_initial_rho +  0.1*(max_initial_rho-min_initial_rho) )
    axDensity.legend()
    lineMomentum, = axMomentum.plot([xmin+0.5*dx + i*dx for i in range(nx)], q_field, label='Momentum')
    axMomentum.set(xlabel='x (m)', ylabel='Momentum')
    axMomentum.set_xlim(xmin,xmax)
    axMomentum.set_ylim(min_initial_q - 0.1*(max_initial_q-min_initial_q), max_initial_q +  0.1*(max_initial_q-min_initial_q) )
    axMomentum.legend()
    lineVelocity, = axVelocity.plot([xmin+0.5*dx + i*dx for i in range(nx)], v_field, label='Velocity')
    axVelocity.set(xlabel='x (m)', ylabel='Velocity')
    axVelocity.set_xlim(xmin,xmax)
    axVelocity.set_ylim(min_initial_v - 0.25*abs(min_initial_v), max_initial_v +  0.05*abs(max_initial_v) )
    axVelocity.legend()
    linePressure, = axPressure.plot([xmin+0.5*dx + i*dx for i in range(nx)], p_field, label='Pressure')
    axPressure.set(xlabel='x (m)', ylabel='Pressure')
    axPressure.set_xlim(xmin,xmax)
    axPressure.set_ylim(min_initial_p - 0.05*abs(min_initial_p), max_initial_p +  0.05*abs(max_initial_p) )
    axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axPressure.legend()
 
    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Upwind (Roe) scheme for the 1D isothermal Euler System", artist = "CEA Saclay", comment="Shock tube")
    writer=FFMpegWriter(fps=10, metadata=metadata, codec='h264')
    with writer.saving(fig, "1DEuler_System_Upwind"+".mp4", ntmax):
        writer.grab_frame()
        plt.savefig("EulerSystem"+str(dim)+"UpwindExplicit"+meshName+"_0"+".png")

        # loop in time
        while (it<ntmax and time <= tmax ):
            # time step
            dt=compTimeStep(UU,dx,CFL);
            # Neumann boundary condition
            UU_limL[0] = UU[0,0];
            UU_limL[1] = UU[0,1];
            UU_limR[0] = UU[M.getNumberOfCells()-1,0];
            UU_limR[1] = UU[M.getNumberOfCells()-1,1];
            flux_iminus = flux(UU_limL[0],UU_limL[1]);
            #Main loop on cells
            for i in range(0,M.getNumberOfCells()):
                if (i<M.getNumberOfCells()-1):
                    Del_UU_RL[0]=UU[i+1,0]-UU[i,0];
                    Del_UU_RL[1]=UU[i+1,1]-UU[i,1];
                    AbsRoeMa = AbsRoeMatrix(UU[i,0],UU[i,1],UU[i+1,0],UU[i+1,1]);
                    flux_iplus = (flux(UU[i,0],UU[i,1])+flux(UU[i+1,0],UU[i+1,1]))-AbsRoeMa*Del_UU_RL;
                    flux_iplus[0]*=0.5;
                    flux_iplus[1]*=0.5;
                else:
                    flux_iplus= flux(UU_limR[0],UU_limR[1]);
                UU_new[i,0] = UU[i,0] - dt/dx*(flux_iplus[0]-flux_iminus[0]);
                UU_new[i,1] = UU[i,1] - dt/dx*(flux_iplus[1]-flux_iminus[1]);
                flux_iminus = flux_iplus;
                pass
            UU = UU_new;
            time+=dt;
            it+=1;

            rho_field = [ UU[i,0]  for i in range(nx)]
            q_field   = [ UU[i,1]  for i in range(nx)]
            v_field   = [ UU[i,1]/UU[i,0]  for i in range(nx)]
            p_field   = [ UU[i,0]*(c0*c0)  for i in range(nx)]

            lineDensity.set_ydata(rho_field)
            lineMomentum.set_ydata(q_field)
            lineVelocity.set_ydata(v_field)
            linePressure.set_ydata(p_field)
            writer.grab_frame()

            if (it%freqSortie==0):
                print( "-- Iter : ", it," Time : ",time," dt : ",dt)
                plt.savefig("EulerSystem"+str(dim)+"UpwindExplicit"+meshName+"_"+str(it)+".png")
            pass

    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    if(it>=ntmax):
        print("Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
        return
    elif(isStationary):
        print("Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print("------------------------------------------------------------------------------------")
        plt.savefig("EulerSystem"+str(dim)+"UpwindExplicit"+meshName+"_Stat.png")
        return
    else:
        print("Temps maximum Tmax= ", tmax, " atteint")

    return

if __name__ == '__main__':
    main("SquareRegularSquares")
