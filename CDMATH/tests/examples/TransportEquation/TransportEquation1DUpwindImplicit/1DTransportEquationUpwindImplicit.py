#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du transport 1D \partial_t u + c \partial_x u = 0 avec conditions aux limites périodiques
# Author      : Michaël Ndjinga, Katia Ait Ameur
# Copyright   : CEA Saclay 2018
# Description : Utilisation du schéma upwind implicite sur un maillage 1D régulier
#		        Création et sauvegarde du champ résultant et des figures
#               Génération d'une video sauvegardée dans un fichier .mp4
#================================================================================================================================


import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import sys
from math import sin, pi, ceil
import cdmath

precision=1.e-5

def upwindSchemeMatrix(nx,cfl):
    upwindMat=cdmath.SparseMatrixPetsc(nx,nx,2)
    for i in range(nx):
        upwindMat.setValue(i,i,1+cfl)
        upwindMat.setValue(i,(i-1)%nx,-cfl)

    return upwindMat
    
def Transport1DUpwindImplicit(nx,cfl):
    print( "Simulation of 1D transport equation with implicit upwind scheme" )

    ##################### Simulation parameters
    a = 0.0 # space domain :  a <= x <= b
    b = 1.0
    dx = (b - a) / nx #space step

    c = 0.25 # advection velocity
    tmax = (b-a)/c # runs the simulation for 0 <= t <= tMax
    dt = cfl * dx / c
    ntmax = ceil(tmax/dt)

    if(cfl>nx):
        raise("Impossible to run this simulation with cfl>nx. Choose another value for nx or cfl.")
        
    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    ########################## Initial data
    
    u_initial = [ 0.5*(1+sin(4*pi*xi-pi*.5))*int(xi<0.5)*int(0<xi) + int(0.6<xi)*int(xi<0.85)  for xi in x];# to be used with a=0, b=1
    u         = [ 0.5*(1+sin(4*pi*xi-pi*.5))*int(xi<0.5)*int(0<xi) + int(0.6<xi)*int(xi<0.85)  for xi in x];# to be used with a=0, b=1

    max_initial=max(u_initial)
    min_initial=min(u_initial)
    total_var_initial = np.sum([abs(u_initial[i] - u_initial[(i-1)%nx]) for i in range(nx)])

    time = 0.
    it = 0
    output_freq = 10
        
    #Linear system initialisation
    systemMat=upwindSchemeMatrix(nx,cfl)
    iterGMRESMax=50
    Un =cdmath.Vector(nx)
    for i in range(nx):
        Un[i]=u[i]
    LS=cdmath.LinearSolver(systemMat,Un,iterGMRESMax, precision, "GMRES","ILU")

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Upwind implicit scheme for transport equation, "+"CFL="+str(cfl), artist = "CEA Saclay", comment="Stable for any CFL>0")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DTransportEquation_UpwindImplicit_"+str(nx)+"CELLS_CFL"+str(cfl)+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.1*(max_initial-min_initial) )
        plt.title("Upwind implicit scheme for transport equation, "+"CFL="+str(cfl))
        line1, = plt.plot(x, u, label='u') #new picture for video # Returns a tuple of line objects, thus the comma

        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
        np.savetxt( "TransportEquation_UpwindImplicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_0.txt", u, delimiter="\n")
        writer.grab_frame()
        plt.savefig("TransportEquation_UpwindImplicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_0.png")
    
        ############################# Time loop
        while (it < ntmax and time <= tmax):
            # Solve linear system
            for i in range(nx):
                Un[i]=u[i]
            LS.setSndMember(Un)
            Un=LS.solve()
            if(not LS.getStatus()):
                print( "Linear system did not converge ", iterGMRES, " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            for i in range(nx):
                u[i]=Un[i]
    
            time += dt
            it += 1

            assert max(u) <= max_initial
            assert min(u) >= min_initial
            assert np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) <= total_var_initial

            # Postprocessing
            line1.set_ydata(u)
            writer.grab_frame()
            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "TransportEquation_UpwindImplicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_"+str(it)+".txt", u, delimiter="\n")
                plt.savefig("TransportEquation_UpwindImplicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_"+str(it)+".png")
                #plt.show()

    print( "Exact solution minimum   : ", min(u_initial), "Numerical solution minimum   : ",  min(u) )
    print( "Exact solution maximum   : ", max(u_initial), "Numerical solution maximum   : ",  max(u) )
    print( "Exact solution variation : ", np.sum([abs(u_initial[i] - u_initial[(i-1)%nx]) for i in range(nx)]), "Numerical solution variation : ",  np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) )
    print( "l1 numerical error       : ", dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)]) )       
    
    print("Simulation of transport equation with implicit upwind scheme done.")
    
    #return min, max, total variation and l1 error
    return min(u), max(u), np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]), dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)])


if __name__ == """__main__""":
    if len(sys.argv) >2 :
        nx = int(sys.argv[1])
        cfl = float(sys.argv[2])
        Transport1DUpwindImplicit(nx,cfl)
    else :
        nx = 50 # number of cells
        cfl = 0.99 # c*dt/dx <= CFL
        Transport1DUpwindImplicit(nx,cfl)
    
