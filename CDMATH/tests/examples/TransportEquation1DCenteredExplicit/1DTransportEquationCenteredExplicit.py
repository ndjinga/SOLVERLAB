#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du transport 1D \partial_t u + c \partial_x u = 0 avec conditions aux limites périodiques
# Author      : Michaël Ndjinga, Katia Ait Ameur
# Copyright   : CEA Saclay 2018
# Description : Utilisation du schéma centré explicite sur un maillage 1D régulier
#		        Création et sauvegarde du champ résultant et des figures
#               Génération d'une video sauvegardée dans un fichier .mp4
#================================================================================================================================


from math import sin, pi, ceil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import sys
from copy import deepcopy

def Transport1DCenteredExplicit(nx,cfl, isSmooth):
    print( "Simulation of 1D transport equation with explicit centered scheme" )

    ##################### Simulation parameters
    a = 0.0 # space domain :  a <= x <= b
    b = 1.0
    dx = (b - a) / nx #space step

    c = 0.25 # advection velocity
    dt = cfl * dx / c

    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    ########################## Initial data
    
    if(isSmooth):
        print( "Smooth initial data" )
        u_initial = [ sin(2*pi*xi)  for xi in x];# to be used with a=0, b=1
        tmax = 3*(b-a)/c # runs the simulation for 0 <= t <= tMax
    else:
        print( "Stiff initial data" )
        u_initial = [ int(1./3<xi)*int(xi<2./3)  for xi in x];# to be used with a=0, b=1
        tmax = (b-a)/c # runs the simulation for 0 <= t <= tMax
    ntmax = ceil(tmax/dt)

    u = deepcopy(u_initial)
    
    max_initial=max(u_initial)
    min_initial=min(u_initial)

    time = 0.
    it = 0
    output_freq = int(10/cfl)

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Centered explicit scheme for transport equation", artist = "CEA Saclay", comment="CFL="+str(cfl)+", Stable if CFL<1")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DTransportEquation_CenteredExplicit_"+str(nx)+"Cells_Smoothness"+str(isSmooth)+"_CFL"+str(cfl)+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.5*(max_initial-min_initial), max_initial +  0.5*(max_initial-min_initial) )
        plt.title('Centered explicit scheme for transport equation')
        line1, = plt.plot(x, u, label='u') #new picture for video # Returns a tuple of line objects, thus the comma
    
        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
        np.savetxt( "TransportEquation_CenteredExplicit_"+str(nx)+"Cells_Smoothness"+str(isSmooth)+"_CFL"+str(cfl)+"_ResultField_0.txt", u, delimiter="\n")
        writer.grab_frame()
        plt.savefig("TransportEquation_CenteredExplicit_"+str(nx)+"Cells_Smoothness"+str(isSmooth)+"_CFL"+str(cfl)+"_ResultField_0.png")

        ############################# Time loop
        while (it < ntmax and time <= tmax):
            un=deepcopy(u)
            for i in range(nx):
                u[i] = un[i] - c * dt / dx * (un[(i+1)%nx] - un[(i-1)%nx])/2
    
            time += dt
            it += 1

            # Postprocessing
            line1.set_ydata(u)
            writer.grab_frame()
            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "TransportEquation_CenteredExplicit_"+str(nx)+"Cells_Smoothness"+str(isSmooth)+"_CFL"+str(cfl)+"_ResultField_"+str(it)+".txt", u, delimiter="\n")
                plt.savefig("TransportEquation_CenteredExplicit_"+str(nx)+"Cells_Smoothness"+str(isSmooth)+"_CFL"+str(cfl)+"_ResultField_"+str(it)+".png")
                #plt.show()

    print( "Exact solution minimum   : ", min(u_initial), "Numerical solution minimum   : ",  min(u) )
    print( "Exact solution maximum   : ", max(u_initial), "Numerical solution maximum   : ",  max(u) )
    print( "Exact solution variation : ", np.sum([abs(u_initial[i] - u_initial[(i-1)%nx]) for i in range(nx)]), "Numerical solution variation : ",  np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) )
    print( "l1 numerical error       : ", dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)])        
    
    print( "Simulation of 1D transport equation with explicit centered scheme done" )
    
    #return min, max, total variation and l1 error
    return min(u), max(u), np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]), dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)])


if __name__ == """__main__""":
    if len(sys.argv) >3 :
        nx = int(sys.argv[1])
        cfl = float(sys.argv[2])
        isSmooth=bool(int(sys.argv[3]))
        Transport1DCenteredExplicit(nx,cfl,isSmooth)
    else :
        nx = 50 # number of cells
        cfl = 0.99 # c*dt/dx <= CFL
        isSmooth=True
        Transport1DCenteredExplicit(nx,cfl,isSmooth)
    
