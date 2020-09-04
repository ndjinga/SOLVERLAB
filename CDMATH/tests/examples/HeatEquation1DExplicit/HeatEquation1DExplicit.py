#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation de la chaleur 1D \partial_t u = d \partial_xx u avec conditions aux limites périodiques
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Maillage 1D régulier
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

def HeatEquation1DExplicit(nx,cfl):
    print( "Simulation of 1D heat equation with an explicit scheme")

    ##################### Simulation parameters
    a = 0.0 # space domain :  a <= x <= b
    b = 1.0
    dx = (b - a) / nx #space step

    d = 1. # thermal diffusivity
    tmax = (b-a)/d # runs the simulation for 0 <= t <= tMax
    dt = cfl * dx *dx/ (2*d)
    ntmax = 100

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

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Explicit scheme for heat equation", artist = "CEA Saclay", comment="CFL="+str(cfl)+", Stable if CFL<1")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "HeatEquation1D_Explicit_"+str(nx)+"Cells_CFL"+str(cfl)+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.1*(max_initial-min_initial) )
        plt.title("Explicit scheme for the heat equation, CFL="+str(cfl))
        line1, = plt.plot(x, u, label='u') #new picture for video # Returns a tuple of line objects, thus the comma
    
        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
        np.savetxt( "HeatEquation1D_Explicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_0.txt", u, delimiter="\n")
        writer.grab_frame()
        plt.savefig("HeatEquation1D_Explicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_0.png")

        ############################# Time loop
        while (it < ntmax and time <= tmax):
            un=deepcopy(u)
            for i in range(nx):
                u[i] = un[i] + d * dt / (dx*dx) * (un[(i+1)%nx] -2*un[i] + un[(i-1)%nx])/2
    
            time += dt
            it += 1

            if cfl<1 :
                assert max(u) <= max_initial
                assert min(u) >= min_initial
                assert np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) <= total_var_initial

            # Postprocessing
            line1.set_ydata(u)
            writer.grab_frame()
            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "HeatEquation1D_Explicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_"+str(it)+".txt", u, delimiter="\n")
                plt.savefig("HeatEquation1D_Explicit_"+str(nx)+"Cells_CFL"+str(cfl)+"_ResultField_"+str(it)+".png")
                #plt.show()

    print( "Exact solution minimum   : ", min(u_initial), "Numerical solution minimum   : ",  min(u) )
    print( "Exact solution maximum   : ", max(u_initial), "Numerical solution maximum   : ",  max(u) )
    print( "Exact solution variation : ", np.sum([abs(u_initial[i] - u_initial[(i-1)%nx]) for i in range(nx)]), "Numerical solution variation : ",  np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) )
    print( "l1 numerical error       : ", dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)])        )
    
    print("Simulation of heat equation with explicit scheme done.")
    
    #return min, max, total variation and l1 error
    return min(u), max(u), np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]), dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)])


if __name__ == """__main__""":
    if len(sys.argv) >2 :
        nx = int(sys.argv[1])
        cfl = float(sys.argv[2])
        HeatEquation1DExplicit(nx,cfl)
    else :
        nx = 50 # number of cells
        cfl = 1 # c*dt/(2*dx) <= CFL
        HeatEquation1DExplicit(nx,cfl)
    
