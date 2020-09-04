#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du Burgers 1D \partial_t u + 1/2 \partial_x u^2 = 0 
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Comparaison des schémas Upwind et de Godunov explicite sur un maillage 1D régulier
#               Donnée initiale discontinue générant une onde de raréfaction
#               Conditions aux limites de Neumann
#               Conclusion : le schéma Upwind n'est pas entropique, contrairement au schéma de Godunov
#		        Création et sauvegarde du champ résultant et des figures
#               Génération d'une video sauvegardée dans un fichier .mp4
#================================================================================================================================


import numpy as np
from copy import deepcopy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

def Flux_upwind(u_l, u_r):
    if ((u_l+u_r)/2>0):
        flux = 0.5*u_l*u_l
    elif ((u_l+u_r)/2<0):
        flux = 0.5*u_r*u_r
    else:
        flux = (0.5*u_l*u_l + 0.5*u_r*u_r)/2
    
    return flux

def Flux_Godunov(u_l, u_r):
    if (u_l==u_r):
        flux = 0.5*u_l*u_l
    elif (u_l<0 and 0<u_r):
        flux = 0.;
    elif (u_l<u_r):
        flux = min(0.5*u_l*u_l,0.5*u_r*u_r);
    elif (u_l>u_r):
        flux = max(0.5*u_l*u_l,0.5*u_r*u_r);
    
    return flux

def Burgers1D():
    ##################### Simulation parameters
    a = -1. # space domain :  a <= x <= b
    b = 1.
    nx=50
    dx = (b - a) / nx #space step

    tmax = 4 # runs the simulation for 0 <= t <= tMax
    ntmax=150
    cfl=0.95

    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    ########################## Initial data
    u_initial = [ (xi<(a+b)/2)*(-1) + (xi>(a+b)/2)*(1.)  for xi in x]
    u_upwind  = [ (xi<(a+b)/2)*(-1) + (xi>(a+b)/2)*(1.)  for xi in x]
    u_godunov = [ (xi<(a+b)/2)*(-1) + (xi>(a+b)/2)*(1.)  for xi in x]
    u_exact   = [ (xi<(a+b)/2)*(-1) + (xi>(a+b)/2)*(1.)  for xi in x]

    Unp1      = [0]*nx
    Unp1_godunov = [0]*nx
    
    max_initial=max(u_initial)
    min_initial=min(u_initial)

    time = 0.
    it = 0
    output_freq = 10

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Finite volumes schemes for Burgers equation", artist = "CEA Saclay", comment="Rarefaction wave")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DBurgersEquation_FV"+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.5*(max_initial-min_initial) )
        plt.title('Finite volume schemes for Burgers equation')
        line1, = plt.plot(x, u_upwind,  label='Conservative (Upwind) scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line2, = plt.plot(x, u_godunov, label='Conservative (Godunov) scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line3, = plt.plot(x, u_exact,   label='Exact solution') #new picture for video # Returns a tuple of line objects, thus the comma
        plt.legend()
    
        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) )
        np.savetxt("BurgersEquation_FV_Upwind_ResultField_0" +".txt", u_upwind,  delimiter="\n")
        np.savetxt("BurgersEquation_FV_Godunov_ResultField_0"+".txt", u_godunov, delimiter="\n")
        np.savetxt("BurgersEquation_FV_Exact_ResultField_0"+".txt",   u_exact,   delimiter="\n")
        writer.grab_frame()
        plt.savefig("BurgersEquation_FV_Rarefaction_ResultField_0"+".png")

        ############################# Time loop
        while (it < ntmax and time <= tmax):
            upw = max(np.abs(u_upwind))
            dt = cfl * dx / upw
            # Loop on all cells
            for i in range(0,nx):
                if (i==0):
                    flux_iminus = 0.5*u_upwind[0]*u_upwind[0]#Flux at the left Neumann boundary
                if (i==nx-1):
                    flux_iplus  = 0.5*u_upwind[nx-1]*u_upwind[nx-1]#Flux at the right Neumann boundary
                else:
                    flux_iplus = Flux_upwind(u_upwind[i],u_upwind[i+1])
                pass
                Unp1[i] = u_upwind[i] - dt/dx*(flux_iplus-flux_iminus);
                flux_iminus = flux_iplus;
            u_upwind = Unp1
    
            for i in range(0,nx):
                if (i==0):
                    flux_iminus = 0.5*u_godunov[0]*u_godunov[0]#Flux at the left Neumann boundary
                if (i==nx-1):
                    flux_iplus  = 0.5*u_godunov[nx-1]*u_godunov[nx-1]#Flux at the right Neumann boundary
                else:
                    flux_iplus = Flux_Godunov(u_godunov[i],u_godunov[i+1])
                pass
                Unp1_godunov[i] = u_godunov[i] - dt/dx*(flux_iplus-flux_iminus);
                flux_iminus = flux_iplus;
            u_godunov = Unp1_godunov
    
            time += dt
            it += 1

            u_exact = [ (xi<-time)*(-1) + (xi>time)*(1.) + (-time<=xi and xi<=time)*(xi/time) for xi in x]

            # Postprocessing
            line1.set_ydata(u_upwind)
            line2.set_ydata(u_godunov)
            line3.set_ydata(u_exact)
            writer.grab_frame()

            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "BurgersEquation_Upwind_ResultField_" +str(it)+".txt", u_upwind,  delimiter="\n")
                np.savetxt( "BurgersEquation_Godunov_ResultField_"+str(it)+".txt", u_godunov, delimiter="\n")
                np.savetxt( "BurgersEquation_Exact_ResultField_"+str(it)+".txt", u_exact, delimiter="\n")
                plt.savefig("BurgersEquation_FV_Rarefaction_ResultField_"+str(it)+".png")
                #plt.show()
    
    print("Simulation of rarefaction wave on Burgers' equation with finite volume schemes done.")
    
    return

if __name__ == """__main__""":
    Burgers1D()
