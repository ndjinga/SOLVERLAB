#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du Burgers 1D \partial_t u + 1/2 \partial_x u^2 = 0 
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Comparaison du schéma de Godunov explicite avec des schémas alternatifs sur un maillage 1D régulier
#               Premier schéma alternatif = Godunov sur l'équation conservative \partial_t u^3 + 1/3 \partial_x u^3 = 0 
#               Deuxième schéma alternatif = Upwind sur l'équation non conservative \partial_t u + u \partial_x u = 0 
#               Donnée initiale continue générant une onde de choc
#               Conditions aux limites de Neumann
#               Conclusion : 1) le choc ne peut être capturé qu'avec un schéma conservatif 2) la vitesse du choc dépend de la formulation conservative employée
#		        Création et sauvegarde du champ résultant et des figures
#               Génération d'une video sauvegardée dans un fichier .mp4
#================================================================================================================================


from math import sin, cos, pi, sqrt
import numpy as np
from copy import deepcopy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

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

def Flux_Godunov2(u_l, u_r):
    if (u_l==u_r):
        flux = 1./3.*u_l*u_l*u_l
    elif (u_l*u_l<u_r*u_r):
        flux = min(1./3.*u_l*u_l*u_l,1./3.*u_r*u_r*u_r)
    elif (u_l*u_l>u_r*u_r):
        flux = max(1./3.*u_l*u_l*u_l,1./3.*u_r*u_r*u_r)
    else:
        print( "u_l ",u_l, " u_r ", u_r )
    return flux

def Du_ncsv(u_l, u_i, u_r):
    if (u_i<0):
        Du= u_r-u_i
    elif (0<u_i):
        Du= u_i-u_l
    else:
        Du= u_r-u_l/2
    
    return Du

def Burgers1D():
    print( "Simulation of 1D Burgers' equation with various explicit schemes" )

    ##################### Simulation parameters
    a = 0.0 # space domain :  a <= x <= b
    b = 1.
    nx=100
    dx = (b - a) / nx #space step

    tmax = 1.5 # runs the simulation for 0 <= t <= tMax
    ntmax=100
    cfl=0.95

    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    ########################## Initial data
    xa=0.1
    xb=0.5
    u_initial = [ (xi<xa)+(xa<=xi)*(xi<=xb)*(np.cos(np.pi*(xa-xi)/(xa-xb))+1.0)*0.5  for xi in x];# to be used with a=0, b=1
    u_godunov = [ (xi<xa)+(xa<=xi)*(xi<=xb)*(np.cos(np.pi*(xa-xi)/(xa-xb))+1.0)*0.5  for xi in x];# to be used with a=0, b=1
    u_ncsv    = [ (xi<xa)+(xa<=xi)*(xi<=xb)*(np.cos(np.pi*(xa-xi)/(xa-xb))+1.0)*0.5  for xi in x];# to be used with a=0, b=1
    u_csv2    = deepcopy(u_initial)# to be used with a=0, b=1
    Unp1      = [0]*nx
    Unp1_ncsv = [0]*nx
    Unp1_csv2 = [0]*nx
    
    max_initial=max(u_initial)
    min_initial=min(u_initial)

    time = 0.
    it = 0
    output_freq = 10

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Finite volumes schemes for Burgers equation", artist = "CEA Saclay", comment="Shock formation")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DBurgersEquation_FV"+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.3*(max_initial-min_initial) )
        plt.title("Finite volume schemes for Burgers' equation on " + str(nx) + " cells")
        line1, = plt.plot(x, u_godunov, label='Conservative (Godunov) scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line2, = plt.plot(x, u_ncsv,    label='Non conservative scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        line3, = plt.plot(x, u_csv2,    label='Alternative conservative scheme') #new picture for video # Returns a tuple of line objects, thus the comma
        plt.legend()
    
        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) )
        np.savetxt("BurgersEquation_FVGodunov_ResultField_0"+".txt", u_godunov, delimiter="\n")
        np.savetxt("BurgersEquation_FVNonCons_ResultField_0"+".txt", u_ncsv, delimiter="\n")
        np.savetxt("BurgersEquation_FVAltCons_ResultField_0"+".txt", u_csv2, delimiter="\n")
        writer.grab_frame()
        plt.savefig("BurgersEquation_FV_Shock_ResultField_0"+".png")

        ############################# Time loop
        while (it < ntmax and time <= tmax):
            upw = max(np.abs(u_godunov))
            dt = cfl * dx / upw
            # Loop on all cells
            for i in range(0,nx):
                if (i==0):
                    flux_iminus = 0.5*u_godunov[0]*u_godunov[0]#Flux at the left Neumann boundary
                if (i==nx-1):
                    flux_iplus  = 0.5*u_godunov[nx-1]*u_godunov[nx-1]#Flux at the right Neumann boundary
                else:
                    flux_iplus = Flux_Godunov(u_godunov[i],u_godunov[i+1])
                pass
                Unp1[i] = u_godunov[i] - dt/dx*(flux_iplus-flux_iminus);
                flux_iminus = flux_iplus;
            u_godunov = Unp1
    
            for i in range(0,nx):
                if (i==0):
                    Du = Du_ncsv(u_ncsv[0],    u_ncsv[0],    u_ncsv[1])# Neumann boundary condition
                elif (i==nx-1):
                    Du = Du_ncsv(u_ncsv[nx-2], u_ncsv[nx-1], u_ncsv[nx-1])# Neumann boundary condition
                else:
                    Du = Du_ncsv(u_ncsv[i-1],    u_ncsv[i],    u_ncsv[i+1]) 
                pass
                Unp1_ncsv[i] = u_ncsv[i] - dt/dx*u_ncsv[i]*Du
            u_ncsv = Unp1_ncsv

            for i in range(0,nx):
                if (i==0):
                    flux_iminus = 1./3.*u_csv2[0]*u_csv2[0]*u_csv2[0]#Flux at the left Neumann boundary
                if (i==nx-1):
                    flux_iplus  = 1./3.*u_csv2[nx-1]*u_csv2[nx-1]*u_csv2[nx-1]#Flux at the right Neumann boundary
                else:
                    flux_iplus = Flux_Godunov2(u_csv2[i],u_csv2[i+1])
                pass
                Unp1_csv2[i] = sqrt(u_csv2[i]*u_csv2[i] - dt/dx*(flux_iplus-flux_iminus))
                flux_iminus = flux_iplus;
            u_csv2 = Unp1_csv2
    
            time += dt
            it += 1

            # Postprocessing
            line1.set_ydata(u_godunov)
            line2.set_ydata(u_ncsv)
            line3.set_ydata(u_csv2)
            writer.grab_frame()
            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "BurgersEquation_FVGodunov_ResultField_"+str(it)+".txt", u_godunov, delimiter="\n")
                np.savetxt( "BurgersEquation_FVNonCons_ResultField_"+str(it)+".txt", u_ncsv,    delimiter="\n")
                np.savetxt( "BurgersEquation_FVAltCons_ResultField_"+str(it)+".txt", u_csv2,    delimiter="\n")
                plt.savefig("BurgersEquation_FV_Shock_ResultField_"+str(it)+".png")
                #plt.show()
    
    print("Simulation of shock formation on Burgers' equation with finite volume schemes done.")
    
    return

if __name__ == """__main__""":
    Burgers1D()
