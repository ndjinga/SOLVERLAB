#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF des équations d'Euler isotherme 2D sans terme source
#                \partial_t rho + \div q = 0
#                \partial_t q   + \div q \otimes q/rho   \grad p = 0
# Author      : Michaël Ndjinga, Coraline Mounier
# Copyright   : CEA Saclay 2020
# Description : Propagation d'une onde de choc droite
#               Utilisation du schéma upwind explicite ou implicite sur un maillage général
#               Initialisation par une discontinuité verticale
#               Conditions aux limites de Neumann
#                Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


import cdmath
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from math import sqrt
import sys
import PV_routines

p0=155.e5 #reference pressure in a pressurised nuclear vessel
c0=700.   #reference sound speed for water at 155 bars
rho0=p0/(c0*c0)   #reference density
precision=1e-5

def initial_conditions_Riemann_problem(my_mesh):
    print( "Initial data : Riemann problem" )
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    xcentre = 0.5

    density_field = cdmath.Field("Density",    cdmath.CELLS, my_mesh, 1)
    q_x_field     = cdmath.Field("Momentum x", cdmath.CELLS, my_mesh, 1)
    q_y_field     = cdmath.Field("Momentum y", cdmath.CELLS, my_mesh, 1)
    #Velocity field with 3 components should be created for streamlines to be drawn
    velocity_field = cdmath.Field("Velocity", cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()

        # Initial momentum is zero
        q_x_field[i] = 0
        q_y_field[i] = 0

        # Initial velocity is zero
        velocity_field[i,0] = 0
        velocity_field[i,1] = 0
        velocity_field[i,2] = 0

        if x < xcentre:
            density_field[i] = rho0
            pass
        else:
            density_field[i] = rho0/2
            pass
        pass

    return density_field, q_x_field, q_y_field, velocity_field


def jacobianMatricesm_x(coeff,rho_l,q_lx,q_ly,rho_r,q_rx,q_ry):
    RoeMatx = cdmath.Matrix(3,3);
    Dmacx = cdmath.Matrix(3,3);

    u_lx=q_lx/rho_l
    u_rx=q_rx/rho_r
    u_ly=q_ly/rho_l
    u_ry=q_ry/rho_r
    if rho_l<0 or rho_r<0:
        print("rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    ux = (u_lx*sqrt(rho_l)+u_rx*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));
    uy = (u_ly*sqrt(rho_l)+u_ry*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));

    RoeMatx[0,0] = 0
    RoeMatx[0,1] = 1
    RoeMatx[0,2] = 0
    RoeMatx[1,0] = c0*c0 - ux*ux
    RoeMatx[1,1] = 2*ux
    RoeMatx[1,2] = 0
    RoeMatx[2,0] = -ux*uy
    RoeMatx[2,1] = uy
    RoeMatx[2,2] = ux

    Dmacx[0,0] = abs(ux)-ux
    Dmacx[0,1] = 1
    Dmacx[0,2] = 0
    Dmacx[1,0] = -c0*c0-ux*ux
    Dmacx[1,1] = abs(ux)+ux
    Dmacx[1,2] = 0
    Dmacx[2,0] = -ux*uy
    Dmacx[2,1] = uy
    Dmacx[2,2] = abs(ux)

    return (RoeMatx-Dmacx)*coeff*0.5

def jacobianMatricesp_x(coeff,rho_l,q_lx,q_ly,rho_r,q_rx,q_ry):
    RoeMatx   = cdmath.Matrix(3,3)
    Dmacx = cdmath.Matrix(3,3)

    u_lx=q_lx/rho_l
    u_rx=q_rx/rho_r
    u_ly=q_ly/rho_l
    u_ry=q_ry/rho_r
    if rho_l<0 or rho_r<0:
        print("rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    ux = (u_lx*sqrt(rho_l)+u_rx*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));
    uy = (u_ly*sqrt(rho_l)+u_ry*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));

    RoeMatx[0,0] = 0
    RoeMatx[0,1] = 1
    RoeMatx[0,2] = 0
    RoeMatx[1,0] = c0*c0 - ux*ux
    RoeMatx[1,1] = 2*ux
    RoeMatx[1,2] = 0
    RoeMatx[2,0] = -ux*uy
    RoeMatx[2,1] = uy
    RoeMatx[2,2] = ux

    Dmacx[0,0] = abs(ux)-ux
    Dmacx[0,1] = 1
    Dmacx[0,2] = 0
    Dmacx[1,0] = -c0*c0-ux*ux
    Dmacx[1,1] = abs(ux)+ux
    Dmacx[1,2] = 0
    Dmacx[2,0] = -ux*uy
    Dmacx[2,1] = uy
    Dmacx[2,2] = abs(ux)

    return (RoeMatx+Dmacx)*coeff*0.5

def jacobianMatricesm_y(coeff,rho_l,q_lx,q_ly,rho_r,q_rx,q_ry):
    RoeMaty   = cdmath.Matrix(3,3);
    Dmacy = cdmath.Matrix(3,3);

    u_lx=q_lx/rho_l
    u_rx=q_rx/rho_r
    u_ly=q_ly/rho_l
    u_ry=q_ry/rho_r
    if rho_l<0 or rho_r<0:
        print("rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    ux = (u_lx*sqrt(rho_l)+u_rx*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));
    uy = (u_ly*sqrt(rho_l)+u_ry*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));

    RoeMaty[0,0] = 0
    RoeMaty[0,1] = 0
    RoeMaty[0,2] = 1
    RoeMaty[1,0] = -ux*uy
    RoeMaty[1,1] = uy
    RoeMaty[1,2] = ux
    RoeMaty[2,0] = c0*c0-uy*uy
    RoeMaty[2,1] = 0
    RoeMaty[2,2] = 2*uy

    Dmacy[0,0] = abs(uy)-uy
    Dmacy[0,1] = 0
    Dmacy[0,2] = 1
    Dmacy[1,0] = -ux*uy
    Dmacy[1,1] = abs(uy)
    Dmacy[1,2] = ux
    Dmacy[2,0] = -c0*c0-uy*uy
    Dmacy[2,1] = 0
    Dmacy[2,2] = abs(uy)+uy

    return (RoeMaty-Dmacy)*coeff*0.5

def jacobianMatricesp_y(coeff,rho_l,q_lx,q_ly,rho_r,q_rx,q_ry):
    RoeMaty   = cdmath.Matrix(3,3);
    Dmacy = cdmath.Matrix(3,3);

    u_lx=q_lx/rho_l
    u_rx=q_rx/rho_r
    u_ly=q_ly/rho_l
    u_ry=q_ry/rho_r
    if rho_l<0 or rho_r<0:
        print("rho_l=",rho_l, " rho_r= ",rho_r)
        raise ValueError("Negative density")
    ux = (u_lx*sqrt(rho_l)+u_rx*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));
    uy = (u_ly*sqrt(rho_l)+u_ry*sqrt(rho_r))/(sqrt(rho_l)+sqrt(rho_r));

    RoeMaty[0,0] = 0
    RoeMaty[0,1] = 0
    RoeMaty[0,2] = 1
    RoeMaty[1,0] = -ux*uy
    RoeMaty[1,1] = uy
    RoeMaty[1,2] = ux
    RoeMaty[2,0] = c0*c0-uy*uy
    RoeMaty[2,1] = 0
    RoeMaty[2,2] = 2*uy

    Dmacy[0,0] = abs(uy)-uy
    Dmacy[0,1] = 0
    Dmacy[0,2] = 1
    Dmacy[1,0] = -ux*uy
    Dmacy[1,1] = abs(uy)
    Dmacy[1,2] = ux
    Dmacy[2,0] = -c0*c0-uy*uy
    Dmacy[2,1] = 0
    Dmacy[2,2] = abs(uy)+uy

    return (RoeMaty+Dmacy)*coeff*0.5

def EulerSystemStaggered(ntmax, tmax, cfl,output_freq, my_mesh, meshName):

    if not my_mesh.isStructured() :
        raise ValueError("Euler_ConservativeStaggered2D_RiemannProblem.py expects a structured mesh")

    dim=my_mesh.getMeshDimension()
    if dim != 2 :
        raise ValueError("Euler_ConservativeStaggered2D_RiemannProblem.py expects a 2D mesh")

    nbComp=dim+1
    # Mesh parameters
    nbCells = my_mesh.getNumberOfCells()
    nbCells_x = my_mesh.getNx()
    nbCells_y = my_mesh.getNy()
    dx,dy = my_mesh.getDXYZ()
    nbVoisinsMax=my_mesh.getMaxNbNeighbours(cdmath.CELLS)
    dx_min=my_mesh.minRatioVolSurf()
    
    # Time evolution parameters
    time = 0.
    it=0;
    isStationary=False
    iterGMRESMax = 50
    newton_max   = 50

    #iteration vectors
    Un =cdmath.Vector(nbCells*nbComp)
    dUn=cdmath.Vector(nbCells*nbComp)
    dUk=cdmath.Vector(nbCells*nbComp)
    Rhs=cdmath.Vector(nbCells*nbComp)

    dUi_x=cdmath.Vector(nbComp)
    dUi_y=cdmath.Vector(nbComp)
    dUi1_x=cdmath.Vector(nbComp)
    dUi2_x=cdmath.Vector(nbComp)
    dUi1_y=cdmath.Vector(nbComp)
    dUi2_y=cdmath.Vector(nbComp)

    # Initial conditions #
    print("Construction of the initial condition …")
    rho_field, q_x_field, q_y_field, velocity_field= initial_conditions_Riemann_problem(my_mesh)

    for i in range(nbCells):
        Un[nbComp*i]   = rho_field[i]
        Un[nbComp*i+1] = q_x_field[i]
        Un[nbComp*i+2] = q_y_field[i]

    #Sauvegarde de la donnée initiale
    rho_field.setTime(time,it);
    rho_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_density");
    q_x_field.setTime(time,it);
    q_x_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumX");
    q_y_field.setTime(time,it);
    q_y_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumY");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_velocity");

    print("Starting computation of the isothermal Euler system with a conservative staggered scheme …")
    divMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        dUn = Un.deepCopy()
        Uk  = Un.deepCopy()
        residu = 1e10
        k=0
        while (k<newton_max and residu > 1/precision ):

            dt = cfl * dx_min / c0# This choice should be improved when possible by using the actual eigenvalue abs(u)+c0, that is the time step should be determined avec the computation of the jacobian matrices
            #DEBUT BOUCLE NEWTON
            for j in range(nbCells_y):
                for i in range(nbCells_x):
                    #Traitement des coins
                    if ( j==0 and i==0) :
                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Am_x*(-1.)-Am_y)

                        dUi_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]

                        temp_x = Am_x*dUi_x
                        temp_y = Am_y*dUi_y
                        #print("Bloc 0 matrix  :   ", Am)
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0] - temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1] - temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2] - temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif(i==0 and j==nbCells_y-1):
                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Am_x*(-1.)+Ap_y)

                        dUi_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp_x = Am_x*dUi_x
                        temp_y = Ap_y*dUi_y
                        #print("Bloc 0 matrix  :   ", Am)
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0]- temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1]- temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2]- temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif(i==nbCells_x-1 and j==0):

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_x-Am_y)

                        dUi_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]
                        dUi_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]

                        temp_x = Ap_x*dUi_x
                        temp_y= Am_y*dUi_y
                        #print("Bloc 0 matrix  :   ", Am)
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0]- temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1]- temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2]- temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif ( j==nbCells_y-1 and i==nbCells_x-1) :

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_x+Ap_y)

                        dUi_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]
                        dUi_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp_x = Ap_x*dUi_x
                        temp_y = Ap_y*dUi_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0]-temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1]-temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2]-temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    #Traitement des bords restants
                    elif i==0 :

                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #remplissage de la divMat
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_y-Am_x-Am_y)

                        #remplissage du membre de droite
                        dUi_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi1_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi1_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi2_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi2_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi2_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp_x  = Am_x*dUi_x
                        temp1_y = Am_y*dUi1_y
                        temp2_y = Ap_y*dUi2_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0]-temp1_y[0]-temp2_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1]-temp1_y[1]-temp2_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2]-temp1_y[2]-temp2_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif i==nbCells_x-1:

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #remplissage de la divMat
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_y+Ap_x-Am_y)

                        #remplissage du membre de droite
                        dUi_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]
                        dUi1_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi1_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi2_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi2_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi2_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp_x  = Ap_x*dUi_x
                        temp1_y = Am_y*dUi1_y
                        temp2_y = Ap_y*dUi2_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp_x[0]-temp1_y[0]-temp2_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp_x[1]-temp1_y[1]-temp2_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp_x[2]-temp1_y[2]-temp2_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif j==0:

                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #remplissage de la divMat
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_x-Am_x-Am_y)

                        #remplissage du membre de droite
                        dUi1_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi1_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi2_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi2_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi2_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]
                        dUi_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]

                        temp1_x = Am_x*dUi1_x
                        temp2_x = Ap_x*dUi2_x
                        temp_y = Am_y*dUi_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp1_x[0]-temp2_x[0]-temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp1_x[1]-temp2_x[1]-temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp1_x[2]-temp2_x[2]-temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    elif j==nbCells_y-1:

                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #remplissage de la divMat
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_x+Ap_y-Am_x)

                        #remplissage du membre de droite
                        dUi1_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi1_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]
                        dUi2_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi2_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi2_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]
                        dUi_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp1_x = Am_x*dUi1_x
                        temp2_x = Ap_x*dUi2_x
                        temp_y = Ap_y*dUi_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp1_x[0]-temp2_x[0]-temp_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp1_x[1]-temp2_x[1]-temp_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp1_x[2]-temp2_x[2]-temp_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

                    #Traitement des autres cellules
                    else:

                        #Calcul de Am_x
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*(i+1)]
                        qx_r =Uk[3*nbCells_x*j+3*(i+1)+1]
                        qy_r =Uk[3*nbCells_x*j+3*(i+1)+2]
                        Am_x= jacobianMatricesm_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #Calcul de Ap_x
                        rho_l=Uk[3*nbCells_x*j+3*(i-1)]
                        qx_l =Uk[3*nbCells_x*j+3*(i-1)+1]
                        qy_l =Uk[3*nbCells_x*j+3*(i-1)+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_x= jacobianMatricesp_x(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Am_y
                        rho_l=Uk[3*nbCells_x*j+3*i]
                        qx_l =Uk[3*nbCells_x*j+3*i+1]
                        qy_l =Uk[3*nbCells_x*j+3*i+2]
                        rho_r=Uk[3*nbCells_x*(j+1)+3*i]
                        qx_r =Uk[3*nbCells_x*(j+1)+3*i+1]
                        qy_r =Uk[3*nbCells_x*(j+1)+3*i+2]
                        Am_y= jacobianMatricesm_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #calcul de Ap_y
                        rho_l=Uk[3*nbCells_x*(j-1)+3*i]
                        qx_l =Uk[3*nbCells_x*(j-1)+3*i+1]
                        qy_l =Uk[3*nbCells_x*(j-1)+3*i+2]
                        rho_r=Uk[3*nbCells_x*j+3*i]
                        qx_r =Uk[3*nbCells_x*j+3*i+1]
                        qy_r =Uk[3*nbCells_x*j+3*i+2]
                        Ap_y= jacobianMatricesp_y(dt/dx,rho_l,qx_l,qy_l,rho_r,qx_r,qy_r)

                        #remplissage de la divMat
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i-1),Ap_x*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j-1)+3*i,Ap_y*(-1.))
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*(i+1),Am_x)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*(j+1)+3*i,Am_y)
                        divMat.addValue(3*nbCells_x*j+3*i,3*nbCells_x*j+3*i,Ap_x+Ap_y-Am_x-Am_y)

                        #remplissage du membre de droite

                        dUi1_x[0] = Uk[3*nbCells_x*j+3*(i+1)]-Uk[3*nbCells_x*j+3*i]
                        dUi1_x[1] = Uk[3*nbCells_x*j+3*(i+1)+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_x[2] = Uk[3*nbCells_x*j+3*(i+1)+2]-Uk[3*nbCells_x*j+3*i+2]

                        dUi2_x[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*j+3*(i-1)]
                        dUi2_x[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*j+3*(i-1)+1]
                        dUi2_x[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*j+3*(i-1)+2]

                        dUi1_y[0] = Uk[3*nbCells_x*(j+1)+3*i]-Uk[3*nbCells_x*j+3*i]
                        dUi1_y[1] = Uk[3*nbCells_x*(j+1)+3*i+1]-Uk[3*nbCells_x*j+3*i+1]
                        dUi1_y[2] = Uk[3*nbCells_x*(j+1)+3*i+2]-Uk[3*nbCells_x*j+3*i+2]

                        dUi2_y[0] = Uk[3*nbCells_x*j+3*i]-Uk[3*nbCells_x*(j-1)+3*i]
                        dUi2_y[1] = Uk[3*nbCells_x*j+3*i+1]-Uk[3*nbCells_x*(j-1)+3*i+1]
                        dUi2_y[2] = Uk[3*nbCells_x*j+3*i+2]-Uk[3*nbCells_x*(j-1)+3*i+2]

                        temp1_x = Am_x*dUi1_x
                        temp2_x = Ap_x*dUi2_x
                        temp1_y = Am_y*dUi1_y
                        temp2_y = Ap_y*dUi2_y
                        Rhs[3*nbCells_x*j+3*i+0] = -temp1_x[0]-temp2_x[0]-temp1_y[0]-temp2_y[0]-(Uk[3*nbCells_x*j+3*i+0]-Un[3*nbCells_x*j+3*i+0])
                        Rhs[3*nbCells_x*j+3*i+1] = -temp1_x[1]-temp2_x[1]-temp1_y[1]-temp2_y[1]-(Uk[3*nbCells_x*j+3*i+1]-Un[3*nbCells_x*j+3*i+1])
                        Rhs[3*nbCells_x*j+3*i+2] = -temp1_x[2]-temp2_x[2]-temp1_y[2]-temp2_y[2]-(Uk[3*nbCells_x*j+3*i+2]-Un[3*nbCells_x*j+3*i+2])

            divMat.diagonalShift(1)  #only after  filling all coefficients
            LS=cdmath.LinearSolver(divMat,Rhs,iterGMRESMax, precision, "GMRES","LU")
            dUk=LS.solve();
            residu = dUk.norm()
            Uk+=dUk
            #print("Newton iteration number ",k, "residu = ",residu)
            if(not LS.getStatus()):
                print("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
                raise ValueError("Pas de convergence du système linéaire");
            k=k+1
        Un = Uk.deepCopy()
        dUn-=Un
        isStationary = dUn.norm()<precision
        
        for i in range(nbCells):
            rho_field[i] = Un[nbComp*i]
            q_x_field[i] = Un[nbComp*i+1]
            q_y_field[i] = Un[nbComp*i+2]

        time=time+dt;
        it=it+1;
        #Sauvegardes
        if(it==1 or it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt) + ", Newton iterations: "+str(k)+", ||dUn|| = "+str(dUn.norm()))
            print("Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations, résidu = ", residu)
            for k in range(nbCells):
                rho  =rho_field[k]
                velocity_field[k,0]=q_x_field[i]/rho
                if(dim>1):
                    velocity_field[k,1]=q_y_field[k]/rho
            print
            rho_field.setTime(time,it);
            rho_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_density",False);
            q_x_field.setTime(time,it);
            q_x_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumX",False);
            q_y_field.setTime(time,it);
            q_y_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumY",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_velocity",False);
            #Postprocessing : save 2D picture
            PV_routines.Save_PV_data_to_picture_file("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_density_"  +str(it)+'.vtu',"Density",   'CELLS',"EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_density"  +str(it))
            PV_routines.Save_PV_data_to_picture_file("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumX_"+str(it)+'.vtu',"Momentum x",'CELLS',"EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumX"+str(it))
            PV_routines.Save_PV_data_to_picture_file("EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumY_"+str(it)+'.vtu',"Momentum y",'CELLS',"EulerIsothermal_"+str(dim)+"DConservativeStaggered_"+meshName+"_momentumY"+str(it))
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    if(it>=ntmax):
        print("Nombre de pas de temps maximum ntmax= ", ntmax, " atteint")
        return
    elif(isStationary):
        print("Régime stationnaire atteint au pas de temps ", it, ", t= ", time)
        print("------------------------------------------------------------------------------------")
        return
    else:
        print("Temps maximum Tmax= ", tmax, " atteint")
        return

def solve( my_mesh,filename, resolution):
    print("Resolution of the Isothermal Euler system in dimension 2 on "+str(my_mesh.getNumberOfCells())+ " cells")
    print("Initial data : ", "Riemann problem")
    print("Boundary conditions : ", "Neumann")
    print("Mesh name : ",filename )
    # Problem data
    tmax = 10.
    ntmax = 10
    cfl=1
    output_freq = 1
    EulerSystemStaggered(ntmax, tmax, cfl, output_freq, my_mesh, filename)
    return

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        filename=sys.argv[1]
        my_mesh = cdmath.Mesh(filename)
    else :
        filename = "CartesianGrid"
        ax=0;bx=1;nx=20;
        ay=0;by=1;ny=10;        
        my_mesh = cdmath.Mesh(ax,bx,nx,ay,by,ny)

    solve(my_mesh,filename,100)
