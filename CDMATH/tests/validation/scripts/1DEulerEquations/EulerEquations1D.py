#!/usr/bin/env python3
# -*-coding:utf-8 -*

"""
Created on Mon Aug 30 2021
@author: Michael NDJINGA, Katia Ait Ameur, Coraline Mounier

Euler system without heating source term (phi) in one dimension on regular domain [a,b]
Ghost cell (Neumann) boundary condition
Roe scheme or conservative MAC scheme
Regular square mesh

State law Stiffened gaz : p = (gamma - 1) * rho * (e - q) - gamma * pinf
4 choices of parameters gamma and pinf are available : 
  - Lagrange interpolation (with q=0)
  - Hermite interpolation with reference point at 575K (with q=0)
  - Hermite interpolation with reference point at 590K (with q=0)
  - Hermite interpolation with reference point at 617.94K (saturation at 155 bar)  with q=0
  
Linearized enthalpy : h = h_sat + cp * (T - T_sat)
Values for cp and T_sat parameters are taken at the reference point chosen for the state law

To do correct the computation of the time step : lambda_max (maximum eigenvalue) should be computed first)
"""

import cdmath
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, time, json
from math import sqrt
from numpy import sign


precision = 1.e-5

def p_to_e_StiffenedGaz(p, rho, gamma, pinf):
	e_field = (p + gamma*pinf) / (gamma - 1.) / rho
	return e_field

def initial_conditions_Riemann_problem(a, b, nx,p_L,v_L,rho_L,p_R,v_R,rho_R,pos_disc, gamma, pinf):
	print("Initial data Riemann problem")
	dx = (b - a) / nx  # space step
	x = [a + 0.5 * dx + i * dx for i in range(nx)]  # array of cell center (1D mesh)	

	p_initial   = np.array([ (xi < pos_disc) * p_L   + (xi >= pos_disc) * p_R   for xi in x])
	v_initial   = np.array([ (xi < pos_disc) * v_L   + (xi >= pos_disc) * v_R   for xi in x])
	rho_initial = np.array([ (xi < pos_disc) * rho_L + (xi >= pos_disc) * rho_R for xi in x])

	e_initial = p_to_e_StiffenedGaz(p_initial, rho_initial, gamma, pinf)
	q_initial = rho_initial * v_initial
	rho_E_initial = rho_initial*e_initial + 0.5 * q_initial**2/rho_initial

	return rho_initial, q_initial, rho_E_initial, p_initial, v_initial, e_initial

def p_to_rho_StiffenedGaz(p_field, T_field, gamma, pinf):
	rho_field = (p_field + pinf) * gamma / (gamma - 1) * 1. / (h_sat + cp * (T_field - T_sat))
	return rho_field
	
def T_to_rhoE_StiffenedGaz(T_field, rho_field, q_field, gamma, pinf):
	rho_E_field = 1. / 2. * (q_field) ** 2 / rho_field + pinf + rho_field / gamma * (h_sat + cp * (T_field- T_sat))
	return rho_E_field

def rho_to_p_StiffenedGaz(rho_field, q_field, rho_E_field, gamma, pinf):
	p_field = (gamma - 1) * (rho_E_field - 1. / 2 * q_field ** 2 / rho_field) - gamma * pinf
	return p_field

def T_to_E_StiffenedGaz(p_field, T_field, v_field, gamma, pinf):
	rho_field = p_to_rho_StiffenedGaz(p_field, T_field)
	E_field = (p_field + gamma * pinf) / ((gamma-1) * rho_field) + 0.5 * v_field **2
	return E_field
		
def dp_drho_e_const_StiffenedGaz( e , gamma, pinf):
	return (gamma-1)*(e)

def dp_de_rho_const_StiffenedGaz( rho, gamma, pinf ):
	return (gamma-1)*rho

def sound_speed_StiffenedGaz( h, gamma, pinf ):
	return np.sqrt((gamma-1)* h)

def rho_h_to_p_StiffenedGaz( rho, h, gamma, pinf ):
	return (gamma - 1) * rho * h / gamma - pinf

def rhoE_to_e_StiffenedGaz(rho_field, q_field, rho_E_field):
	e_field = rho_E_field / rho_field - 1. / 2. * (q_field / rho_field)**2
	return e_field

def MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf):
	RoeMat = cdmath.Matrix(3, 3)

	u_l = q_l / rho_l
	u_r = q_r / rho_r
	p_l = rho_to_p_StiffenedGaz(rho_l, q_l, rho_E_l, gamma, pinf)
	p_r = rho_to_p_StiffenedGaz(rho_r, q_r, rho_E_r, gamma, pinf)
	H_l = rho_E_l / rho_l + p_l / rho_l
	H_r = rho_E_r / rho_r + p_r / rho_r

	# Roe averages
	rho = np.sqrt(rho_l * rho_r)
	u   = (u_l * sqrt(rho_l) + u_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))
	H   = (H_l * sqrt(rho_l) + H_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))

	p = rho_h_to_p_StiffenedGaz( rho, H - u**2/2., gamma, pinf )
	e = H - p / rho - 1./2 * u**2
	dp_drho = dp_drho_e_const_StiffenedGaz( e, gamma, pinf )
	dp_de   = dp_de_rho_const_StiffenedGaz( rho, gamma, pinf )

	RoeMat[0, 0] = 0
	RoeMat[0, 1] = 1
	RoeMat[0, 2] = 0
	RoeMat[1, 0] = dp_drho - u ** 2 + dp_de / rho * (u**2/2 - e)
	RoeMat[1, 1] = 2 * u - u * dp_de / rho
	RoeMat[1, 2] = dp_de / rho
	RoeMat[2, 0] = -u * ( -dp_drho + H - dp_de / rho * (u**2/2 - e) )
	RoeMat[2, 1] = H - dp_de / rho * u ** 2
	RoeMat[2, 2] = (dp_de / rho +1) * u
	
	return(RoeMat)

	
def Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf):
	Dmac   = cdmath.Matrix(3, 3)

	u_l = q_l / rho_l
	u_r = q_r / rho_r
	p_l = rho_to_p_StiffenedGaz(rho_l, q_l, rho_E_l, gamma, pinf)
	p_r = rho_to_p_StiffenedGaz(rho_r, q_r, rho_E_r, gamma, pinf)
	H_l = rho_E_l / rho_l + p_l / rho_l
	H_r = rho_E_r / rho_r + p_r / rho_r

	# Roe averages
	rho = np.sqrt(rho_l * rho_r)
	u = (u_l * sqrt(rho_l) + u_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))
	H = (H_l * sqrt(rho_l) + H_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))

	p = rho_h_to_p_StiffenedGaz( rho, H - u**2/2., gamma, pinf )
	e = H - p / rho - 1./2 * u**2
	dp_drho = dp_drho_e_const_StiffenedGaz( e, gamma, pinf )
	dp_de   = dp_de_rho_const_StiffenedGaz( rho, gamma, pinf )
 
	#Third choice for Dstag
	Dmac[0, 0] = 0
	Dmac[0, 1] = 1
	Dmac[0, 2] = 0
	Dmac[1, 0] = -dp_drho - u ** 2 - dp_de / rho * (u**2/2 - e)
	Dmac[1, 1] = 2 * u + u * dp_de / rho
	Dmac[1, 2] = -dp_de / rho
	Dmac[2, 0] = -u * ( dp_drho + H + dp_de / rho * (u**2/2 - e) )
	Dmac[2, 1] = H + dp_de / rho * u ** 2
	Dmac[2, 2] = (-dp_de / rho +1) * u
	
	#return Dmac * sign(u)
	
	#Fifth choice for Dstag
	Dmac[0, 0] = abs(u)-u
	Dmac[0, 1] = 1
	Dmac[0, 2] = 0
	Dmac[1, 0] = -dp_drho - u ** 2 - dp_de / rho * (u**2/2 - e)
	Dmac[1, 1] = abs(u) + u + u * dp_de / rho
	Dmac[1, 2] = -dp_de / rho
	Dmac[2, 0] = -u * ( dp_drho + H + u*(u-abs(u)) + dp_de / rho * (u**2/2 - e) )
	Dmac[2, 1] = H + dp_de / rho * u ** 2
	Dmac[2, 2] = -u*dp_de / rho + abs(u)
	
	return Dmac 

def Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf):
    Droe   = cdmath.Matrix(3, 3)

    u_l = q_l / rho_l
    u_r = q_r / rho_r
    p_l = rho_to_p_StiffenedGaz(rho_l, q_l, rho_E_l, gamma, pinf)
    p_r = rho_to_p_StiffenedGaz(rho_r, q_r, rho_E_r, gamma, pinf)
    H_l = rho_E_l / rho_l + p_l / rho_l
    H_r = rho_E_r / rho_r + p_r / rho_r

    # Roe averages
    rho = np.sqrt(rho_l * rho_r)
    u = (u_l * sqrt(rho_l) + u_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))
    H = (H_l * sqrt(rho_l) + H_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))

    c = sound_speed_StiffenedGaz( H - u**2/2., gamma, pinf )
    
    lamb  = cdmath.Vector(3)
    lamb[0] = u-c
    lamb[1] = u
    lamb[2] = u+c    

    r   = cdmath.Matrix(3, 3)
    r[0,0] = 1.
    r[1,0] = u-c
    r[2,0] = H-u*c    
    r[0,1] = 1.
    r[1,1] = u   
    r[2,1] = H-c**2/(gamma-1)    
    r[0,2] = 1.
    r[1,2] = u+c
    r[2,2] = H+u*c         

    l   = cdmath.Matrix(3, 3)
    l[0,0] = (1./(2*c**2))*(0.5*(gamma-1)*u**2+u*c)
    l[1,0] = (1./(2*c**2))*(-u*(gamma-1)-c)
    l[2,0] = (1./(2*c**2))*(gamma-1)
    l[0,1] = ((gamma-1)/c**2)*(H-u**2)
    l[1,1] = ((gamma-1)/c**2)*u   
    l[2,1] = -((gamma-1)/c**2)    
    l[0,2] = (1./(2*c**2))*(0.5*(gamma-1)*u**2-u*c)
    l[1,2] = (1./(2*c**2))*(c-u*(gamma-1))
    l[2,2] = (1./(2*c**2))*(gamma-1)

    M1 = cdmath.Matrix(3, 3) #abs(lamb[0])*r[:,0].tensProduct(l[:,0])
    M2 = cdmath.Matrix(3, 3) #abs(lamb[1])*r[:,1].tensProduct(l[:,1])   
    M3 = cdmath.Matrix(3, 3) #abs(lamb[2])*r[:,2].tensProduct(l[:,2])
    for i in range(3):
        for j in range(3):
            M1[i,j] = abs(lamb[0])*r[i,0]*l[j,0]
            M2[i,j] = abs(lamb[1])*r[i,1]*l[j,1]            
            M3[i,j] = abs(lamb[2])*r[i,2]*l[j,2]
            
    Droe = M1+M2+M3 
    
    return(Droe)    


def jacobianMatricesm(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf):

	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)
	
	if scheme == 'Stag':
		Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)
		return (RoeMat - Dmac) * coeff * 0.5
	
	elif scheme == 'Roe':
		Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)    
		return (RoeMat - Droe) * coeff * 0.5
	else:
		raise ValueError("Wrong scheme given : ", scheme, ". Accepted scheme are Roe and Stag")


def jacobianMatricesp(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf):
	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)

	if scheme == 'Stag':
		Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)
		return (RoeMat - Dmac) * coeff * 0.5
	
	elif scheme == 'Roe':
		Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, gamma, pinf)    
		return (RoeMat - Droe) * coeff * 0.5

	else:
		raise ValueError("Wrong scheme given : ", scheme, ". Accepted scheme are Roe and Stag")

def FillEdges(nx, j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit, scheme, gamma, pinf):
	dUi1 = cdmath.Vector(3)
	dUi2 = cdmath.Vector(3)
	temp1 = cdmath.Vector(3)
	temp2 = cdmath.Vector(3)

	if (j == 0):
		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j + 1) * nbComp + 0]
		q_r     = Uk[(j + 1) * nbComp + 1]
		rho_E_r = Uk[(j + 1) * nbComp + 2]	

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf)
		divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
		divMat.addValue(j * nbComp,       j * nbComp, Am * (-1.))

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_l, q_l, rho_E_l, scheme, gamma, pinf)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
	
		if(isImplicit):
			dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
			dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
			dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
			temp1 = Am * dUi1
		else:
			dUi2[0] = rho_l   
			dUi2[1] = q_l     
			dUi2[2] = rho_E_l 
			temp2 = Ap * dUi2

	elif (j == nx - 1):
		rho_l   = Uk[(j - 1) * nbComp + 0]
		q_l     = Uk[(j - 1) * nbComp + 1]
		rho_E_l = Uk[(j - 1) * nbComp + 2]
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
		divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))

		Am = jacobianMatricesm(dt / dx, rho_r, q_r, rho_E_r, rho_r, q_r, rho_E_r, scheme, gamma, pinf)	
		divMat.addValue(j * nbComp, j * nbComp, Am * (-1.))

		if(isImplicit):
			dUi2[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
			dUi2[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
			dUi2[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]
			temp2 = Ap * dUi2
		else:
			dUi1[0] = rho_r   
			dUi1[1] = q_r     
			dUi1[2] = rho_E_r 
			temp1 = Am * dUi1
	
	if(isImplicit):#implicit scheme, contribution from the Newton scheme residual
		Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
		Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
		Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])
	else:#explicit scheme, contribution from the boundary data the right hand side
		Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] 
		Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] 
		Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] 

def FillInnerCell(j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit, scheme, gamma, pinf):

	rho_l   = Uk[(j - 1) * nbComp + 0]
	q_l     = Uk[(j - 1) * nbComp + 1]
	rho_E_l = Uk[(j - 1) * nbComp + 2]
	rho_r   = Uk[j * nbComp + 0]
	q_r     = Uk[j * nbComp + 1]
	rho_E_r = Uk[j * nbComp + 2]
	Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf)
	
	rho_l   = Uk[j * nbComp + 0]
	q_l     = Uk[j * nbComp + 1]
	rho_E_l = Uk[j * nbComp + 2]
	rho_r   = Uk[(j + 1) * nbComp + 0]
	q_r     = Uk[(j + 1) * nbComp + 1]
	rho_E_r = Uk[(j + 1) * nbComp + 2]
	Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme, gamma, pinf)

	divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
	divMat.addValue(j * nbComp, j * nbComp, Am * (-1.))
	divMat.addValue(j * nbComp, j * nbComp, Ap)
	divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))

	if(isImplicit):
		dUi1 = cdmath.Vector(3)
		dUi2 = cdmath.Vector(3)
		
		dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
	
		dUi2[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi2[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi2[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]
		
		temp1 = Am * dUi1
		temp2 = Ap * dUi2

		Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
		Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
		Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])
	else:
		Rhs[j * nbComp + 0] = 0
		Rhs[j * nbComp + 1] = 0
		Rhs[j * nbComp + 2] = 0


def EulerSystem(ntmax, tmax, cfl, a, b, nbCells, output_freq, meshName, isImplicit, scheme,fig,p_L,v_L,rho_L,p_R,v_R,rho_R, lam_max, gamma, pinf,pos_disc, test_desc):
	dim = 1
	nbComp = 3
	dt = 0.
	time = 0.
	it = 0
	isStationary = False
	dx = (b - a) / nbCells
	dt = cfl * dx / lam_max
	nbVoisinsMax = 2
	
	e_L=p_to_e_StiffenedGaz(p_L, rho_L, gamma, pinf)
	
	# iteration vectors
	Un  = cdmath.Vector(nbCells * (nbComp))
	dUn = cdmath.Vector(nbCells * (nbComp))
	dUk = cdmath.Vector(nbCells * (nbComp))
	Rhs = cdmath.Vector(nbCells * (nbComp))

	# Initial conditions
	print("Construction of the initial condition …")

	rho_field, q_field, rho_E_field, p_field, v_field, e_field = initial_conditions_Riemann_problem(a, b, nbCells, p_L,v_L,rho_L,p_R,v_R,rho_R,pos_disc, gamma, pinf)
	h_field = (rho_E_field + p_field) / rho_field - 0.5 * (q_field / rho_field) **2
	p_field = p_field * 10 ** (-5)
	

	for k in range(nbCells):
		Un[k * nbComp + 0] = rho_field[k]
		Un[k * nbComp + 1] = q_field[k]
		Un[k * nbComp + 2] = rho_E_field[k]

	divMat = cdmath.SparseMatrixPetsc(nbCells * nbComp, nbCells * nbComp, (nbVoisinsMax + 1) * nbComp)
	
	iterGMRESMax = 50
	newton_max = 100

	print("Starting computation of the non linear Euler non isentropic system with scheme : ", scheme)
	# STARTING TIME LOOP
	while (it < ntmax and time <= tmax and not isStationary):
		dUn = Un.deepCopy()
		Uk  = Un.deepCopy()
		residu = 1.
		
		k = 0
		while (k < newton_max and residu > precision):
			# STARTING NEWTON LOOP
			divMat.zeroEntries()  #sets the matrix coefficients to zero
			for j in range(nbCells):
				
				# traitements des bords
				if (j == 0):
					FillEdges(nbCells,j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit, scheme, gamma, pinf)
				elif (j == nbCells - 1):
					FillEdges(nbCells,j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit, scheme, gamma, pinf)

				# traitement des cellules internes
				else:
					FillInnerCell(    j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit, scheme, gamma, pinf)
					
			if(isImplicit):
				#solving the linear system divMat * dUk = Rhs
				divMat.diagonalShift(1.)
				LS = cdmath.LinearSolver(divMat, Rhs, iterGMRESMax, precision, "GMRES", "LU")
				dUk = LS.solve()
				vector_residu = dUk.maxVector(nbComp)
				residu = max(abs(vector_residu[0])/rho_L, abs(vector_residu[1])/(rho_L*v_L), abs(vector_residu[2])/rho_L*e_L )
			else:
				dUk=Rhs - divMat*Un
				residu = 0.#Convergence schéma Newton
			
			if (isImplicit and (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax)):
				print("Residu Newton at iteration ",k, " :   ", residu)
				print("Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations")

			#updates for Newton loop
			Uk += dUk
			k = k + 1
			if (isImplicit and not LS.getStatus()):
				print("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
				raise ValueError("No convergence of the linear system")
			
			if k == newton_max:
				raise ValueError("No convergence of Newton algorithm for " + scheme + " scheme")

		#updating fields
		Un = Uk.deepCopy()
		dUn -= Un

		#Testing stationarity
		residu_stat = dUn.maxVector(nbComp)#On prend le max de chaque composante
		if (it % output_freq == 0 ):
			print("Test de stationarité : Un+1-Un= ", max(abs(residu_stat[0])/rho_L, abs(residu_stat[1])/(rho_L*v_L), abs(residu_stat[2])/rho_L*e_L ))

		if ( it>1 and abs(residu_stat[0])/rho_L<precision  and abs(residu_stat[1])/(rho_L*v_L)<precision and abs(residu_stat[2])/rho_L*e_L<precision):
				isStationary = True
		
		for k in range(nbCells):
			rho_field[k]   = Un[k * nbComp + 0]
			q_field[k]     = Un[k * nbComp + 1]
			rho_E_field[k] = Un[k * nbComp + 2]

		v_field = q_field / rho_field
		p_field = rho_to_p_StiffenedGaz(rho_field, q_field, rho_E_field, gamma, pinf)
		e_field = rhoE_to_e_StiffenedGaz(rho_field, q_field, rho_E_field )
		h_field = (rho_E_field + p_field) / rho_field - 0.5 * (q_field / rho_field) **2
		p_field = p_field * 10 ** (-5)
		
		if( min(p_field)<0) :
			raise ValueError("Negative pressure, stopping calculation")

		time = time + dt
		it = it + 1

		# Savings
		if (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax):
			print("-- Time step : " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

	print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)+"\n")
	
	if(isImplicit):
		test_desc["Linear_solver_algorithm"]=LS.getNameOfMethod()
		test_desc["Linear_solver_preconditioner"]=LS.getNameOfPc()
		test_desc["Linear_solver_precision"]=LS.getTolerance()
		test_desc["Linear_solver_maximum_iterations"]=LS.getNumberMaxOfIter()
		test_desc["Linear_system_max_actual_iterations_number"]=LS.getNumberOfIter()
		test_desc["Linear_system_max_actual_error"]=LS.getResidu()

	if (it >= ntmax):
		print("Maximum number of time steps ntmax= ", ntmax, " reached")

	elif (isStationary):
		print("Stationary regime reached at time step ", it, ", t= ", time)
		print("------------------------------------------------------------------------------------")
		np.savetxt("EulerComplet_" + str(dim) + "_" + scheme + meshName + "_rho_Stat.txt", rho_field, delimiter="\n")
		np.savetxt("EulerComplet_" + str(dim) + "_" + scheme + meshName + "_q_Stat.txt", q_field, delimiter="\n")
		np.savetxt("EulerComplet_" + str(dim) + "_" + scheme + meshName + "_rhoE_Stat.txt", rho_E_field, delimiter="\n")
		np.savetxt("EulerComplet_" + str(dim) + "_" + scheme + meshName + "_p_Stat.txt", p_field, delimiter="\n")
	else:
		print("Maximum time Tmax= ", tmax, " reached")

	return 	p_field, v_field, e_field, rho_field, q_field, h_field


def solve(a, b, nx, meshName, meshType, isImplicit, scheme, simulation_name, testColor,fig,p_L,v_L,rho_L,p_R,v_R,rho_R, pos_disc, tmax, lam_max, gamma, pinf):

	print("Simulation of Euler equations in dimension 1 on " + str(nx) + " cells")
	print("Problem name ", simulation_name)
	print("State Law Stiffened Gaz, gamma= " , gamma, " pinf= ", pinf)
	print("Initial data : ", "Riemann data")
	print("Boundary conditions : ", "Neumann")
	print("Mesh name : ", meshName, ", ", nx, " cells")

	# Problem data
	ntmax = 100000
	output_freq = 1000
	cfl=0.95
	
	test_desc={}
	test_desc["Initial_data"]=simulation_name
	test_desc["Boundary_conditions"]="Neumann"
	test_desc["Global_name"]="FV simulation of "+simulation_name
	test_desc["Global_comment"]="Regular 1D mesh"
	test_desc["PDE_model"]="Euler equations"
	test_desc["PDE_is_stationary"]=False
	test_desc["PDE_search_for_stationary_solution"]=False
	test_desc["Numerical_method_name"]=scheme
	test_desc["Numerical_method_space_discretization"]="Finite volumes"
	test_desc["Numerical_method_time_discretization"]="Explicit"
	test_desc["Mesh_is_unstructured"]=False
	test_desc["Mesh_type"]=meshType
	test_desc["Test_color"]=testColor
	test_desc["Geometry"]="Segment"
	test_desc["Part_of_mesh_convergence_analysis"]=True
	test_desc["Space_dimension"]=1
	test_desc["Mesh_dimension"]=1
	test_desc["Mesh_number_of_elements"]=nx
	test_desc["Mesh_cell_type"]='SEG'
	test_desc["Mesh_max_number_of_neighbours"]=2

	start = time.time()
	p_field, v_field, e_field, rho_field, q_field, h_field = EulerSystem(ntmax, tmax, cfl, a, b, nx, output_freq, meshName,  isImplicit,scheme,fig,p_L,v_L,rho_L,p_R,v_R,rho_R, pos_disc,lam_max, gamma, pinf, test_desc)
	end = time.time()

	test_desc["Computational_time_taken_by_run"]=end-start
	with open('test_EulerEquations_1D_VF_'+simulation_name+"_"+str(nx)+ "Cells.json", 'w') as outfile:  
		json.dump(test_desc, outfile)

	
	return p_field, v_field, e_field, rho_field, q_field, h_field, end - start
	

if __name__ == """__main__""":
	nbComp=3 # number of equations 
	a = 0.# domain is interval [a,b]
	b = 4.2# domain is interval [a,b]
	nx = 10# number of cells
	dx = (b - a) / nx  # space step
	x = [a + 0.5 * dx + i * dx for i in range(nx)]  # array of cell center (1D mesh)
	rho0=p_to_rho_StiffenedGaz(p_0, T_0)
	rhoE0=T_to_rhoE_StiffenedGaz(T_0, rho0, rho0*v_0)


#### initial condition (T in K, v in m/s, p in Pa)
	p_initial   = np.array([ p_outlet      for xi in x])
	v_initial   = np.array([ v_inlet       for xi in x])
	T_initial   = np.array([ T_inlet       for xi in x])
	
	rho_field = p_to_rho_StiffenedGaz(p_initial, T_initial)
	q_field = rho_field * v_initial
	rho_E_field = rho_field * T_to_E_StiffenedGaz(p_initial, T_initial, v_initial)

	U = cdmath.Vector(nx * (nbComp))#Inutile à terme mais nécessaire pour le moment

	for k in range(nx):
		U[k * nbComp + 0] = rho_field[k]
		U[k * nbComp + 1] = q_field[k]
		U[k * nbComp + 2] = rho_E_field[k]

	print("\n #############   Testing Roe scheme #############\n")
	scheme='Roe'

	print("\n Testing function solve (Implicit Roe scheme) \n")
	isImplicit=True
	cfl = 1000.
	solve(a, b, nx, "RegularSquares", "", cfl,  isImplicit, scheme, gamma, pinf)  #Implicit Roe simulation
	
	print("\n Testing function solve (Explicit Roe scheme) \n")
	isImplicit=False
	cfl = 0.5
	solve(a, b, nx, "RegularSquares", "", cfl, isImplicit, scheme, gamma, pinf)  #Explicit Roe simulation

	print("\n #############   Testing MAC scheme #############\n")
	scheme='Mac'
	isImplicit=True #Scheme must be implicit (explicit version is unstable)

	print("\n Testing function solve (Implicit MAC scheme) \n")
	isImplicit=True
	cfl = 1000.
	solve(a, b, nx, "RegularSquares", "", cfl, isImplicit, scheme, gamma, pinf)  #Implicit MAC simulation
	
