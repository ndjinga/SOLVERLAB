#!/usr/bin/env python3
# -*-coding:utf-8 -*

"""
Created on Mon Aug 30 2021
@author: Michael NDJINGA, Katia Ait Ameur, Coraline Mounier

Euler system with heating source term (phi) in one dimension on regular domain [a,b]
Riemann problemn with ghost cell boundary condition
Left : Inlet boundary condition (velocity and temperature imposed)
Right : Outlet boundary condition (pressure imposed)
Roe scheme
Regular square mesh

State law Stiffened gaz : p = (gamma - 1) * rho * (e - q) - gamma * p0
4 choices of parameters gamma and p0 are available : 
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
import matplotlib.animation as manimation
import sys
from math import sqrt, atan, pi
from numpy import sign


#### Initial and boundary condition (T in K, v in m/s, p in Pa)
T_inlet  = 565.
v_inlet  = 5.
p_outlet = 155.0 * 10**5

#initial parameters are determined from boundary conditions
p_0   = p_outlet       #initial pressure
v_0   = v_inlet        #initial velocity
T_0   = T_inlet        #initial temperature
### Heating source term
phi=1.e8

## Numerical parameter
precision = 1e-6

#state law parameter : can be 'Lagrange', 'Hermite590K', 'Hermite617K', or 'FLICA'
state_law = "Hermite575K"

def state_law_parameters(state_law):
	#state law Stiffened Gaz : p = (gamma - 1) * rho * e - gamma * p0
	global gamma
	global p0
	global q
	global c0
	global cp
	global h_sat
	global T_sat
	
	if state_law == "Lagrange":
		# reference values for Lagrange interpolation
		p_ref = 155. * 10**5     #Reference pressure in a REP 900 nuclear power plant     
		p1    = 153. * 10**5     # value of pressure at inlet of a 900 MWe PWR vessel
		rho_ref = 594.38        #density of water at saturation temperature of 617.94K and 155 bars
		rho1 = 742.36           # value of density at inlet of a 900 MWe PWR vessel (T1 = 565K)
		e_ref = 1603.8 * 10**3  #internal energy of water at saturation temperature of 617.94K and 155 bars
		e1 = 1273.6 * 10**3     # value of internal energy at inlet of a 900 MWe PWR vessel
		
		gamma = (p1 - p_ref) / (rho1 * e1 - rho_ref *e_ref) + 1.
		p0 = - 1. / gamma * ( - (gamma - 1) * rho_ref * e_ref + p_ref)
		q=0.
		c0 = sqrt((gamma - 1) * (e_ref + p_ref / rho_ref))
		
		cp = 8950.
		h_sat = 1.63 * 10 ** 6     # saturation enthalpy of water at 155 bars
		T_sat = 617.94 

	elif state_law == "Hermite617K":
		# reference values for Hermite interpolation
		p_ref = 155. * 10**5     #Reference pressure in a REP 900 nuclear power plant
		T_ref = 617.94          #Reference temperature for interpolation at 617.94K
		rho_ref = 594.38        #density of water at saturation temperature of 617.94K and 155 bars
		e_ref = 1603.8 * 10**3  #internal energy of water at saturation temperature of 617.94K and 155 bars
		h_ref   = e_ref + p_ref / rho_ref
		c_ref = 621.43          #sound speed for water at 155 bars and 617.94K

		gamma = 1. + c_ref * c_ref / (e_ref + p_ref / rho_ref)           # From the sound speed formula
		p0 = 1. / gamma * ( (gamma - 1) * rho_ref * e_ref - p_ref)       
		q=0.
		c0 = sqrt((gamma - 1) * (e_ref + p_ref / rho_ref))
		
		cp = 8950.                  # value at 155 bar and 617.94K
		h_sat = 1.63 * 10 ** 6     # saturation enthalpy of water at 155 bars
		T_sat = 617.94  
	
	elif state_law == 'Hermite590K':
		# reference values for Hermite interpolation
		p_ref = 155. * 10**5     #Reference pressure  in a REP 900 nuclear power plant
		T_ref = 590.             #Reference temperature for interpolation at 590K
		rho_ref = 688.3         #density of water at 590K and 155 bars
		e_ref = 1411.4 * 10**3  #internal energy of water at 590K and 155 bars
		h_ref   = e_ref + p_ref / rho_ref
		c_ref = 866.29          #sound speed for water at 155 bars and 590K
		
		gamma = 1. + c_ref * c_ref / (e_ref + p_ref / rho_ref)           # From the sound speed formula
		p0 = 1. / gamma * ( (gamma - 1) * rho_ref * e_ref - p_ref)       
		q=0.
		c0 = sqrt((gamma - 1) * (e_ref + p_ref / rho_ref))
		
		cp = 5996.8                  # value at 155 bar and 590K
		h_sat = 1433.9 * 10 ** 3     # saturation enthalpy of water at 155 bars
		T_sat = 590.  

	elif state_law == 'Hermite575K':
		# reference values for Hermite interpolation
		p_ref = 155 * 10**5     #Reference pressure  in a REP 900 nuclear power plant
		T_ref = 575             #Reference temperature at inlet in a REP 900 nuclear power plant
		#Remaining values determined using iapws python package
		rho_ref = 722.66        #density of water at 575K and 155 bars
		e_ref = 1326552.66  #internal energy of water at 575K and 155 bars
		h_ref   = e_ref + p_ref / rho_ref
		c_ref = 959.28          #sound speed for water at 155 bars and 575K
		
		gamma = 1 + c_ref * c_ref / (e_ref + p_ref / rho_ref)           # From the sound speed formula
		p0 = 1 / gamma * ( (gamma - 1) * rho_ref * e_ref - p_ref)       
		q=0.
		c0 = sqrt((gamma - 1) * (e_ref + p_ref / rho_ref))#This is actually c_ref
		
		cp = 5504.05                 # value at 155 bar and 590K
		h_sat = h_ref     # saturation enthalpy of water at 155 bars
		T_sat = T_ref
	else:
		raise ValueError("Incorrect value for parameter state_law")
		
def initial_conditions_Riemann_problem(a, b, nx):
	print("Initial data Riemann problem")
	dx = (b - a) / nx  # space step
	x = [a + 0.5 * dx + i * dx for i in range(nx)]  # array of cell center (1D mesh)

	p_initial = np.array([ p_0 for xi in x])
	v_initial = np.array([ v_0 for xi in x])
	T_initial = np.array([ T_0 for xi in x])

	rho_initial = p_to_rho_StiffenedGaz(p_initial, T_initial)
	q_initial = rho_initial * v_initial
	rho_E_initial = T_to_rhoE_StiffenedGaz(T_initial, rho_initial, q_initial)

	return rho_initial, q_initial, rho_E_initial, p_initial, v_initial, T_initial

def p_to_rho_StiffenedGaz(p_field, T_field):
	rho_field = (p_field + p0) * gamma / (gamma - 1) * 1. / (h_sat + cp * (T_field - T_sat))
	return rho_field
	
def T_to_rhoE_StiffenedGaz(T_field, rho_field, q_field):
	rho_E_field = 1. / 2. * (q_field) ** 2 / rho_field + p0 + rho_field / gamma * (h_sat + cp * (T_field- T_sat))
	return rho_E_field

def rhoE_to_T_StiffenedGaz(rho_field, q_field, rho_E_field):
	T_field = T_sat + 1 / cp * (gamma * (rho_E_field / rho_field - 1 / 2 * (q_field / rho_field) ** 2) - gamma * p0 / rho_field - h_sat)
	return T_field

def rho_to_p_StiffenedGaz(rho_field, q_field, rho_E_field):
	p_field = (gamma - 1) * (rho_E_field - 1. / 2 * q_field ** 2 / rho_field) - gamma * p0
	return p_field

def T_to_E_StiffenedGaz(p_field, T_field, v_field):
	rho_field = p_to_rho_StiffenedGaz(p_field, T_field)
	E_field = (p_field + gamma * p0) / ((gamma-1) * rho_field) + 0.5 * v_field **2
	return E_field
		
def dp_drho_e_const_StiffenedGaz( e ):
	return (gamma-1)*(e-q)

def dp_de_rho_const_StiffenedGaz( rho ):
	return (gamma-1)*rho

def sound_speed_StiffenedGaz( h ):
	return np.sqrt((gamma-1)*(h-q))

def rho_h_to_p_StiffenedGaz( rho, h ):
	return (gamma - 1) * rho * ( h - q ) / gamma - p0

def MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):
	RoeMat = cdmath.Matrix(3, 3)

	u_l = q_l / rho_l
	u_r = q_r / rho_r
	p_l = rho_to_p_StiffenedGaz(rho_l, q_l, rho_E_l)
	p_r = rho_to_p_StiffenedGaz(rho_r, q_r, rho_E_r)
	H_l = rho_E_l / rho_l + p_l / rho_l
	H_r = rho_E_r / rho_r + p_r / rho_r

	# Roe averages
	rho = np.sqrt(rho_l * rho_r)
	u   = (u_l * sqrt(rho_l) + u_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))
	H   = (H_l * sqrt(rho_l) + H_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))

	p = rho_h_to_p_StiffenedGaz( rho, H - u**2/2. )
	e = H - p / rho - 1./2 * u**2
	dp_drho = dp_drho_e_const_StiffenedGaz( e )
	dp_de   = dp_de_rho_const_StiffenedGaz( rho )

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

	
def Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):
    Droe   = cdmath.Matrix(3, 3)

    u_l = q_l / rho_l
    u_r = q_r / rho_r
    p_l = rho_to_p_StiffenedGaz(rho_l, q_l, rho_E_l)
    p_r = rho_to_p_StiffenedGaz(rho_r, q_r, rho_E_r)
    H_l = rho_E_l / rho_l + p_l / rho_l
    H_r = rho_E_r / rho_r + p_r / rho_r

    # Roe averages
    rho = np.sqrt(rho_l * rho_r)
    u = (u_l * sqrt(rho_l) + u_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))
    H = (H_l * sqrt(rho_l) + H_r * sqrt(rho_r)) / (sqrt(rho_l) + sqrt(rho_r))

    c = sound_speed_StiffenedGaz( H - u**2/2. )
    
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


def jacobianMatricesm(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):

	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    
	return (RoeMat - Droe) * coeff * 0.5


def jacobianMatricesp(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):
	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    

	return (RoeMat + Droe) * coeff * 0.5


def FillEdges(j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit):
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

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
		divMat.addValue(j * nbComp,       j * nbComp, Am * (-1.))

		p_inlet = rho_to_p_StiffenedGaz(Uk[j * nbComp + 0], Uk[j * nbComp + 1], Uk[j * nbComp + 2])# We take p from inside the domain
		rho_l=p_to_rho_StiffenedGaz(p_inlet, T_inlet) # rho is computed from the temperature BC and the inner pressure
		#rho_l   = Uk[j * nbComp + 0]                            # We take rho from inside the domain
		q_l     = rho_l * v_inlet                               # q is imposed by the boundary condition v_inlet
		rho_E_l = T_to_rhoE_StiffenedGaz(T_inlet, rho_l, q_l)   #rhoE is obtained  using the two boundary conditions v_inlet and e_inlet
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
	
		if(isImplicit):
			dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
			dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
			dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
			temp1 = Am * dUi1
		
			dUi2[0] = rho_l   -  Uk[(j ) * nbComp + 0]
			dUi2[1] = q_l     -  Uk[(j ) * nbComp + 1]
			dUi2[2] = rho_E_l -  Uk[(j ) * nbComp + 2]
			temp2 = Ap * dUi2
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

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
		divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))

		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j ) * nbComp + 0]                               # We take rho inside the domain
		q_r     = Uk[(j ) * nbComp + 1]                               # We take q from inside the domain
		rho_E_r = (p_outlet+gamma*p0)/(gamma-1) + 0.5*q_r**2/rho_r    # rhoE is obtained using the boundary condition p_outlet

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)	
		divMat.addValue(j * nbComp, j * nbComp, Am * (-1.))

		if(isImplicit):
			dUi1[0] = rho_r   - Uk[j * nbComp + 0]
			dUi1[1] = q_r     - Uk[j * nbComp + 1]
			dUi1[2] = rho_E_r - Uk[j * nbComp + 2]
			temp1 = Am * dUi1

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

def FillInnerCell(j, Uk, nbComp, divMat, Rhs, Un, dt, dx, isImplicit):

	rho_l   = Uk[(j - 1) * nbComp + 0]
	q_l     = Uk[(j - 1) * nbComp + 1]
	rho_E_l = Uk[(j - 1) * nbComp + 2]
	rho_r   = Uk[j * nbComp + 0]
	q_r     = Uk[j * nbComp + 1]
	rho_E_r = Uk[j * nbComp + 2]
	Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	rho_l   = Uk[j * nbComp + 0]
	q_l     = Uk[j * nbComp + 1]
	rho_E_l = Uk[j * nbComp + 2]
	rho_r   = Uk[(j + 1) * nbComp + 0]
	q_r     = Uk[(j + 1) * nbComp + 1]
	rho_E_r = Uk[(j + 1) * nbComp + 2]
	Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)

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

def SetPicture(rho_field_Roe, q_field_Roe, h_field_Roe, p_field_Roe, v_field_Roe, T_field_Roe, dx):
	max_initial_rho = max(rho_field_Roe)
	min_initial_rho = min(rho_field_Roe)
	max_initial_q = max(q_field_Roe)
	min_initial_q = min(q_field_Roe)
	min_initial_h = min(h_field_Roe)
	max_initial_h = max(h_field_Roe)
	max_initial_p = max(p_field_Roe)
	min_initial_p = min(p_field_Roe)
	max_initial_v = max(v_field_Roe)
	min_initial_v = min(v_field_Roe)
	max_initial_T = max(T_field_Roe)
	min_initial_T = min(T_field_Roe)

	fig, ([axDensity, axMomentum, axh],[axPressure, axVitesse, axTemperature]) = plt.subplots(2, 3,sharex=True, figsize=(14,10))
	plt.gcf().subplots_adjust(wspace = 0.5)

	lineDensity_Roe, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field_Roe, label='Roe')
	axDensity.set(xlabel='x (m)', ylabel='Densité (kg/m3)')
	axDensity.set_xlim(a,b)
	axDensity.set_ylim(680, 800)
	axDensity.legend()

	lineMomentum_Roe, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field_Roe, label='Roe')
	axMomentum.set(xlabel='x (m)', ylabel='Momentum (kg/(m2.s))')
	axMomentum.set_xlim(a,b)
	axMomentum.set_ylim(3500, 	4000)
	axMomentum.legend()

	lineh_Roe, = axh.plot([a+0.5*dx + i*dx for i in range(nx)], h_field_Roe, label='Roe')
	axh.set(xlabel='x (m)', ylabel='h (J/m3)')
	axh.set_xlim(a,b)
	axh.set_ylim(1.2 * 10**6, 1.5*10**6)
	axh.legend()
	
	linePressure_Roe, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field_Roe, label='Roe')
	axPressure.set(xlabel='x (m)', ylabel='Pression (bar)')
	axPressure.set_xlim(a,b)
	axPressure.set_ylim(155, 155.5)
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axPressure.legend()

	lineVitesse_Roe, = axVitesse.plot([a+0.5*dx + i*dx for i in range(nx)], v_field_Roe, label='Roe')
	axVitesse.set(xlabel='x (m)', ylabel='Vitesse (m/s)')
	axVitesse.set_xlim(a,b)
	axVitesse.set_ylim(v_0-1, v_0+1)
	axVitesse.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axVitesse.legend()

	lineTemperature_Roe, = axTemperature.plot([a+0.5*dx + i*dx for i in range(nx)], T_field_Roe, label='Roe')
	axTemperature.set(xlabel='x (m)', ylabel='Température (K)')
	axTemperature.set_xlim(a,b)
	axTemperature.set_ylim(T_0-10, T_0+30)
	axTemperature.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axTemperature.legend()
	
	return(fig, lineDensity_Roe, lineMomentum_Roe, lineh_Roe, linePressure_Roe, lineVitesse_Roe, lineTemperature_Roe)


def EulerSystemRoe(ntmax, tmax, cfl, a, b, nbCells, output_freq, meshName, state_law, isImplicit):
	state_law_parameters(state_law)
	dim = 1
	nbComp = 3
	dt = 0.
	time = 0.
	it = 0
	isStationary = False
	dx = (b - a) / nx
	dt = cfl * dx / c0
	#dt = 5*10**(-6)
	nbVoisinsMax = 2

	# iteration vectors
	Un_Roe  = cdmath.Vector(nbCells * (nbComp))
	dUn_Roe = cdmath.Vector(nbCells * (nbComp))
	dUk_Roe = cdmath.Vector(nbCells * (nbComp))
	Rhs_Roe = cdmath.Vector(nbCells * (nbComp))

	# Initial conditions
	print("Construction of the initial condition …")

	rho_field_Roe, q_field_Roe, rho_E_field_Roe, p_field_Roe, v_field_Roe, T_field_Roe = initial_conditions_Riemann_problem(a, b, nx)
	h_field_Roe = (rho_E_field_Roe + p_field_Roe) / rho_field_Roe - 0.5 * (q_field_Roe / rho_field_Roe) **2
	p_field_Roe = p_field_Roe * 10 ** (-5)
	

	for k in range(nbCells):
		Un_Roe[k * nbComp + 0] = rho_field_Roe[k]
		Un_Roe[k * nbComp + 1] = q_field_Roe[k]
		Un_Roe[k * nbComp + 2] = rho_E_field_Roe[k]

	divMat_Roe = cdmath.SparseMatrixPetsc(nbCells * nbComp, nbCells * nbComp, (nbVoisinsMax + 1) * nbComp)
	
	# Picture settings
	fig, lineDensity_Roe, lineMomentum_Roe, lineRhoE_Roe, linePressure_Roe, lineVitesse_Roe, lineTemperature_Roe = SetPicture( rho_field_Roe, q_field_Roe, h_field_Roe, p_field_Roe, v_field_Roe, T_field_Roe, dx)

	plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "0" + ".png")
	iterGMRESMax = 50
	newton_max = 100

	print("Starting computation of the non linear Euler non isentropic system with Roe scheme …")
	# STARTING TIME LOOP
	while (it < ntmax and time <= tmax and not isStationary):
		dUn_Roe = Un_Roe.deepCopy()
		Uk_Roe  = Un_Roe.deepCopy()
		residu_Roe = 1.
		
		k_Roe = 0
		while (k_Roe < newton_max and residu_Roe > precision):
			# STARTING NEWTON LOOP
			divMat_Roe.zeroEntries()  #sets the matrix coefficients to zero
			for j in range(nbCells):
				
				# traitements des bords
				if (j == 0):
					FillEdges(j, Uk_Roe, nbComp, divMat_Roe, Rhs_Roe, Un_Roe, dt, dx, isImplicit)
				elif (j == nbCells - 1):
					FillEdges(j, Uk_Roe, nbComp, divMat_Roe, Rhs_Roe, Un_Roe, dt, dx, isImplicit)

				# traitement des cellules internes
				else:
					FillInnerCell(j, Uk_Roe, nbComp, divMat_Roe, Rhs_Roe, Un_Roe, dt, dx, isImplicit)
					
				Rhs_Roe[j * nbComp + 2] += phi*dt
			
			if(isImplicit):
				#solving the linear system divMat * dUk = Rhs
				divMat_Roe.diagonalShift(1.)
				LS_Roe = cdmath.LinearSolver(divMat_Roe, Rhs_Roe, iterGMRESMax, precision, "GMRES", "LU")
				dUk_Roe = LS_Roe.solve()
				vector_residu_Roe = dUk_Roe.maxVector(nbComp)
				residu_Roe = max(abs(vector_residu_Roe[0])/rho0, abs(vector_residu_Roe[1])/(rho0*v_0), abs(vector_residu_Roe[2])/rhoE0 )
			else:
				dUk_Roe=Rhs_Roe - divMat_Roe*Un_Roe
				residu_Roe = 0.#Convergence schéma Newton
			
			if (isImplicit and (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax)):
				print("Residu Newton at iteration ",k_Roe, " :   ", residu_Roe)
				print("Linear system converged in ", LS_Roe.getNumberOfIter(), " GMRES iterations")

			#updates for Newton loop
			Uk_Roe += dUk_Roe
			k_Roe = k_Roe + 1
			if (isImplicit and not LS_Roe.getStatus()):
				print("Linear system did not converge ", LS_Roe.getNumberOfIter(), " GMRES iterations")
				raise ValueError("No convergence of the linear system")
			
			if k_Roe == newton_max:
				raise ValueError("No convergence of Newton Roe Scheme")

		#updating fields
		Un_Roe = Uk_Roe.deepCopy()
		dUn_Roe -= Un_Roe

		#Testing stationarity
		residu_stat = dUn_Roe.maxVector(nbComp)#On prend le max de chaque composante
		if (it % output_freq == 0 ):
			print("Test de stationarité : Un+1-Un= ", max(abs(residu_stat[0])/rho0, abs(residu_stat[1])/(rho0*v_0), abs(residu_stat[2])/rhoE0 ))

		if ( it>1 and abs(residu_stat[0])/rho0<precision  and abs(residu_stat[1])/(rho0*v_0)<precision and abs(residu_stat[2])/rhoE0<precision):
				isStationary = True
		
		for k in range(nbCells):
			rho_field_Roe[k]   = Un_Roe[k * nbComp + 0]
			q_field_Roe[k]     = Un_Roe[k * nbComp + 1]
			rho_E_field_Roe[k] = Un_Roe[k * nbComp + 2]

		v_field_Roe = q_field_Roe / rho_field_Roe
		p_field_Roe = rho_to_p_StiffenedGaz(rho_field_Roe, q_field_Roe, rho_E_field_Roe)
		T_field_Roe = rhoE_to_T_StiffenedGaz(rho_field_Roe, q_field_Roe, rho_E_field_Roe)
		h_field_Roe = (rho_E_field_Roe + p_field_Roe) / rho_field_Roe - 0.5 * (q_field_Roe / rho_field_Roe) **2
		p_field_Roe = p_field_Roe * 10 ** (-5)
		
		if( min(p_field_Roe)<0) :
			raise ValueError("Negative pressure, stopping calculation")

		#picture and video updates
		lineDensity_Roe.set_ydata(rho_field_Roe)
		lineMomentum_Roe.set_ydata(q_field_Roe)
		lineRhoE_Roe.set_ydata(h_field_Roe)
		linePressure_Roe.set_ydata(p_field_Roe)
		lineVitesse_Roe.set_ydata(v_field_Roe)
		lineTemperature_Roe.set_ydata(T_field_Roe)
		
		time = time + dt
		it = it + 1

		# Savings
		if (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax):
	
			print("-- Time step : " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

			print("Temperature gain between inlet and outlet is ", T_field_Roe[nbCells-1]-T_field_Roe[0],"\n")

			plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "Implicit"+str(isImplicit)+ str(it) + '_time' + str(time) + ".png")

	print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)+"\n")
	
	if (it >= ntmax):
		print("Maximum number of time steps ntmax= ", ntmax, " reached")
		return

	elif (isStationary):
		print("Stationary regime reached at time step ", it, ", t= ", time)
		print("------------------------------------------------------------------------------------")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "_rho_Stat.txt", rho_field_Roe, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "_q_Stat.txt", q_field_Roe, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "_rhoE_Stat.txt", rho_E_field_Roe, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "_p_Stat.txt", p_field_Roe, delimiter="\n")
		plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_Roe" + meshName + "Implicit"+str(isImplicit)+"_Stat.png")
		return
	else:
		print("Maximum time Tmax= ", tmax, " reached")
		return


def solve(a, b, nx, meshName, meshType, cfl, state_law, isImplicit):
	print("Simulation of a heated channel in dimension 1 on " + str(nx) + " cells")
	print("State Law Stiffened Gaz, " + state_law)
	print("Initial data : ", "constant fields")
	print("Boundary conditions : ", "Inlet (Left), Outlet (Right)")
	print("Mesh name : ", meshName, ", ", nx, " cells")
	# Problem data
	tmax = 10.
	ntmax = 100000
	output_freq = 1000
	EulerSystemRoe(ntmax, tmax, cfl, a, b, nx, output_freq, meshName, state_law, isImplicit)
	return

def FillMatrixFromEdges(j, Uk, nbComp, divMat, dt, dx):

	if (j == 0):
		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j + 1) * nbComp + 0]
		q_r     = Uk[(j + 1) * nbComp + 1]
		rho_E_r = Uk[(j + 1) * nbComp + 2]	

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
		divMat.addValue(j * nbComp,       j * nbComp, Am * (-1.))
	
		p_inlet = rho_to_p_StiffenedGaz(Uk[j * nbComp + 0], Uk[j * nbComp + 1], Uk[j * nbComp + 2])# We take p from inside the domain
		rho_l=p_to_rho_StiffenedGaz(p_inlet, T_inlet) # rho is computed from the temperature BC and the inner pressure
		#rho_l   = Uk[j * nbComp + 0]                            # We take rho from inside the domain
		q_l     = rho_l * v_inlet                               # q is imposed by the boundary condition v_inlet
		rho_E_l = T_to_rhoE_StiffenedGaz(T_inlet, rho_l, q_l)   #rhoE is obtained  using the two boundary conditions v_inlet and e_inlet
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, j * nbComp, Ap)

	elif (j == nx - 1):
		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j ) * nbComp + 0]                               # We take rho inside the domain
		q_r     = Uk[(j ) * nbComp + 1]                               # We take q from inside the domain
		rho_E_r = (p_outlet+gamma*p0)/(gamma-1) + 0.5*q_r**2/rho_r    # rhoE is obtained using the boundary condition p_outlet

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)	
		divMat.addValue(j * nbComp, j * nbComp, Am * (-1.))
	
		rho_l   = Uk[(j - 1) * nbComp + 0]
		q_l     = Uk[(j - 1) * nbComp + 1]
		rho_E_l = Uk[(j - 1) * nbComp + 2]
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
		divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))


def FillMatrixFromInnerCell(j, Uk, nbComp, divMat, dt, dx):

	rho_l   = Uk[(j - 1) * nbComp + 0]
	q_l     = Uk[(j - 1) * nbComp + 1]
	rho_E_l = Uk[(j - 1) * nbComp + 2]
	rho_r   = Uk[j * nbComp + 0]
	q_r     = Uk[j * nbComp + 1]
	rho_E_r = Uk[j * nbComp + 2]
	Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	rho_l   = Uk[j * nbComp + 0]
	q_l     = Uk[j * nbComp + 1]
	rho_E_l = Uk[j * nbComp + 2]
	rho_r   = Uk[(j + 1) * nbComp + 0]
	q_r     = Uk[(j + 1) * nbComp + 1]
	rho_E_r = Uk[(j + 1) * nbComp + 2]
	Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)

	divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
	divMat.addValue(j * nbComp, j * nbComp, Am * (-1.))
	divMat.addValue(j * nbComp, j * nbComp, Ap)
	divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))
			
def FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx, isImplicit):
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

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)

		p_inlet = rho_to_p_StiffenedGaz(Uk[j * nbComp + 0], Uk[j * nbComp + 1], Uk[j * nbComp + 2])# We take p from inside the domain
		rho_l=p_to_rho_StiffenedGaz(p_inlet, T_inlet) # rho is computed from the temperature BC and the inner pressure
		#rho_l   = Uk[j * nbComp + 0]                            # We take rho from inside the domain
		q_l     = rho_l * v_inlet                               # q is imposed by the boundary condition v_inlet
		rho_E_l = T_to_rhoE_StiffenedGaz(T_inlet, rho_l, q_l)   #rhoE is obtained  using the two boundary conditions v_inlet and e_inlet
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
		if(isImplicit):
			dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
			dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
			dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
			temp1 = Am * dUi1
		
			dUi2[0] = rho_l   -  Uk[(j ) * nbComp + 0]
			dUi2[1] = q_l     -  Uk[(j ) * nbComp + 1]
			dUi2[2] = rho_E_l -  Uk[(j ) * nbComp + 2]
			temp2 = Ap * dUi2
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

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)

		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j ) * nbComp + 0]                               # We take rho inside the domain
		q_r     = Uk[(j ) * nbComp + 1]                               # We take q from inside the domain
		rho_E_r = (p_outlet+gamma*p0)/(gamma-1) + 0.5*q_r**2/rho_r    # rhoE is obtained using the boundary condition p_outlet

		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)	

		if(isImplicit):
			dUi1[0] = rho_r   - Uk[j * nbComp + 0]
			dUi1[1] = q_r     - Uk[j * nbComp + 1]
			dUi1[2] = rho_E_r - Uk[j * nbComp + 2]
			temp1 = Am * dUi1
	
			dUi2[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
			dUi2[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
			dUi2[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]
			temp2 = Ap * dUi2
		else:
			dUi1[0] = rho_r   
			dUi1[1] = q_r     
			dUi1[2] = rho_E_r 
			temp1 = Am * dUi1
	
	if(isImplicit):
		Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
		Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
		Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])
	else:#explicit scheme, contribution from the boundary data the right hand side
		Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] 
		Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] 
		Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] 

def FillRHSFromInnerCell(j, Uk, nbComp, Rhs, Un, dt, dx, isImplicit):

	rho_l   = Uk[(j - 1) * nbComp + 0]
	q_l     = Uk[(j - 1) * nbComp + 1]
	rho_E_l = Uk[(j - 1) * nbComp + 2]
	rho_r   = Uk[j * nbComp + 0]
	q_r     = Uk[j * nbComp + 1]
	rho_E_r = Uk[j * nbComp + 2]
	Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	rho_l   = Uk[j * nbComp + 0]
	q_l     = Uk[j * nbComp + 1]
	rho_E_l = Uk[j * nbComp + 2]
	rho_r   = Uk[(j + 1) * nbComp + 0]
	q_r     = Uk[(j + 1) * nbComp + 1]
	rho_E_r = Uk[(j + 1) * nbComp + 2]
	Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)

	if(isImplicit):#Contribution to the right hand side if te scheme is implicit
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


def computeSystemMatrix(a,b,nx, cfl, Uk, isImplicit):
	dim = 1
	nbComp = 3
	dx = (b - a) / nx
	dt = cfl * dx / c0
	nbVoisinsMax = 2

	nbCells = nx
	divMat = cdmath.SparseMatrixPetsc(nbCells * nbComp, nbCells * nbComp, (nbVoisinsMax + 1) * nbComp)

	divMat.zeroEntries()  #sets the matrix coefficients to zero
	for j in range(nbCells):
		
		# traitements des bords
		if (j == 0):
			FillMatrixFromEdges(j, Uk, nbComp, divMat, dt, dx)
		elif (j == nbCells - 1):
			FillMatrixFromEdges(j, Uk, nbComp, divMat, dt, dx)
		# traitement des cellules internes
		else:
			FillMatrixFromInnerCell(j, Uk, nbComp, divMat, dt, dx)
	
	if(isImplicit):		
		divMat.diagonalShift(1.)  # add one on the diagonal
	else:
		divMat*=-1.
		divMat.diagonalShift(1.)  # add one on the diagonal

	return divMat

def computeRHSVector(a,b,nx, cfl, Uk, Un, isImplicit):
	dim = 1
	nbComp = 3
	dx = (b - a) / nx
	dt = cfl * dx / c0
	nbVoisinsMax = 2

	nbCells = nx
	Rhs = cdmath.Vector(nbCells * (nbComp))

	for j in range(nbCells):
		
		# traitements des bords
		if (j == 0):
			FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx, isImplicit)
		elif (j == nbCells - 1):
			FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx, isImplicit)
		# traitement des cellules internes
		else:
			FillRHSFromInnerCell(j, Uk, nbComp, Rhs, Un, dt, dx, isImplicit)
			
	return Rhs


if __name__ == """__main__""":
	nbComp=3 # number of equations 
	a = 0.# domain is interval [a,b]
	b = 4.2# domain is interval [a,b]
	nx = 10# number of cells
	dx = (b - a) / nx  # space step
	x = [a + 0.5 * dx + i * dx for i in range(nx)]  # array of cell center (1D mesh)
	state_law = "Hermite575K"
	state_law_parameters(state_law)
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
	print("\n Testing function computeSystemMatrix \n")
	cfl = 0.5
	computeSystemMatrix(a, b, nx, cfl, U,True)  #Implicit matrix
	computeSystemMatrix(a, b, nx, cfl, U,False) #Explicit matrix

	print("\n Testing function computeRHSVector \n")
	cfl = 0.5
	computeRHSVector(a, b, nx, cfl, U, U,True)  #Implicit RHS
	computeRHSVector(a, b, nx, cfl, U, U,False) #Explicit RHS

	print("\n Testing function solve (Implicit scheme) \n")
	isImplicit=True
	cfl = 1000.
	solve(a, b, nx, "RegularSquares", "", cfl, state_law, isImplicit)
	
	print("\n Testing function solve (Explicit scheme) \n")
	isImplicit=False
	cfl = 0.5
	solve(a, b, nx, "RegularSquares", "", cfl, state_law, isImplicit)

