#!/usr/bin/env python3
# -*-coding:utf-8 -*

"""
Created on Mon Sep 27 2021
@author: Michael NDJINGA, Katia Ait Ameur, Coraline Mounier

Euler system with heating source term (phi) in one dimension on regular domain [a,b]
Riemann problemn with ghost cell boundary condition
Left : Inlet boundary condition (velocity and temperature imposed)
Right : Outlet boundary condition (pressure imposed)
Staggered scheme
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
from math import sqrt
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

	
def Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):
	Dmac   = cdmath.Matrix(3, 3)

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

	p = rho_h_to_p_StiffenedGaz( rho, H - u**2/2. )
	e = H - p / rho - 1./2 * u**2
	dp_drho = dp_drho_e_const_StiffenedGaz( e )
	dp_de   = dp_de_rho_const_StiffenedGaz( rho )
 
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
	
	return Dmac * sign(u)
	
def jacobianMatricesm(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):

	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    
	return (RoeMat - Dmac) * coeff * 0.5


def jacobianMatricesp(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r):
	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    

	return (RoeMat + Dmac) * coeff * 0.5


def FillEdges(j, Uk, nbComp, divMat, Rhs, Un, dt, dx):
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
		q_l     = rho_l * v_inlet                               # q is imposed by the boundary condition v_inlet
		rho_E_l = T_to_rhoE_StiffenedGaz(T_inlet, rho_l, q_l)   #rhoE is obtained  using the two boundary conditions v_inlet and e_inlet
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
	
		dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
		temp1 = Am * dUi1
	
		dUi2[0] = rho_l   -  Uk[(j ) * nbComp + 0]
		dUi2[1] = q_l     -  Uk[(j ) * nbComp + 1]
		dUi2[2] = rho_E_l -  Uk[(j ) * nbComp + 2]
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

		dUi1[0] = rho_r   - Uk[j * nbComp + 0]
		dUi1[1] = q_r     - Uk[j * nbComp + 1]
		dUi1[2] = rho_E_r - Uk[j * nbComp + 2]
		temp1 = Am * dUi1

		dUi2[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi2[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi2[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]
		temp2 = Ap * dUi2
	
	Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
	Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
	Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])

def FillInnerCell(j, Uk, nbComp, divMat, Rhs, Un, dt, dx):

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

def SetPicture(rho_field, q_field, h_field, p_field, v_field, T_field, dx):
	max_initial_rho = max(rho_field)
	min_initial_rho = min(rho_field)
	max_initial_q = max(q_field)
	min_initial_q = min(q_field)
	min_initial_h = min(h_field)
	max_initial_h = max(h_field)
	max_initial_p = max(p_field)
	min_initial_p = min(p_field)
	max_initial_v = max(v_field)
	min_initial_v = min(v_field)
	max_initial_T = max(T_field)
	min_initial_T = min(T_field)

	fig, ([axDensity, axMomentum, axh],[axPressure, axVitesse, axTemperature]) = plt.subplots(2, 3,sharex=True, figsize=(14,10))
	plt.gcf().subplots_adjust(wspace = 0.5)

	lineDensity, = axDensity.plot([a+0.5*dx + i*dx for i in range(nx)], rho_field, label='MAC scheme')
	axDensity.set(xlabel='x (m)', ylabel='Densité (kg/m3)')
	axDensity.set_xlim(a,b)
	axDensity.set_ylim(680, 800)
	axDensity.legend()

	lineMomentum, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nx)], q_field, label='MAC scheme')
	axMomentum.set(xlabel='x (m)', ylabel='Momentum (kg/(m2.s))')
	axMomentum.set_xlim(a,b)
	axMomentum.set_ylim(3760, 	3780)
	axMomentum.legend()

	lineh, = axh.plot([a+0.5*dx + i*dx for i in range(nx)], h_field, label='MAC scheme')
	axh.set(xlabel='x (m)', ylabel='h (J/m3)')
	axh.set_xlim(a,b)
	axh.set_ylim(1.25 * 10**6, 1.45*10**6)
	axh.legend()
	
	linePressure, = axPressure.plot([a+0.5*dx + i*dx for i in range(nx)], p_field, label='MAC scheme')
	axPressure.set(xlabel='x (m)', ylabel='Pression (bar)')
	axPressure.set_xlim(a,b)
	axPressure.set_ylim(154.99, 155.02)
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axPressure.legend()

	lineVitesse, = axVitesse.plot([a+0.5*dx + i*dx for i in range(nx)], v_field, label='MAC scheme')
	axVitesse.set(xlabel='x (m)', ylabel='Vitesse (m/s)')
	axVitesse.set_xlim(a,b)
	axVitesse.set_ylim(v_0-1, v_0+1)
	axVitesse.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axVitesse.legend()

	lineTemperature, = axTemperature.plot([a+0.5*dx + i*dx for i in range(nx)], T_field, label='MAC scheme')
	axTemperature.set(xlabel='x (m)', ylabel='Température (K)')
	axTemperature.set_xlim(a,b)
	axTemperature.set_ylim(T_0-10, T_0+30)
	axTemperature.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axTemperature.legend()
	
	return(fig, lineDensity, lineMomentum, lineh, linePressure, lineVitesse, lineTemperature)


def EulerSystemMAC(ntmax, tmax, cfl, a, b, nbCells, output_freq, meshName, state_law):
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
	Un  = cdmath.Vector(nbCells * (nbComp))
	dUn = cdmath.Vector(nbCells * (nbComp))
	dUk = cdmath.Vector(nbCells * (nbComp))
	Rhs = cdmath.Vector(nbCells * (nbComp))

	# Initial conditions
	print("Construction of the initial condition …")

	rho_field, q_field, rho_E_field, p_field, v_field, T_field = initial_conditions_Riemann_problem(a, b, nx)
	h_field = (rho_E_field + p_field) / rho_field - 0.5 * (q_field / rho_field) **2
	p_field = p_field * 10 ** (-5)
	

	for k in range(nbCells):
		Un[k * nbComp + 0] = rho_field[k]
		Un[k * nbComp + 1] = q_field[k]
		Un[k * nbComp + 2] = rho_E_field[k]

	divMat = cdmath.SparseMatrixPetsc(nbCells * nbComp, nbCells * nbComp, (nbVoisinsMax + 1) * nbComp)
	
	# Picture settings
	fig, lineDensity, lineMomentum, lineRhoE, linePressure, lineVitesse, lineTemperature = SetPicture( rho_field, q_field, h_field, p_field, v_field, T_field, dx)

	plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + "0" + ".png")
	iterGMRESMax = 50
	newton_max = 100

	print("Starting computation of the non linear Euler non isentropic system with MAC scheme …")
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
					FillEdges(j, Uk, nbComp, divMat, Rhs, Un, dt, dx)
				elif (j == nbCells - 1):
					FillEdges(j, Uk, nbComp, divMat, Rhs, Un, dt, dx)

				# traitement des cellules internes
				else:
					FillInnerCell(j, Uk, nbComp, divMat, Rhs, Un, dt, dx)
					
				Rhs[j * nbComp + 2] += phi*dt
			
			#solving the linear system divMat * dUk = Rhs
			divMat.diagonalShift(1.)
			LS = cdmath.LinearSolver(divMat, Rhs, iterGMRESMax, precision, "GMRES", "LU")
			dUk = LS.solve()
			vector_residu = dUk.maxVector(nbComp)
			residu = max(abs(vector_residu[0])/rho0, abs(vector_residu[1])/(rho0*v_0), abs(vector_residu[2])/rhoE0 )
			
			if (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax):
				print("Residu Newton at iteration ",k, " :   ", residu)
				print("Linear system converged in ", LS.getNumberOfIter(), " GMRES iterations")

			#updates for Newton loop
			Uk += dUk
			k = k + 1
			if (not LS.getStatus()):
				print("Linear system did not converge ", LS.getNumberOfIter(), " GMRES iterations")
				raise ValueError("No convergence of the linear system")
			
			if k == newton_max:
				raise ValueError("No convergence of Newton MAC Scheme")

		#updating fields
		Un = Uk.deepCopy()
		dUn -= Un

		#Testing stationarity
		residu_stat = dUn.maxVector(nbComp)#On prend le max de chaque composante
		if (it % output_freq == 0 ):
			print("Test de stationarité : Un+1-Un= ", max(abs(residu_stat[0])/rho0, abs(residu_stat[1])/(rho0*v_0), abs(residu_stat[2])/rhoE0 ))

		if ( it>1 and abs(residu_stat[0])/rho0<precision  and abs(residu_stat[1])/(rho0*v_0)<precision and abs(residu_stat[2])/rhoE0<precision):
				isStationary = True
		
		for k in range(nbCells):
			rho_field[k]   = Un[k * nbComp + 0]
			q_field[k]     = Un[k * nbComp + 1]
			rho_E_field[k] = Un[k * nbComp + 2]

		v_field = q_field / rho_field
		p_field = rho_to_p_StiffenedGaz(rho_field, q_field, rho_E_field)
		T_field = rhoE_to_T_StiffenedGaz(rho_field, q_field, rho_E_field)
		h_field = (rho_E_field + p_field) / rho_field - 0.5 * (q_field / rho_field) **2
		p_field = p_field * 10 ** (-5)
		
		if( min(p_field)<0) :
			raise ValueError("Negative pressure, stopping calculation")

		#picture and video updates
		lineDensity.set_ydata(rho_field)
		lineMomentum.set_ydata(q_field)
		lineRhoE.set_ydata(h_field)
		linePressure.set_ydata(p_field)
		lineVitesse.set_ydata(v_field)
		lineTemperature.set_ydata(T_field)
		
		time = time + dt
		it = it + 1

		# Savings
		if (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax):
	
			print("-- Time step : " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

			print("Temperature gain between inlet and outlet is ", T_field[nbCells-1]-T_field[0],"\n")

			plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + str(it) + '_time' + str(time) + ".png")

	print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)+"\n")
	
	if (it >= ntmax):
		print("Maximum number of time steps ntmax= ", ntmax, " reached")
		return

	elif (isStationary):
		print("Stationary regime reached at time step ", it, ", t= ", time)
		print("------------------------------------------------------------------------------------")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + "_rho_Stat.txt", rho_field, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + "_q_Stat.txt", q_field, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + "_rhoE_Stat.txt", rho_E_field, delimiter="\n")
		np.savetxt("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName + "_p_Stat.txt", p_field, delimiter="\n")
		plt.savefig("EulerComplet_HeatedChannel_" + str(dim) + "D_MAC" + meshName +"_Stat.png")
		return
	else:
		print("Maximum time Tmax= ", tmax, " reached")
		return


def solve(a, b, nx, meshName, meshType, cfl, state_law):
	print("Simulation of a heated channel in dimension 1 on " + str(nx) + " cells")
	print("State Law Stiffened Gaz, " + state_law)
	print("Initial data : ", "constant fields")
	print("Boundary conditions : ", "Inlet (Left), Outlet (Right)")
	print("Mesh name : ", meshName, ", ", nx, " cells")
	# Problem data
	tmax = 10.
	ntmax = 100000
	output_freq = 1000
	EulerSystemMAC(ntmax, tmax, cfl, a, b, nx, output_freq, meshName, state_law)
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
			
def FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx):
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
	
		dUi1[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi1[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi1[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
		temp1 = Am * dUi1
	
		dUi2[0] = rho_l   -  Uk[(j ) * nbComp + 0]
		dUi2[1] = q_l     -  Uk[(j ) * nbComp + 1]
		dUi2[2] = rho_E_l -  Uk[(j ) * nbComp + 2]
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

		dUi1[0] = rho_r   - Uk[j * nbComp + 0]
		dUi1[1] = q_r     - Uk[j * nbComp + 1]
		dUi1[2] = rho_E_r - Uk[j * nbComp + 2]
		temp1 = Am * dUi1

		dUi2[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi2[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi2[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]
		temp2 = Ap * dUi2
	
	Rhs[j * nbComp + 0] = -temp1[0] + temp2[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
	Rhs[j * nbComp + 1] = -temp1[1] + temp2[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
	Rhs[j * nbComp + 2] = -temp1[2] + temp2[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])

def FillRHSFromInnerCell(j, Uk, nbComp, Rhs, Un, dt, dx):

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


def computeSystemMatrix(a,b,nx, cfl, Uk):
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
	
	divMat.diagonalShift(1.)  # add one on the diagonal

	return divMat

def computeRHSVector(a,b,nx, cfl, Uk, Un):
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
			FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx)
		elif (j == nbCells - 1):
			FillRHSFromEdges(j, Uk, nbComp, Rhs, Un, dt, dx)
		# traitement des cellules internes
		else:
			FillRHSFromInnerCell(j, Uk, nbComp, Rhs, Un, dt, dx)
			
	return Rhs


if __name__ == """__main__""":
	nbComp=3 # number of equations 
	a = 0.# domain is interval [a,b]
	b = 4.2# domain is interval [a,b]
	nx = 50# number of cells
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
	computeSystemMatrix(a, b, nx, cfl, U)

	print("\n Testing function computeRHSVector \n")
	cfl = 0.5
	computeRHSVector(a, b, nx, cfl, U, U) 

	print("\n Testing function solve \n")
	cfl = 1000.
	solve(a, b, nx, "RegularSquares", "", cfl, state_law)
	
