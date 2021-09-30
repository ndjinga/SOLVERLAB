#!/usr/bin/env python3
# -*-coding:utf-8 -*



"""
Created on Mon Apr 26 2021
@author: Michael NDJINGA

Toro 4 shock tube benchmarck from book 'Riemann Solvers and Numerical Methods for Fluid Dynamics, Third edition, Springer, 2009' (chapter 10)

Euler system in one dimension on domain [a,b]
Riemann problemn with ghost cell boundary condition
Stag scheme compared with exact solution
Regular mesh


"""

import cdmath
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from math import sqrt
from numpy import sign

#initial parameters for Riemann problem (p in Pa, v in m/s, rho in Kg/m**3)

p_L = 460.894
p_R = 46.095
v_L = 19.5975
v_R = -6.19633
rho_L = 5.99924
rho_R = 5.99242	

pos_disc = 0.4

precision = 1.e-5

gamma = 1.4
p0 = 0.
q = 0

tmax=0.035
lam_max= max(abs(v_L)+sqrt(gamma*p_L/rho_L), abs(v_R)+sqrt(gamma*p_R/rho_R) )

def initial_conditions_Riemann_problem(a, b, nx):
	print("Initial data Riemann problem")
	dx = (b - a) / nx  # space step
	x = [a + 0.5 * dx + i * dx for i in range(nx)]  # array of cell center (1D mesh)
	

	p_initial   = np.array([ (xi < pos_disc) * p_L   + (xi >= pos_disc) * p_R   for xi in x])
	v_initial   = np.array([ (xi < pos_disc) * v_L   + (xi >= pos_disc) * v_R   for xi in x])
	rho_initial = np.array([ (xi < pos_disc) * rho_L + (xi >= pos_disc) * rho_R for xi in x])

	e_initial = p_to_e_StiffenedGaz(p_initial, rho_initial)
	q_initial = rho_initial * v_initial
	rho_E_initial = rho_initial*e_initial + 0.5 * q_initial**2/rho_initial

	return rho_initial, q_initial, rho_E_initial, p_initial, v_initial, e_initial

def rho_to_p_StiffenedGaz(rho_field, q_field, rho_E_field):
	p_field = (gamma - 1.) * ( rho_E_field - 1. / 2 * q_field ** 2 / rho_field - rho_field * q) - gamma * p0
	return p_field	

def p_to_e_StiffenedGaz(p_field, rho_field):
	e_field = (p_field + gamma*p0) / (gamma - 1.) / rho_field + q
	return e_field
	
#Does not involve stiffened gas
def e_to_rhoE_StiffenedGaz(e_field, rho_field, q_field):
	rho_E_field = 1 / 2 * (q_field) ** 2 / rho_field + rho_field * e_field
	return rho_E_field

#Does not involve stiffened gas		
def rhoE_to_e_StiffenedGaz(rho_field, q_field, rho_E_field):
	e_field = rho_E_field / rho_field - 1. / 2. * (q_field / rho_field)**2
	return e_field

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
	# Dmac[0, 0] = 0
	# Dmac[0, 1] = 1
	# Dmac[0, 2] = 0
	# Dmac[1, 0] = -dp_drho - u ** 2 - dp_de / rho * (u**2/2 - e)
	# Dmac[1, 1] = 2 * u + u * dp_de / rho
	# Dmac[1, 2] = -dp_de / rho
	# Dmac[2, 0] = -u * ( dp_drho + H + dp_de / rho * (u**2/2 - e) )
	# Dmac[2, 1] = H + dp_de / rho * u ** 2
	# Dmac[2, 2] = (-dp_de / rho +1) * u
	
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
	

def jacobianMatricesm(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme):

	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	
	if scheme == 'Stag':
		Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		return (RoeMat - Dmac) * coeff * 0.5
	
	elif scheme == 'Roe':
		Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    
		return (RoeMat - Droe) * coeff * 0.5


def jacobianMatricesp(coeff, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme):
	if rho_l < 0 or rho_r < 0:
		print("rho_l=", rho_l, " rho_r= ", rho_r)
		raise ValueError("Negative density")
	if rho_E_l < 0 or rho_E_r < 0:
		print("rho_E_l=", rho_E_l, " rho_E_r= ", rho_E_r)
		raise ValueError("Negative total energy")

	RoeMat = MatRoe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
	if scheme == 'Stag':
		Dmac = Dmac_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)
		return (RoeMat + Dmac) * coeff * 0.5
	
	elif scheme == 'Roe':
		Droe = Droe_StiffenedGaz( rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r)    
		return (RoeMat + Droe) * coeff * 0.5


def FillEdges(j, Uk, nbComp, divMat, Rhs, Un, dt, dx, scheme):
	if (j == 0):
		rho_l   = Uk[j * nbComp + 0]
		q_l     = Uk[j * nbComp + 1]
		rho_E_l = Uk[j * nbComp + 2]
		rho_r   = Uk[(j + 1) * nbComp + 0]
		q_r     = Uk[(j + 1) * nbComp + 1]
		rho_E_r = Uk[(j + 1) * nbComp + 2]
		
		Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme)
		divMat.addValue(j * nbComp, (j + 1) * nbComp, Am)
		divMat.addValue(j * nbComp,       j * nbComp, Am * (-1.))

		dUi = cdmath.Vector(3)
		dUi[0] = Uk[(j + 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi[1] = Uk[(j + 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi[2] = Uk[(j + 1) * nbComp + 2] - Uk[j * nbComp + 2]
		temp = Am * dUi

		Rhs[j * nbComp + 0] = -temp[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
		Rhs[j * nbComp + 1] = -temp[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
		Rhs[j * nbComp + 2] = -temp[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])

	elif (j == Un.size()/nbComp - 1):
		rho_l   = Uk[(j - 1) * nbComp + 0]
		q_l     = Uk[(j - 1) * nbComp + 1]
		rho_E_l = Uk[(j - 1) * nbComp + 2]
		rho_r   = Uk[j * nbComp + 0]
		q_r     = Uk[j * nbComp + 1]
		rho_E_r = Uk[j * nbComp + 2]

		Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme)
		divMat.addValue(j * nbComp, j * nbComp, Ap)
		divMat.addValue(j * nbComp, (j - 1) * nbComp, Ap * (-1.))

		dUi = cdmath.Vector(3)
		dUi[0] = Uk[(j - 1) * nbComp + 0] - Uk[j * nbComp + 0]
		dUi[1] = Uk[(j - 1) * nbComp + 1] - Uk[j * nbComp + 1]
		dUi[2] = Uk[(j - 1) * nbComp + 2] - Uk[j * nbComp + 2]

		temp = Ap * dUi
		Rhs[j * nbComp + 0] = temp[0] - (Uk[j * nbComp + 0] - Un[j * nbComp + 0])
		Rhs[j * nbComp + 1] = temp[1] - (Uk[j * nbComp + 1] - Un[j * nbComp + 1])
		Rhs[j * nbComp + 2] = temp[2] - (Uk[j * nbComp + 2] - Un[j * nbComp + 2])

def FillInnerCell(j, Uk, nbComp, divMat, Rhs, Un, dt, dx, scheme):

	rho_l   = Uk[(j - 1) * nbComp + 0]
	q_l     = Uk[(j - 1) * nbComp + 1]
	rho_E_l = Uk[(j - 1) * nbComp + 2]
	rho_r   = Uk[j * nbComp + 0]
	q_r     = Uk[j * nbComp + 1]
	rho_E_r = Uk[j * nbComp + 2]
	Ap = jacobianMatricesp(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme)
	
	rho_l   = Uk[j * nbComp + 0]
	q_l     = Uk[j * nbComp + 1]
	rho_E_l = Uk[j * nbComp + 2]
	rho_r   = Uk[(j + 1) * nbComp + 0]
	q_r     = Uk[(j + 1) * nbComp + 1]
	rho_E_r = Uk[(j + 1) * nbComp + 2]
	Am = jacobianMatricesm(dt / dx, rho_l, q_l, rho_E_l, rho_r, q_r, rho_E_r, scheme)

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

def SetPicture():
	max_initial_rho = max(rho_L,rho_R)
	min_initial_rho = min(rho_L,rho_R)
	max_initial_p = max(p_L,p_R)
	min_initial_p = min(p_L,p_R)
	max_initial_v = max(v_L,v_R)
	min_initial_v = min(v_L,v_R)
	max_initial_q = max_initial_rho*max_initial_v
	min_initial_q = min_initial_rho*min_initial_v

	e_L=p_to_e_StiffenedGaz(p_L, rho_L)
	e_R=p_to_e_StiffenedGaz(p_R, rho_R)
	h_L=e_L+p_L/rho_L
	h_R=e_R+p_R/rho_R
	max_initial_e = max(e_L, e_R)
	min_initial_e = min(e_L, e_R)
	min_initial_h = min(h_L,h_R)
	max_initial_h = max(h_L,h_R)

	fig, ([axDensity, axMomentum, axEnthalpie],[axPressure, axVitesse, axEinterne]) = plt.subplots(2, 3,sharex=True, figsize=(14,10))
	#fig.suptitle('Staggered scheme')
	plt.gcf().subplots_adjust(wspace = 0.5)

	axDensity.set(xlabel='x (m)', ylabel='Densité (Kg/m3)')
	axDensity.set_xlim(a,b)
	axDensity.set_ylim(0.9*min_initial_rho , 6*max_initial_rho )
	#axDensity.set_ylim(0.125, 1.01)
	axDensity.legend()

	axMomentum.set(xlabel='x (m)', ylabel='Momentum (Kg/m2/s)')
	axMomentum.set_xlim(a,b)
	axMomentum.set_ylim(0.9*min_initial_q , 2.5*max_initial_q )
	#axMomentum.set_ylim(0., 1.)
	axMomentum.legend()
	
	axEnthalpie.set(xlabel='x (m)', ylabel='h (J/Kg)')
	axEnthalpie.set_xlim(a,b)
	axEnthalpie.set_ylim(0.9*min_initial_h , 1.75*max_initial_h )
	#axEnthalpie.set_ylim(2.5, 5.)
	axEnthalpie.legend()
	
	axPressure.set(xlabel='x (m)', ylabel='Pression (Pa)')
	axPressure.set_xlim(a,b)
	axPressure.set_ylim(0.9*min_initial_p , 4*max_initial_p)
	#axPressure.set_ylim(0.1, 1)
	axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axPressure.legend()

	axVitesse.set(xlabel='x (m)', ylabel='Vitesse (m/s)')
	axVitesse.set_xlim(a,b)
	axVitesse.set_ylim(0.9*min_initial_v , 1.1*max_initial_v)
	#axVitesse.set_ylim(0., 1.5)
	axVitesse.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axVitesse.legend()

	axEinterne.set(xlabel='x (m)', ylabel='Energie interne (J/Kg)')
	axEinterne.set_xlim(a,b)
	axEinterne.set_ylim(0.9*min_initial_e , 1.75*max_initial_e)
	#axEinterne.set_ylim(2., 3.5)
	axEinterne.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
	axEinterne.legend()
	
	return(fig)

def addCurves(rho_field_Stag, q_field_Stag, h_field_Stag, p_field_Stag, v_field_Stag, e_field_Stag, dx, nbCells, fig):

	[axDensity, axMomentum, axEnthalpie,axPressure, axVitesse, axEinterne] = fig.get_axes()

	lineDensity_Stag, = axDensity.plot([a+0.5*dx + i*dx for i in range(nbCells)], rho_field_Stag, label='Stag on '+str(nbCells)+' cells') #new picture for video # Returns a tuple of line objects, thus the comma

	lineMomentum_Stag, = axMomentum.plot([a+0.5*dx + i*dx for i in range(nbCells)], q_field_Stag, label='Stag on '+str(nbCells)+' cells')
	
	lineEnthalpie_Stag, = axEnthalpie.plot([a+0.5*dx + i*dx for i in range(nbCells)], h_field_Stag, label='Stag on '+str(nbCells)+' cells')
	
	linePressure_Stag, = axPressure.plot([a+0.5*dx + i*dx for i in range(nbCells)], p_field_Stag, label='Stag on '+str(nbCells)+' cells')

	lineVitesse_Stag, = axVitesse.plot([a+0.5*dx + i*dx for i in range(nbCells)], v_field_Stag, label='Stag on '+str(nbCells)+' cells')

	lineEinterne_Stag, = axEinterne.plot([a+0.5*dx + i*dx for i in range(nbCells)], e_field_Stag, label='Stag on '+str(nbCells)+' cells')
	
	return

def addSolutionToro(fig) :

	x_exact, rho_exact, v_exact, p_exact, e_exact, h_exact, Mach_exact, left_or_right = np.loadtxt(os.path.dirname(__file__)+'/TTC'+str(4) + '.dat', unpack = True )

	[axDensity, axMomentum, axEnthalpie,axPressure, axVitesse, axEinterne] = fig.get_axes()
	lineRhoToro, = axDensity.plot(x_exact, rho_exact , label='Solution exacte')
	axDensity.legend()
	lineVToro, = axVitesse.plot(x_exact, v_exact, label='Solution exacte')
	axVitesse.legend()
	lineQToro, = axMomentum.plot(x_exact, v_exact*rho_exact, label='Solution exacte')
	axMomentum.legend()
	lineEToro, = axEinterne.plot(x_exact, e_exact, label='Solution exacte')
	axEinterne.legend()
	linePToro, = axPressure.plot(x_exact, p_exact, label='Solution exacte')
	axPressure.legend()
	lineEnthalpieToro, = axEnthalpie.plot(x_exact, h_exact, label='Solution exacte')
	axEnthalpie.legend()


def EulerSystemStaggered(ntmax, tmax, cfl, a, b, nbCells, output_freq, meshName,fig):
	nbComp = 3
	time = 0.
	it = 0
	isStationary = False
	dx = (b - a) / nbCells
	dt = cfl * dx / lam_max
	nbVoisinsMax = 2

	# iteration vectors
	Un_Stag  = cdmath.Vector(nbCells * (nbComp))
	dUn_Stag = cdmath.Vector(nbCells * (nbComp))
	dUk_Stag = cdmath.Vector(nbCells * (nbComp))
	Rhs_Stag = cdmath.Vector(nbCells * (nbComp))
	
	# Initial conditions
	print("Construction of the initial condition …")
	rho_field_Stag, q_field_Stag, rho_E_field_Stag, p_field_Stag, v_field_Stag, e_field_Stag = initial_conditions_Riemann_problem(a, b, nbCells)
	h_field_Stag = (rho_E_field_Stag + p_field_Stag) / rho_field_Stag - 0.5 * (q_field_Stag / rho_field_Stag) **2
	

	for k in range(nbCells):
		Un_Stag[k * nbComp + 0] = rho_field_Stag[k]
		Un_Stag[k * nbComp + 1] = q_field_Stag[k]
		Un_Stag[k * nbComp + 2] = rho_E_field_Stag[k]
		
	divMat_Stag = cdmath.SparseMatrixPetsc(nbCells * nbComp, nbCells * nbComp, (nbVoisinsMax + 1) * nbComp)
		
	iterGMRESMax = 50
	newton_max = 100

	print("Starting computation of the non linear Euler non isentropic system with staggered scheme …")
	# STARTING TIME LOOP
	while (it < ntmax and time <= tmax and not isStationary):
		dUn_Stag = Un_Stag.deepCopy()
		Uk_Stag  = Un_Stag.deepCopy()
		residu_Stag = 1.
		
		k_Stag = 0
		while (k_Stag < newton_max and residu_Stag > precision):
			# STARTING NEWTON LOOP
			#print("Iteration k=", k_Stag, " residu = ",residu_Stag)
			divMat_Stag.zeroEntries()  #sets the matrix coefficients to zero
			for j in range(nbCells):
				# traitements des bords
				if (j == 0):
					FillEdges(j, Uk_Stag, nbComp, divMat_Stag, Rhs_Stag, Un_Stag, dt, dx, 'Stag')
					
				elif (j == nbCells - 1):
					FillEdges(j, Uk_Stag, nbComp, divMat_Stag, Rhs_Stag, Un_Stag, dt, dx, 'Stag')

				# traitement des cellules internes
				else:
					FillInnerCell(j, Uk_Stag, nbComp, divMat_Stag, Rhs_Stag, Un_Stag, dt, dx, 'Stag')
					
			#print(divMat_Stag)
			divMat_Stag.diagonalShift(1)  # only after  filling all coefficients

			#solving the linear system divMat * dUk = Rhs
			LS_Stag = cdmath.LinearSolver(divMat_Stag, Rhs_Stag, iterGMRESMax, precision, "GMRES", "LU")
			dUk_Stag = LS_Stag.solve()
			residu_Stag = dUk_Stag.norm()
			
			#updates for Newton loop
			Uk_Stag += dUk_Stag
			k_Stag = k_Stag + 1
			if (not LS_Stag.getStatus()):
				print("Linear system did not converge ", LS_Stag.getNumberOfIter(), " GMRES iterations")
				raise ValueError("No convergence of the linear system")
				
			if k_Stag == newton_max:
				raise ValueError("No convergence of Newton Staggered Scheme")
		

		#updating fields
		Un_Stag = Uk_Stag.deepCopy()
		dUn_Stag -= Un_Stag
		
		if (dUn_Stag.norm()<precision):
				isStationary = True
		
		time = time + dt
		it = it + 1

		# Savings
		if (it == 1 or it % output_freq == 0 or it >= ntmax or isStationary or time >= tmax):
	
			print("-- Time step : " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
			print("Linear system converged in ", LS_Stag.getNumberOfIter(), " GMRES iterations")

			plt.savefig("EulerCompletStaggeredToro4" + meshName + str(it) + '_time' + str(time) + ".png")

	print("##################### Last time step: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

	# Postprocessing
	for k in range(nbCells):
		rho_field_Stag[k]   = Un_Stag[k * nbComp + 0]
		q_field_Stag[k]     = Un_Stag[k * nbComp + 1]
		rho_E_field_Stag[k] = Un_Stag[k * nbComp + 2]
		
	v_field_Stag = q_field_Stag / rho_field_Stag
	p_field_Stag = rho_to_p_StiffenedGaz(rho_field_Stag, q_field_Stag, rho_E_field_Stag)
	e_field_Stag = rhoE_to_e_StiffenedGaz(rho_field_Stag, q_field_Stag, rho_E_field_Stag)
	h_field_Stag = (rho_E_field_Stag + p_field_Stag) / rho_field_Stag - 0.5 * (q_field_Stag / rho_field_Stag) **2
				
	# Picture settings
	addCurves(rho_field_Stag, q_field_Stag, h_field_Stag, p_field_Stag, v_field_Stag, e_field_Stag, dx, nbCells,fig)

	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_rho_" + str(nbCells) + "_cells.txt", rho_field_Stag,    delimiter="\n")
	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_q_"   + str(nbCells) + "_cells.txt",   q_field_Stag,    delimiter="\n")
	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_h_"   + str(nbCells) + "._cellstxt", h_field_Stag, delimiter="\n")
	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_p_"    + str(nbCells) + "_cells.txt",     p_field_Stag, delimiter="\n")
	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_v_"    + str(nbCells) + "_cells.txt",     v_field_Stag, delimiter="\n")
	np.savetxt("EulerCompletStaggeredToro4" + meshName + "_e_"    + str(nbCells) + "_cells.txt",     e_field_Stag, delimiter="\n")
			
	if (it >= ntmax):
		print("Maximum number of time steps ntmax= ", ntmax, " reached on "+ str(nbCells) + " cells")

	elif (time >= tmax):
		print("Maximum time tmax= ", tmax, " reached on "+ str(nbCells) + " cells")
	elif (isStationary):
		print("Stationary regime reached at time step ", it, ", t= ", time,  ", on "+ str(nbCells) + " cells")
	else:
		raise ValueError("Unexpected end of simulation")

	return

def solve(a, b, meshName, meshType, cfl):
	print("Convergence study for the Euler System in dimension 1")
	print("State Law Stiffened Gaz, gamma = ", gamma, " p0= ", p0 )
	print("Initial data : ", "Riemann problem")
	print("Boundary conditions : ", "Neumann")
	print("Mesh name : ", meshName)
	# Problem data
	ntmax = 1000000
	output_freq = 100
	fig=SetPicture()
	for nx in [20, 50, 100]:
		print("Resolution of the Euler System in dimension 1 on " + str(nx) + " cells")
		EulerSystemStaggered(ntmax, tmax, cfl, a, b, nx, output_freq, meshName,fig)

	addSolutionToro(fig)
	plt.savefig("EulerCompletStaggeredToro4" + meshName + "_convergence.png")

	return


if __name__ == """__main__""":
	a = 0.
	b = 1.
	cfl = 0.99
	solve(a, b, "RegularGrid", "", cfl)
