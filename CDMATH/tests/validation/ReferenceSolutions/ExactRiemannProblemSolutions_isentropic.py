#!/usr/bin/env python3
# -*-coding:utf-8 -*

######################################################################################################################
#	This file runs a set of Riemann problems for the one dimensional isentropic Euler equations with isentropic equation of state P=c1*rho^gamma + pinf. 
#	It uses the class exact_rs_stiffenedgas_isentropic to compute the exact solution of the Riemann Problems.
#
#	Author: Michael Ndjinga
#	Date:	15/05/2023
#   Description : Translated from C++ package developped by Murray Cutforth
#######################################################################################################################

import exact_rs_stiffenedgas_isentropic
from math import fabs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def run_Riemann_problems(numsamples = 100):
	# Output test solution for many different Riemann problems
	
	print("")
	print("Determination of the exact solutions of some Riemann problems for the isentropic Euler equations on " + str(numsamples) + " sample points.")

	for TC in range(1, 23):
		WL = [0.]*2
		WR = [0.]*2
	
		gamma = 1.4;
		pinf = 0.0;
		offset = 0.5;#position of the initial discontinuity
		xmin = 0.0;
		xmax = 1.0;


		# TC1 to TC5 are the 7 shock tube problems from Toro
	
		if (TC == 1):
			WL[0] = 1.0;
			WL[1] = 0.75;
			WR[0] = 0.125;
			WR[1] = 0.0;
			t = 0.2;
			offset = 0.3;
			filename = "1";
		elif (TC ==2):
			WL[0] = 1.0;
			WL[1] = -2.0;
			WL[2] = 0.4;
			WR[0] = 1.0;
			WR[1] = 2.0;
			WR[2] = 0.4;
			t = 0.15;
			filename = "TTC2";
		elif (TC ==3):
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1000.0;
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 0.01;
			t = 0.012;
			filename = "TTC3";
		elif (TC ==4):
			WL[0] = 5.99924;
			WL[1] = 19.5975;
			WL[2] = 460.894;
			WR[0] = 5.99242;
			WR[1] = -6.19633;
			WR[2] = 46.0950;
			t = 0.035;
			offset = 0.4;
			filename = "TTC4";
		elif (TC ==5):
			WL[0] = 1.0;
			WL[1] = -19.59745;
			WL[2] = 1000.;
			WR[0] = 1.0;
			WR[1] = -19.59745;
			WR[2] = 0.01;
			t = 0.012;
			offset = 0.8;
			filename = "TTC5";
		elif (TC == 21):
			WL[0] = 1.4;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 1.;
			WR[1] = 0.0;
			WR[2] = 1.;
			t = 2.;
			filename = "TTC6";
		elif (TC == 22):
			WL[0] = 1.4;
			WL[1] = 0.1;
			WL[2] = 1.0;
			WR[0] = 1.;
			WR[1] = 0.1;
			WR[2] = 1.;
			t = 2.;
			filename = "TTC7";
		elif (TC == 6):
			# Air - Helium shock tube from Sambasivan 2009
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 0.1;
			t = 0.25;
			filename = "Samb1";
			gammaR = 1.667;
		else:
			print( "Unknown test case" )
			return 1;
	
	
		RS = exact_rs_stiffenedgas.exact_rs_stiffenedgas(gamma, c1, pinf);
		RS.solve_RP(WL,WR);
	
		print( "")
		print( "Solved Riemann problem for test case = " , TC )
		print( "Star state pressure calculated as " , RS.P_STAR )
		print( "Star state velocity calculated as " , RS.S_STAR )
		print( "Left shock speed calculated as " , RS.S_L )
		print( "Right shock speed calculated as " , RS.S_R )
		print( "Left rarefaction head speed calculated as " , RS.S_HL )
		print( "Left rarefaction tail speed calculated as " , RS.S_TL )
		print( "Right rarefaction head speed calculated as " , RS.S_HR )
		print( "Right rarefaction tail speed calculated as " , RS.S_TR  )
	
		
		dx = (xmax - xmin)/numsamples;

		outfile=open(filename+'.dat', 'w')
	
		rho_field = [0]*(numsamples+1)
		u_field   = [0]*(numsamples+1)
		p_field   = [0]*(numsamples+1)
		q_field = [0]*(numsamples+1)

		x = xmin;
		for i in range(numsamples+1):
		
			S = x/t;
			soln = RS.sample_solution(WL, WR, S - offset/t);

			thisgamma = gamma
			thispinf = pinf
			thisz = 1.0  if S - offset/t < RS.S_STAR else 0.0;

			u_field[i]  =soln[1]
			p_field[i]  =soln[2]
			q_field[i]  =rho_field[i]*u_field[i]
			rho_field[i]=soln[0]

			outfile.write( str( x ) + " " + str( soln[0]) + " " + str( soln[1]) + " " + str( soln[2]) + " " + str(exact_rs_stiffenedgas.stiffenedgas_e(soln[0], soln[2], thisgamma, thispinf)) + " " + str(exact_rs_stiffenedgas.stiffenedgas_h(soln[0], soln[2], thisgamma, thispinf)) + " " + str(fabs(soln[1])/RS.a(soln[0], soln[2], thisgamma, thispinf)) + " " + str(thisz) + "\n")

			x += dx;

		outfile.close()

		#Set picture parameters    
		max_initial_rho = max(WL[0],WR[0])
		min_initial_rho = min(WL[0],WR[0])
		max_initial_p = max(WL[2],WR[2])
		min_initial_p = min(WL[2],WR[2])
		max_initial_v = max(WL[2],WR[2])
		min_initial_v = min(WL[2],WR[2])
		max_initial_q = max_initial_rho*max_initial_v
		min_initial_q = min_initial_rho*min_initial_v
		
		e_L=exact_rs_stiffenedgas.stiffenedgas_e(WL[0], WL[2], gammaL, pinf_L)
		e_R=exact_rs_stiffenedgas.stiffenedgas_e(WR[0], WR[2], gammaR, pinf_R)
		h_L=e_L+WL[2]/WL[0]
		h_R=e_R+WR[2]/WR[0]
		max_initial_e = max(e_L,e_R)
		min_initial_e = min(e_L,e_R)
		min_initial_h = min(h_L,h_R)
		max_initial_h = max(h_L,h_R)
		
		# Create plots
		fig, ([axDensity, axMomentum, axEnthalpie],[axPressure, axVitesse, axEinterne]) = plt.subplots(2, 3,sharex=True, figsize=(14,10))
		plt.gcf().subplots_adjust(wspace = 0.5)
		
		axDensity.plot(  [xmin + i*dx for i in range(numsamples+1)], rho_field, label='Density, '+str(numsamples)+' cells') #new picture for video # Returns a tuple of line objects, thus the comma
		axDensity.set(xlabel='x (m)', ylabel='Densité (Kg/m3)')
		axDensity.set_xlim(xmin,xmax)
		axDensity.set_ylim(0.9*min_initial_rho , 6*max_initial_rho )
		axDensity.legend()
		
		axMomentum.plot( [xmin + i*dx for i in range(numsamples+1)],  q_field, label='Momentum, '+str(numsamples)+' cells')    
		axMomentum.set(xlabel='x (m)', ylabel='Momentum (Kg/m2/s)')
		axMomentum.set_xlim(xmin,xmax)
		axMomentum.set_ylim(0.9*min_initial_q , 2.5*max_initial_q )
		axMomentum.legend()
		
		axEnthalpie.plot([xmin + i*dx for i in range(numsamples+1)],   h_field, label='Enthalpy, '+str(numsamples)+' cells')    
		axEnthalpie.set(xlabel='x (m)', ylabel='Enthalpy (J/Kg)')
		axEnthalpie.set_xlim(xmin,xmax)
		axEnthalpie.set_ylim(0.9*min_initial_h , 1.75*max_initial_h )
		axEnthalpie.legend()
		
		axPressure.plot( [xmin + i*dx for i in range(numsamples+1)],   p_field, label='Pressure, '+str(numsamples)+' cells')
		axPressure.set(xlabel='x (m)', ylabel='Pressure (Pa)')
		axPressure.set_xlim(xmin,xmax)
		axPressure.set_ylim(0.9*min_initial_p , 4*max_initial_p)
		axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
		axPressure.legend()
		
		axVitesse.plot(  [xmin + i*dx for i in range(numsamples+1)],   u_field, label=' Velocity, '+str(numsamples)+' cells')
		axVitesse.set(xlabel='x (m)', ylabel='Velocity (m/s)')
		axVitesse.set_xlim(xmin,xmax)
		axVitesse.set_ylim(0.9*min_initial_v , 1.1*max_initial_v)
		axVitesse.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
		axVitesse.legend()
		
		axEinterne.plot( [xmin + i*dx for i in range(numsamples+1)],   e_field, label='Internal energy, '+str(numsamples)+' cells')
		axEinterne.set(xlabel='x (m)', ylabel='Internal energy (J/Kg)')
		axEinterne.set_xlim(xmin,xmax)
		axEinterne.set_ylim(0.9*min_initial_e , 1.75*max_initial_e)
		axEinterne.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
		axEinterne.legend()

		plt.suptitle('Euler equations : exact solution of Riemann problem ' + filename)
		plt.savefig("EulerRiemannProblemExactSolution_" + filename + ".png")


	return 0.0;


if __name__ == """__main__""":
	run_Riemann_problems()
