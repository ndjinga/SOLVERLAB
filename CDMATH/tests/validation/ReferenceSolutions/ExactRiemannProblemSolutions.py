#!/usr/bin/env python3
# -*-coding:utf-8 -*

######################################################################################################################
#	This file runs a set of Riemann problems for the one dimensional Euler equations with stiffened gas equation of state. 
#	It uses the class exact_rs_stiffenedgas to compute the exact solution of the Riemann Problems.
#
#	Author: Michael Ndjinga
#	Date:	19/02/2021
#   Description : Translated from C++ package developped by Murray Cutforth
#######################################################################################################################

import exact_rs_stiffenedgas
from math import fabs

def stiffenedgas_e (rho, p, gamma, pinf):
	return (p+gamma*pinf)/(rho*(gamma-1));

def stiffenedgas_h (rho, p, gamma, pinf):
	return gamma*(p+pinf)/(rho*(gamma-1));


def run_Riemann_problems(numsamples = 100):
	# Output test solution for many different Riemann problems
	
	print("")
	print("Determination of the exact solutions of some Riemann problems for the Euler equations on " + str(numsamples) + " sample points.")

	for TC in range(1, 23):
		WL = [0.]*3
		WR = [0.]*3
	
		gammaL = 1.4;
		gammaR = 1.4;
		pinf_L = 0.0;
		pinf_R = 0.0;
		offset = 0.5;#position of the initial discontinuity
		xmin = 0.0;
		xmax = 1.0;


		# TC1 to TC5 are the 7 shock tube problems from Toro
	
		if (TC == 1):
			WL[0] = 1.0;
			WL[1] = 0.75;
			WL[2] = 1.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 0.1;
			t = 0.2;
			offset = 0.3;
			filename = "TTC1.dat";
		elif (TC ==2):
			WL[0] = 1.0;
			WL[1] = -2.0;
			WL[2] = 0.4;
			WR[0] = 1.0;
			WR[1] = 2.0;
			WR[2] = 0.4;
			t = 0.15;
			filename = "TTC2.dat";
		elif (TC ==3):
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1000.0;
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 0.01;
			t = 0.012;
			filename = "TTC3.dat";
		elif (TC ==4):
			WL[0] = 5.99924;
			WL[1] = 19.5975;
			WL[2] = 460.894;
			WR[0] = 5.99242;
			WR[1] = -6.19633;
			WR[2] = 46.0950;
			t = 0.035;
			offset = 0.4;
			filename = "TTC4.dat";
		elif (TC ==5):
			WL[0] = 1.0;
			WL[1] = -19.59745;
			WL[2] = 1000.;
			WR[0] = 1.0;
			WR[1] = -19.59745;
			WR[2] = 0.01;
			t = 0.012;
			offset = 0.8;
			filename = "TTC5.dat";
		elif (TC == 21):
			WL[0] = 1.4;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 1.;
			WR[1] = 0.0;
			WR[2] = 1.;
			t = 2.;
			filename = "TTC6.dat";
		elif (TC == 22):
			WL[0] = 1.4;
			WL[1] = 0.1;
			WL[2] = 1.0;
			WR[0] = 1.;
			WR[1] = 0.1;
			WR[2] = 1.;
			t = 2.;
			filename = "TTC7.dat";
		elif (TC == 6):
			# Air - Helium shock tube from Sambasivan 2009
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 0.1;
			t = 0.25;
			filename = "Samb1.dat";
			gammaR = 1.667;
		elif (TC == 7):
			# Gaseous shock tube from rGFM
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 100000.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 10000.0;
			t = 0.0007;
			filename = "NE1.dat";
			gammaR = 1.2;
		elif (TC == 8):
			# Air - water shock from rGFM
			
			WL[0] = 0.00596521;
			WL[1] = 911.8821;
			WL[2] = 1000.0;
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 1.0;
			gammaR = 7.15;
			pinf_R = 3309.0;
			t = 0.0007;
			filename = "rGFM2.dat";
		elif (TC == 9):
			# Air - water jet from rGFM
			
			WL[0] = 1.0;
			WL[1] = 90.0;
			WL[2] = 1.0;
			WR[0] = 1000.0;
			WR[1] = 0.0;
			WR[2] = 1.0;
			gammaR = 7.15;
			pinf_R = 3309.0;
			t = 0.015;
			filename = "rGFM4.dat";
			offset = 0.6;
		elif (TC == 10):
			# Reversed version of TC 9
			
			WR[0] = 1.0;
			WR[1] = -90.0;
			WR[2] = 1.0;
			WL[0] = 1000.0;
			WL[1] = 0.0;
			WL[2] = 1.0;
			gammaL = 7.15;
			pinf_L = 3309.0;
			t = 0.015;
			filename = "rGFM_reversed.dat";
			offset = 0.4;
		elif (TC == 11):
			# Water - air shock from Saurel 1999
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 4.0;
			WR[0] = 0.05;
			WR[1] = 0.0;
			WR[2] = 0.0004;
			gammaL = 4.4;
			pinf_L = 2.4;
			t = 0.12;
			filename = "Saurel1.dat";
			offset = 0.7;
		elif (TC == 12):
			# Reversed version of TC 11
			
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 4.0;
			WL[0] = 0.05;
			WL[1] = 0.0;
			WL[2] = 0.0004;
			gammaR = 4.4;
			pinf_R = 2.4;
			t = 0.12;
			filename = "Saurel1_reversed.dat";
			offset = 0.3;
		elif (TC == 13):
			# Numerical experiment 2
			
			WL[0] = 3.175962;
			WL[1] = 9.434992;
			WL[2] = 100.0;
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 1.0;
			gammaL = 1.667;
			gammaR = 1.2;
			t = 0.045;
			filename = "NE2.dat";
			offset = 0.2;
		elif (TC == 14):
			# Numerical experiment 3
			
			WL[0] = 0.00596521;
			WL[1] = 911.8821;
			WL[2] = 100.0;
			WR[0] = 1.0;
			WR[1] = 0.0;
			WR[2] = 1.0;
			gammaL = 1.4;
			gammaR = 7.15;
			pinf_L = 0.0;
			pinf_R = 3309.0;
			t = 0.0007;
			filename = "NE3.dat";
			offset = 0.5;
		elif (TC == 15):
			# Gaseous shock tube from So/Hu/Adams 2012
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 0.1;
			t = 0.15;
			filename = "ST3.dat";
			gammaR = 1.667;
		elif (TC == 16):
			# Gaseous shock tube from Garrick/Owkes/Regele 2016
			
			WL[0] = 1.0;
			WL[1] = 0.0;
			WL[2] = 1.0;
			WR[0] = 0.125;
			WR[1] = 0.0;
			WR[2] = 0.1;
			t = 0.14;
			filename = "ST1.dat";
			gammaR = 2.4;
		elif (TC == 17):
			# Water - air shock tube from Murrone/Guillard 2004
			
			WL[0] = 1000.0;
			WL[1] = 0.0;
			WL[2] = 1000000000.0;
			WR[0] = 50.0;
			WR[1] = 0.0;
			WR[2] = 100000.0;
			t = 0.0009;
			filename = "ST2.dat";
			offset = 0.7;
			gammaL = 4.4;
			gammaR = 1.4;
			pinf_L = 600000000.0;
			pinf_R = 0.0;
			xmin = -2.0;
			xmax = 2.0;
		elif (TC == 18):
			# Water shock tube from PWR - version 1
			
			WL[0] = 700.0;
			WL[1] = 0.0;
			WL[2] = 15500000.0;
			WR[0] = 700.0;
			WR[1] = 0.0;
			WR[2] = 100000.0;
			t = 3.e-4;
			filename = "PWR-ShockTube1.dat";
			gammaL = 1.58;
			gammaR = 1.58;
			pinf_L = 353637173.0;
			pinf_R = 353637173.0;
		elif (TC == 19):
			# Water shock tube from PWR - version 2
			
			WL[0] = 700.0;
			WL[1] = 0.0;
			WL[2] = 15500000.0;
			WR[0] = 700.0;
			WR[1] = 20.0;
			WR[2] = 15500000.0;
			t = 3.e-4;
			filename = "PWR-ShockTube2.dat";
			gammaL = 1.58;
			gammaR = 1.58;
			pinf_L = 353637173.0;
			pinf_R = 353637173.0;
		elif (TC == 20):
			# Water shock tube from PWR - version 3
			
			WL[0] = 700.0;
			WL[1] = 0.0;
			WL[2] = 15500000.0;
			WR[0] = 650.0;
			WR[1] = 0.0;
			WR[2] = 15500000.0;
			t = 3.e-4;
			filename = "PWR-ShockTube3.dat";
			gammaL = 1.58;
			gammaR = 1.58;
			pinf_L = 353637173.0;
			pinf_R = 353637173.0;
		else:
			print( "Unknown test case" )
			return 1;
	
	
		RS = exact_rs_stiffenedgas.exact_rs_stiffenedgas(gammaL, gammaR, pinf_L, pinf_R);
		RS.solve_RP(WL,WR);
	
		print( "")
		print( "Solved Riemann problem for TC = " , TC )
		print( "Star state pressure calculated as " , RS.P_STAR )
		print( "Star state velocity calculated as " , RS.S_STAR )
		print( "Left star state density calculated as " , RS.rho_star_L )
		print( "Right state state density calculated as " , RS.rho_star_R )
		print( "Left shock speed calculated as " , RS.S_L )
		print( "Right shock speed calculated as " , RS.S_R )
		print( "Left rarefaction head speed calculated as " , RS.S_HL )
		print( "Left rarefaction tail speed calculated as " , RS.S_TL )
		print( "Right rarefaction head speed calculated as " , RS.S_HR )
		print( "Right rarefaction tail speed calculated as " , RS.S_TR  )
	
		
		delx = (xmax - xmin)/numsamples;

		outfile=open(filename, 'w')
	
		x = xmin;
		while (x <= xmax):
		
			S = x/t;
			soln = RS.sample_solution(WL, WR, S - offset/t);
			thisgamma = gammaL if S - offset/t < RS.S_STAR else gammaR;
			thispinf = pinf_L  if S - offset/t < RS.S_STAR else pinf_R;
			thisz = 1.0  if S - offset/t < RS.S_STAR else 0.0;
			outfile.write( str( x ) + " " + str( soln[0]) + " " + str( soln[1]) + " " + str( soln[2]) + " " + str(stiffenedgas_e(soln[0], soln[2], thisgamma, thispinf)) + " " + str(stiffenedgas_h(soln[0], soln[2], thisgamma, thispinf)) + " " + str(fabs(soln[1])/RS.a(soln[0], soln[2], thisgamma, thispinf)) + " " + str(thisz) + "\n")
			x += delx;
		outfile.close()

	return 0.0;


if __name__ == """__main__""":
	run_Riemann_problems()
