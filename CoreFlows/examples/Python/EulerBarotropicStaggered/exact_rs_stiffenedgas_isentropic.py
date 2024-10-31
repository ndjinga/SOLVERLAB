#!/usr/bin/env python3
# -*-coding:utf-8 -*

######################################################################################################################
#	This file contains a class to solve for the exact solution of the Riemann Problem for the one dimensional 
#	isentropic Euler equations with Laplace law equation of state P=c1*rho^gamma + pinf
#
#	Author: Michael Ndjinga
#	Date:	15/09/2022
#   Description : Translated from C++ package developped by Murray Cutforth
#######################################################################################################################

from math import pow, fabs, sqrt, log

# Functions used to generate exact solutions to Riemann problems
def exact_sol_Riemann_problem(xmin, xmax, t, gamma, c1, WL, WR, offset, numsamples = 100):#offset= position of the initial discontinuity 
	print("")
	print("Determination of the exact solution of the Riemann problem for the isentropic Euler equations, p = c1 rho^\gamma, gamma = ", gamma, ", c1 = ", c1, ", t = ", t)

	RS = exact_rs_stiffenedgas_isentropic(); #(gamma, c1, pinf, tol, max_iter);
	RS.solve_RP(WL,WR);

	delx = (xmax - xmin)/numsamples;
	
	velocity  = [0.]*numsamples
	density = [0.]*numsamples

	for i in range(numsamples):
		S = (xmin +i*delx)/t;
		soln = RS.sample_solution(WL, WR, S - offset/t);
		density[i] = RS.rho( soln[0] )
		velocity[i]= soln[1]

	return density, velocity 


class exact_rs_stiffenedgas_isentropic :

	def __init__(self, gamma=2., c1=1.0, tol=1.e-9, max_iter=200):
		self.TOL = tol
		self.MAX_NB_ITER = max_iter
	
		assert gamma > 1, "gamma should be greater than 1" # TO DO : case where gamma ==1
		assert c1    >= 0, "c1 should be positive"
		self.gamma = gamma
		self.c1    = c1
		
		self.S_STAR = 0.
		self.P_STAR = 0.
		
		self.S_L = 0.
		self.S_R = 0.
		self.S_HL = 0.
		self.S_TL = 0.
		self.S_HR = 0.
		self.S_TR = 0.

	
	def rho(self, p):
		return pow(p /self.c1, 1/self.gamma)
	
	
	def solve_RP (self, W_L, W_R):
		assert len(W_L) == 2, "Left state should have two components (p, u)"
		assert len(W_R) == 2, "Right state should have two components (p, u)"
		assert W_L[0] >= 0.0, "Left pressure should be positive"
		assert W_R[0] >= 0.0, "Right pressure should be positive"
		
		print("")
		print("Solving Riemann problem for left state W_L=", W_L, ", and right state W_R=",W_R)
		
		# Calculate p_star, u_star
		self.P_STAR = self.find_p_star_newtonraphson(W_L[0], W_L[1], W_R[0], W_R[1])		
		self.S_STAR = 0.5*(W_L[1]+W_R[1]) + 0.5*(self.f_R(self.P_STAR,W_R[0]) - self.f_L(self.P_STAR,W_L[0]))
	
		# Solution now depends on character of 1st and 2nd waves
		if (self.P_STAR > W_L[0]):
			# Left shock speed
			self.S_L =  (self.rho(W_L[0])* W_L[1] - self.rho(self.P_STAR)*self.S_STAR )/(self.rho(W_L[0]) - self.rho(self.P_STAR))   #To do W_L[0] - self.Q_K( self.P_STAR, W_L[0])/self.rho(W_L[0]) 
		else:
			# Left rarefaction : determining the velocity bounding the rarefaction fan
			self.S_HL = W_L[1] - self.a(W_L[0])
			self.S_TL = self.S_STAR - self.a(self.P_STAR)
	
		if (self.P_STAR > W_R[0]):
			# Right shock speed
			self.S_R = (self.rho(self.P_STAR)*self.S_STAR - self.rho(W_R[0])* W_R[1]  )/(self.rho(self.P_STAR) - self.rho(W_R[0]) ) #To do  W_R[0] + (self.Q_K(self.P_STAR,W_R[0])/self.rho(W_R[0]))
		else:
			# Right rarefaction : determining the velocity bounding the rarefaction fan
			self.S_HR = W_R[1] + self.a(W_R[0])
			self.S_TR = self.S_STAR + self.a(self.P_STAR)
		
		print("self.S_L =",self.S_L, "self.S_R = ", self.S_R, "self.S_STAR = ", self.S_STAR)
		print("P_Star = ", self.P_STAR)

	
	def sample_solution (self, W_L, W_R, S):
		W = [0.]*2
		# Find appropriate part of solution and return (rho, u)
		if (S < self.S_STAR):
			# To the left of the contact
			if (self.P_STAR > W_L[0]): # Left shock
				if (S < self.S_L):
					W = W_L
				else:
					W[1] = self.S_STAR
					W[0] = self.P_STAR
			else: # Left rarefaction
				if (S < self.S_HL):
					W = W_L
				else:
					if (S > self.S_TL):
						W[1] = self.S_STAR
						W[0] = self.P_STAR
					else:
						W = self.set_left_rarefaction_fan_state(W_L, S)
		else:
			# To the right of the contact
			if (self.P_STAR > W_R[0]): # Right shocks
				if (S > self.S_R):
					W = W_R
				else:
					W[1] = self.S_STAR
					W[0] = self.P_STAR
			else: # Right rarefaction
				if (S > self.S_HR):
					W = W_R
				else:
					if (S < self.S_TR):
						W[1] = self.S_STAR
						W[0] = self.P_STAR
					else:
						W = self.set_right_rarefaction_fan_state(W_R, S) 
		return W
	
	# Functions used to solve for p_star iteratively
	def find_p_star_newtonraphson (self, p_L, u_L, p_R, u_R ):
		# First we set the initial guess for p_star using a simple mean-value approximation
		p_star_next = 0.5*(p_L+p_R)
		n = 0
		# Now use the Newton-Raphson algorithm
		while True:#conversion of do ... while by while True... if (...) break
			p_star = p_star_next
			p_star_next = p_star - self.total_pressure_function(p_star,p_L,u_L,p_R,u_R)/self.total_pressure_function_deriv(p_star,p_L,p_R)
			p_star_next = max(p_star_next, self.TOL)
			n+=1
			if not ((fabs(p_star_next - p_star)/(0.5*(p_star+p_star_next)) > self.TOL) and n < self.MAX_NB_ITER):
				break
		if (n == self.MAX_NB_ITER):
			raise ValueError("!!!!!!!!!!Newton algorithm did not converge. Increase tolerance or maximum number of time steps. Current values : tol=" + str(self.TOL) + ", max_iter=" + str(self.MAX_NB_ITER) )
			#p_star_next = 0.5*(p_L+p_R)
		return p_star_next


	def total_pressure_function (self, p_star, p_L, u_L, p_R, u_R):
		return	self.f_L(p_star, p_L)	+ self.f_R(p_star, p_R) + u_R - u_L

	def total_pressure_function_deriv (self, p_star, p_L, p_R ): 
		return 	self.f_L_deriv (p_star, p_L) + self.f_R_deriv (p_star, p_R)



	def f_L (self, p_star, p):
		if (p_star > p):
			return self.RankineHugoniot(p_star, p) 
		else:
			if (self.gamma>1):
				return (2.0*self.a(p)/(self.gamma-1.0))*(pow((p_star)/p,  (self.gamma-1.0)/(2.0*self.gamma)) -1.0 ) 
			else:
				return sqrt(self.c1)*log( p_star/p  )
		
	def f_L_deriv (self, p_star, p):
		rho = pow(p_star/self.c1, 1/self.gamma)
		if (p_star > p):
			""" A = 2.0/((self.gamma+1.0)*rho)
			B = (p)*(self.gamma-1.0)/(self.gamma+1.0)
			return sqrt(A/(B+p_star))*(1.0 - ((p_star-p)/(2.0*(B+p_star)))) # To Do : à vérif """
			return (  -1/(self.gamma * pow(self.c1,1/self.gamma))*(p - p_star) * pow(p_star, -(self.gamma+1))   + (1.0/self.rho(p)- 1.0/self.rho(p_star)) )/( 2.0*  self.f_L(p_star, p) )
		else:
			if (self.gamma>1):
				return (1.0/(rho*self.a(p)))*pow(p_star/p, -(self.gamma+1.0)/(2.0*self.gamma) ) 
			else:
				return sqrt(self.c1)/(p_star)

	def f_R (self, p_star, p):
		if (p_star > p):
			return self.RankineHugoniot(p_star, p) 
		else:
			if (self.gamma>1):
				return (2.0*self.a(p)/(self.gamma-1.0))*(pow((p_star)/p,  (self.gamma-1.0)/(2.0*self.gamma)) -1.0 ) 
			else:
				return sqrt(self.c1)*log( p_star/p  )
		
	def f_R_deriv (self, p_star, p):
		rho = pow(p_star/self.c1, 1/self.gamma)
		if (p_star > p):
			""" A = 2.0/((self.gamma+1.0)*rho)
			B = (p)*(self.gamma-1.0)/(self.gamma+1.0)
			return sqrt(A/(B+p_star))*(1.0 - ((p_star-p)/(2.0*(B+p_star)))) # To Do : à vérif """
			return (  -1/(self.gamma * pow(self.c1,1/self.gamma))*(p - p_star) * pow(p_star, -(self.gamma+1))   + (1.0/self.rho(p)- 1.0/self.rho(p_star)) )/( 2.0*  self.f_R(p_star, p) )
		else:
			if (self.gamma>1):
				return (1.0/(rho*self.a(p)))*pow(p_star/p, -(self.gamma+1.0)/(2.0*self.gamma) ) 
			else:
				return sqrt(self.c1)/(p_star)
		

	# Functions to find the state inside a rarefaction fan

	def set_left_rarefaction_fan_state (self, W_L, S): # To Do : à vérif
		W = [0.]*2
		a_L = self.a(W_L[0])
		W[1] = (2.0/(self.gamma+1.0))*(a_L + S + ((self.gamma-1.0)/2.0)*W_L[1])
		W[0] = (W_L[0] )*pow((2.0/(self.gamma+1.0)) + ((self.gamma-1.0)/(a_L*(self.gamma+1.0)))*(W_L[1] - S), (2.0*self.gamma)/(self.gamma-1.0))  
		return W

	def set_right_rarefaction_fan_state (self, W_R, S): 
		W = [0.]*2
		a_R = self.a(W_R[0])
		W[1] = (2.0/(self.gamma+1.0))*(- a_R + S + ((self.gamma-1.0)/2.0)*W_R[1])
		W[0] = (W_R[0] )*pow((2.0/(self.gamma+1.0)) - ((self.gamma-1.0)/(a_R*(self.gamma+1.0)))*(W_R[1] - S), (2.0*self.gamma)/(self.gamma-1.0)) 
		return W


	# Misc functions
	def RankineHugoniot(self, p_star, p):
		return pow(  -(p_star - p) * (1/self.rho(p_star) - 1/self.rho(p) ) , 1/2.0 ) 

	def Q_K (self, p_star, p): # To Do : à vérif
		return (p_star - p)/(1/self.rho(p) - 1/self.rho(p_star)) 

	# Equation of state functions
	def a (self, p):#sound speed
		return sqrt(self.gamma*( p/ self.rho(p) )) 

















	#Determine the solution value at position x and time t
	def p_u_solution (initialLeftState, initialRightState, x, t, gamma, c1=1, pinf=0, offset=0, tol=1.e-6, max_iter=100):
		RS = exact_rs_stiffenedgas_isentropic(gamma, c1, pinf, tol, max_iter);
		RS.solve_RP(initialLeftState, initialRightState);
		return 	RS.sample_solution(initialLeftState, initialRightState, (x - offset)/t);


