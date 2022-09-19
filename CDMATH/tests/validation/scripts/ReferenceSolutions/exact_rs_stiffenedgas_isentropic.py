#!/usr/bin/env python3
# -*-coding:utf-8 -*

######################################################################################################################
#	This file contains a class to solve for the exact solution of the Riemann Problem for the one dimensional isentropic Euler
#	equations with Laplace law equation of state P=c1*rho^gamma + pinf
#
#	Author: Michael Ndjinga
#	Date:	15/09/2022
#   Description : Translated from C++ package developped by Murray Cutforth
#######################################################################################################################

from math import pow, fabs, sqrt

class exact_rs_stiffenedgas_isentropic :

	def __init__(self, gamma_L, gamma_R, tol=1.e-6, max_iter=100):
		self.TOL = tol
		self.MAX_NB_ITER = max_iter
	
		self.gamma_L = gamma_L
		self.gamma_R = gamma_R
		
		self.S_STAR = 0.
		self.P_STAR = 0.
		
		self.S_L = 0.
		self.S_R = 0.
		self.S_HL = 0.
		self.S_TL = 0.
		self.S_HR = 0.
		self.S_TR = 0.

	
	
	# Functions used to generate exact solutions to Riemann problems

	def solve_RP (self, W_L, W_R):
		assert len(W_L) == 2, "Left state should have two components (p, u)"
		assert len(W_R) == 2, "Right state should have two components (p, u)"
		assert W_L[0] >= 0.0, "Left density should be positive"
		assert W_R[0] >= 0.0, "Right density should be positive"
		
		print("")
		print("Solving Riemann problem for left state W_L=", W_L, ", and right state W_R=",W_R)
		
		# Calculate p_star
	
		self.P_STAR = self.find_p_star_newtonraphson(W_L[0], W_L[1], W_R[0], W_R[1])	
			
		
		# Calculate u_star
	
		self.S_STAR = 0.5*(W_L[1]+W_R[1]) + 0.5*(self.f(self.P_STAR,W_R[0],self.gamma_R) - self.f(self.P_STAR,W_L[0],self.gamma_L))
	
	
	def sample_solution (self, W_L, W_R, S):
		W = [0.]*2
		
		# Find appropriate part of solution and return primitives
	
		if (S < self.S_STAR):
			# To the left of the contact
	
			if (self.P_STAR > W_L[0]):
				# Left shock
				
				if (S < self.S_L):
					W = W_L
				else:
					W[1] = self.S_STAR
					W[0] = self.P_STAR
			else:
				# Left rarefaction
				
				if (S < self.S_HL):
					W = W_L
				else:
					if (S > self.S_TL):
						W[1] = self.S_STAR
						W[0] = self.P_STAR
					else:
						self.set_left_rarefaction_fan_state(W_L, S, W)
		else:
			# To the right of the contact
	
			if (self.P_STAR > W_R[0]):
				# Right shock
				
				if (S > self.S_R):
					W = W_R
				else:
					W[1] = self.S_STAR
					W[0] = self.P_STAR
			else:
				# Right rarefaction
				
				if (S > self.S_HR):
					W = W_R
				else:
					if (S < self.S_TR):
						W[1] = self.S_STAR
						W[0] = self.P_STAR
					else:
						self.set_right_rarefaction_fan_state(W_R, S, W)
	
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

		return	self.f(p_star, p_L, self.gamma_L)	+ self.f(p_star, p_R, self.gamma_R) + u_R - u_L

	def total_pressure_function_deriv (self, p_star, p_L, p_R ):

		return 	self.f_deriv (p_star, p_L, self.gamma_L) + self.f_deriv (p_star, p_R, self.gamma_R)


	def f (self, p_star, p, gamma):
		if (p_star > p):
		
			return (p_star - p)/self.Q_K(p_star, p, gamma)
		
		else:
		
			return (2.0*self.a(p,gamma,pinf)/(gamma-1.0))*(pow((p_star )/(p ), (gamma-1.0)/(2.0*gamma)) - 1.0)
		

	def f_deriv (self, p_star, rho, p, gamma):
		A = 2.0/((gamma+1.0)*rho)
		B = (p+pinf)*(gamma-1.0)/(gamma+1.0)
	
		if (p_star > p):
		
			return sqrt(A/(B+p_star+pinf))*(1.0 - ((p_star-p)/(2.0*(B+p_star+pinf))))
		
		else:
		
			return (1.0/(rho*self.a(p,gamma,pinf)))*pow((p_star+pinf)/(p+pinf), -(gamma+1.0)/(2.0*gamma))
		


	# Functions to find the state inside a rarefaction fan

	def set_left_rarefaction_fan_state (self, W_L, S, W):
		a_L = self.a(W_L[0],self.gamma_L)
		W[1] = (2.0/(self.gamma_L+1.0))*(a_L + S + ((self.gamma_L-1.0)/2.0)*W_L[1])
		W[0] = (W_L[0] + self.pinf_L)*pow((2.0/(self.gamma_L+1.0)) + ((self.gamma_L-1.0)/(a_L*(self.gamma_L+1.0)))*(W_L[1] - S), (2.0*self.gamma_L)/(self.gamma_L-1.0)) 

	def set_right_rarefaction_fan_state (self, W_R, S, W):
		a_R = self.a(W_R[0],self.gamma_R)
		W[1] = (2.0/(self.gamma_R+1.0))*(- a_R + S + ((self.gamma_R-1.0)/2.0)*W_R[1])
		W[0] = (W_R[0] + self.pinf_R)*pow((2.0/(self.gamma_R+1.0)) - ((self.gamma_R-1.0)/(a_R*(self.gamma_R+1.0)))*(W_R[1] - S), (2.0*self.gamma_R)/(self.gamma_R-1.0)) 



	# Misc functions

	def Q_K (self, p_star, rho, p, gamma):
		A = 2.0/((gamma+1.0)*rho)
		B = (p+pinf)*(gamma-1.0)/(gamma+1.0)
		return sqrt((p_star+pinf+B)/A)

	

	# Equation of state functions

	def a (self, rho, p, gamma, pinf):#sound speed
		return sqrt(gamma*((p+pinf)/rho))

	#Determine the solution value at position x and time t
	def rho_u_p_solution (initialLeftState, initialRightState, x, t, gamma, pinf, offset=0):
		RS = exact_rs_stiffenedgas(gamma, gamma, pinf, pinf);
		RS.solve_RP(initialLeftState, initialRightState);
		return 	RS.sample_solution(initialLeftState, initialRightState, (x - offset)/t);


