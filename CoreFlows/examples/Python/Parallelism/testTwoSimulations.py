#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Tests of launching two independent simulations in parallel
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : 
#================================================================================================================================

from mpi4py import MPI
import numpy as np
import solverlab
from math import sin, pi

def StationaryDiffusionEquation_2DEF_StructuredTriangles_par(split_direction, rank):
	spaceDim = 2;
	# Prepare for the mesh
	print("Processor ", rank, " : Building mesh " );
	xinf = 0 ;
	xsup=1.0;
	yinf=0.0;
	ysup=1.0;
	nx=20;
	ny=20; 
	M=solverlab.Mesh(xinf,xsup,nx,yinf,ysup,ny,split_direction)#Regular triangular mesh
	# set the limit field for each boundary
	eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"Bord1")
	M.setGroupAtPlan(xinf,0,eps,"Bord2")
	M.setGroupAtPlan(ysup,1,eps,"Bord3")
	M.setGroupAtPlan(yinf,1,eps,"Bord4")
	
	print("Processor ", rank, " : Built a regular triangular 2D mesh from a square mesh with ", nx,"x" ,ny, " cells.")
	print("Processor ", rank, " : Each square was split in two in direction ",split_direction)
	FEComputation=True
	myProblem = solverlab.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(M);

    # set the limit value for each boundary
	T1=0;
	T2=0;
	T3=0;
	T4=0;
    
	myProblem.setDirichletBoundaryCondition("Bord1",T1)
	myProblem.setDirichletBoundaryCondition("Bord2",T2)
	myProblem.setDirichletBoundaryCondition("Bord3",T3)
	myProblem.setDirichletBoundaryCondition("Bord4",T4)

	#Set the right hand side function
	my_RHSfield = solverlab.Field("RHS_field", solverlab.NODES, M, 1)
	for i in range(M.getNumberOfNodes()):
		Ni= M.getNode(i)
		x = Ni.x()
		y = Ni.y()

		my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
	
	myProblem.setHeatPowerField(my_RHSfield)
	myProblem.setLinearSolver(solverlab.GMRES,solverlab.ILU);

    # name of result file
	fileName = "StationnaryDiffusion_2DEF_StructuredTriangles"+str(rank);

    # computation parameters
	myProblem.setFileName(fileName);

    # Run the computation
	myProblem.initialize();
	print("Processor ", rank, " : Running python "+ fileName );

	ok = myProblem.solveStationaryProblem();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		####################### Postprocessing #########################
		my_ResultField = myProblem.getOutputTemperatureField()
		#The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi
		max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(2*pi*pi)
		max_sol_num=my_ResultField.max()
		min_sol_num=my_ResultField.min()
		erreur_abs=0
		for i in range(M.getNumberOfNodes()) :
			if erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) :
				erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i])
		
		print("Processor ", rank, " : Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Processor ", rank, " : Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
		print("Processor ", rank, " : Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
		
		assert erreur_abs/max_abs_sol_exacte <1.
		pass

	print("Processor ", rank, " : ------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	return erreur_abs/max_abs_sol_exacte

if __name__ == """__main__""":
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()
	
	if(size!=2):
		raise ValueError("Processor ", rank, " : aborting.\n Simulation should done on two processors.\n", size, " processors given")
		
	print("My rank is ", rank, " among ", size, "processors ")
	
	my_relative_error=StationaryDiffusionEquation_2DEF_StructuredTriangles_par(rank, rank)

	if rank == 0:
	    comm.send(my_relative_error, dest=1, tag=11)
	    other_relative_error = comm.recv(source=1, tag=17)
	elif rank == 1:
	    other_relative_error = comm.recv(source=0, tag=11)
	    comm.send(my_relative_error, dest=0, tag=17)
	
	print("Processor ", rank, " : Difference between the two processor relative errors is ", abs(my_relative_error-other_relative_error) )
	#print("Processor ", rank, " : Difference between the two processors is ", (my_ResultField-other_ResultField).normMax()[0] )
	
