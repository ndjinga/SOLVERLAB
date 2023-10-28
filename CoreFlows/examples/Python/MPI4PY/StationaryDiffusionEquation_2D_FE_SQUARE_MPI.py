#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Test launching one simulation on several procs (stationary diffusion quation)
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2023
# Description : Utilisation de la méthode des éléments finis avec champs T et f discrétisés aux noeuds d'un maillage triangulaire structurés
#               -\lambda\Delta T=\Phi 
#               Comparaison de la solution numérique avec la solution exacte T=sin(pi*x)*sin(pi*y)
#================================================================================================================================

from mpi4py import MPI
import numpy as np
import solverlab
from math import sin, pi

# Get communicator
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print("My rank is ", rank, " among ", size, "processors ")
# Prepare for the mesh
print("Processor ", rank, " : Building mesh " );
xinf = 0 ;
xsup=1.0;
yinf=0.0;
ysup=1.0;
nx=20;
ny=20; 
M=solverlab.Mesh(xinf,xsup,nx,yinf,ysup,ny,0)#Regular triangular mesh
# set the limit field for each boundary
eps=1e-6;
M.setGroupAtPlan(xsup,0,eps,"Bord1")
M.setGroupAtPlan(xinf,0,eps,"Bord2")
M.setGroupAtPlan(ysup,1,eps,"Bord3")
M.setGroupAtPlan(yinf,1,eps,"Bord4")

print("Processor ", rank, " : Built a regular triangular 2D mesh from a square mesh with ", nx,"x" ,ny, " cells.")

FEComputation=True
Lambda=1.#Thermal conductivity
spaceDim = 2

myProblem = solverlab.StationaryDiffusionEquation(spaceDim,FEComputation, Lambda, comm);
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

# name of result file
fileName = "StationnaryDiffusion_2DEF_StructuredTriangles";

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
	print("Processor ", rank, " : simulation succeeded ")
	if( rank == 0) :
		my_ResultField = myProblem.getOutputTemperatureField()
		max_sol_num=my_ResultField.max()
		min_sol_num=my_ResultField.min()
		erreur_abs=0
		for i in range(M.getNumberOfCells()) :
			if erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) :
				erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i])
		
		print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
	pass

print("Processor ", rank, " : ------------ !!! End of calculation !!! -----------" );

myProblem.terminate();

