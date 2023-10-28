#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Test launching one simulation on several procs (diffusion quation)
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2023
# Description : Test solving the diffusion of the temperature T in a solid subject to a heat source
#               \rho cp dT/dt-\lambda\Delta T=\Phi 
#               The heat source function is an eigenvector of the Laplacean, so the exact solution is known using separation of variables
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
spaceDim = 2

initialTemperature=0

# Mandatory physical values
solid_specific_heat=300# specific heat capacity, default value 300
solid_density=10000# density, default value 10000
solid_conductivity=5# conductivity, default value 5
diffusivity=solid_conductivity/(solid_density*solid_specific_heat)

myProblem = solverlab.DiffusionEquation(spaceDim,FEComputation,solid_density,solid_specific_heat,solid_conductivity, comm);
myProblem.setInitialFieldConstant(M,[initialTemperature],solverlab.NODES);

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
fileName = "Diffusion_2DEF_StructuredTriangles";

# set the numerical method
myProblem.setTimeScheme( solverlab.Implicit)
max_nb_its_lin_solver = 50
myProblem.setLinearSolver(solverlab.GMRES, solverlab.ILU, max_nb_its_lin_solver );

# computation parameters
myProblem.setFileName(fileName);

MaxNbOfTimeStep = 3 ;# default value is 10
freqSave = 1;# default value is 1
cfl = 100;# default value is 1
maxTime = 100000000;# default value is 10
precision = 1e-6;# default value is 1e-6
result_directory="."# default value = current directory

myProblem.setCFL(cfl);
myProblem.setPrecision(precision);
myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
myProblem.setTimeMax(maxTime);
myProblem.setFreqSave(freqSave);
myProblem.setFileName(fileName);
myProblem.setResultDirectory(result_directory)
#myProblem.setSaveFileFormat(solverlab.MED)#default value is solverlab.VTK

# Run the computation
myProblem.initialize();
print("Processor ", rank, " : Running python "+ fileName );

ok = myProblem.run();
if (not ok):
	print( "Python simulation of " + fileName + "  failed ! " );
	pass
else:
	print("Processor ", rank, " : simulation succeeded ")
	pass

print("Processor ", rank, " : ------------ !!! End of calculation !!! -----------" );

myProblem.terminate();

