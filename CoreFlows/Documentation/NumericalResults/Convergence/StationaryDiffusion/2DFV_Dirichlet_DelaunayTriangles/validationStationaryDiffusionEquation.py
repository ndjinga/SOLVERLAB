#!/usr/bin/env python
# -*-coding:utf-8 -*
#=====================================================================================================================================
# Name        : Résolution EF/VF de l'équation de Poisson 2D/3D -\triangle T = f avec conditions aux limites de Dirichlet ou de Neumann homogènes
# Author      : Michaël Ndjinga, Sédrick Kameni Ngwamou
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des éléménts finis P1 ou des volumes finis avec champs T et f discrétisés aux noeuds ou aux cellules 
#				 d'un maillage triangulaire / carré / cubique /tétrahèdrique.
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte T=-sin(pi*x)*sin(pi*y) (Dirichlet) ou T=-sin(pi*x)*sin(pi*y)*sin(pi*z) (Neumann
#=====================================================================================================================================

import CoreFlows as cf
import cdmath as cm
from math import sin, pi, cos
import time,json
import PV_routines
import VTK_routines

test_desc={}
test_desc["Initial_data"]="None"
test_desc["PDE_model"]="Poisson"
test_desc["PDE_is_stationary"]=True
test_desc["PDE_search_for_stationary_solution"]=False
test_desc["Mesh_is_unstructured"]=True
test_desc["Part_of_mesh_convergence_analysis"]=True
test_desc["Numerical_method_time_discretization"]="None"


def SolveStationaryDiffusionEquation(my_mesh,resolution,MeshType,method,BC):
	start = time.time()
	# test_desc["Mesh_type"]=meshType ****filename,

	spaceDim=my_mesh.getSpaceDimension()
	
	# Prepare for the mesh groups
	print("Building mesh " );
	xinf = 0 ;
	xsup = 1.0;
	yinf = 0.0;
	ysup = 1.0;
	if spaceDim == 3:
		zinf = 0.0;
		zsup = 1.0;
	
	#set the boundary groups
	eps=1e-6;
	my_mesh.setGroupAtPlan(xsup,0,eps,"Right")
	my_mesh.setGroupAtPlan(xinf,0,eps,"Left")
	my_mesh.setGroupAtPlan(ysup,1,eps,"Top")
	my_mesh.setGroupAtPlan(yinf,1,eps,"Bottom")
	if spaceDim == 3:
		my_mesh.setGroupAtPlan(zsup,2,eps,"Front")
		my_mesh.setGroupAtPlan(zinf,2,eps,"Back")
	
	nbCells = my_mesh.getNumberOfCells()	

	test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
	test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
	test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfCells()
	test_desc["Mesh_cell_type"]=my_mesh.getElementTypesNames()
		
	#Define the right hand side function
	if method == 'FE':
		test_desc["Numerical_method_name"]="FE"
		test_desc["Numerical_method_space_discretization"]="Finite elements"
		FEComputation=True
		my_RHSfield = cm.Field("RHS_field", cm.NODES, my_mesh, 1)
		for i in range(my_mesh.getNumberOfNodes()):
			Ni= my_mesh.getNode(i)
			x = Ni.x()
			y = Ni.y()
			if BC == 'Dirichlet' :	
				if spaceDim == 2 : 
					my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
				elif spaceDim == 3 : 
					z = Ni.z()
					my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)#mettre la fonction definie au second membre de l'edp
			elif BC =='Neumann':
				if spaceDim == 2 : 
					my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y)#set the function define in the right hand side 				
				elif spaceDim == 3 : 
					z = Ni.z()
					my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y)*cos(pi*z)#set the function define in the right hand side 				
	elif method == 'FV':
		test_desc["Numerical_method_name"]="FV5"
		test_desc["Numerical_method_space_discretization"]="Finite volumes"
		test_desc["Global_comment"]="2 points FV diffusion scheme"
		FEComputation=False
		my_RHSfield = cm.Field("RHS_field", cm.CELLS, my_mesh, 1)
		for i in range(my_mesh.getNumberOfCells()):
			Ci= my_mesh.getCell(i)
			x = Ci.x()
			y = Ci.y()
			if BC == 'Dirichlet' :	
				if spaceDim == 2 : 
					my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l'edp
				elif spaceDim == 3 : 
					z = Ci.z()
					my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)#mettre la fonction definie au second membre de l'edp
			elif BC =='Neumann':
				if spaceDim == 2 : 
					my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y)#set the function define in the right hand side 
				elif spaceDim == 3 : 
					z = Ci.z()
					my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y)*cos(pi*z)#set the function define in the right hand side 
				
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(my_mesh);
	myProblem.setHeatPowerField(my_RHSfield)

	#Set the boundary conditions
	if BC == 'Dirichlet' :	
		test_desc["Boundary_conditions"]="Dirichlet"
		T1=0;
		T2=0;
		T3=0;
		T4=0;
		myProblem.setDirichletBoundaryCondition("Right",T1)
		myProblem.setDirichletBoundaryCondition("Left",T2)
		myProblem.setDirichletBoundaryCondition("Top",T3)
		myProblem.setDirichletBoundaryCondition("Bottom",T4)
		if spaceDim == 3:
			T5=0;
			T6=0;
			myProblem.setDirichletBoundaryCondition("Front",T5)
			myProblem.setDirichletBoundaryCondition("Back",T6)
	elif BC =='Neumann':
		test_desc["Boundary_conditions"]="Neumann"
		myProblem.setNeumannBoundaryCondition("Right")
		myProblem.setNeumannBoundaryCondition("Left")
		myProblem.setNeumannBoundaryCondition("Top")
		myProblem.setNeumannBoundaryCondition("Bottom")
		if spaceDim == 3:
			myProblem.setNeumannBoundaryCondition("Front")
			myProblem.setNeumannBoundaryCondition("Back")
		if method == 'FE' :
			myProblem.setLinearSolver(cf.GMRES,cf.ILU);#LU solvers breaks down

	if spaceDim == 2 : 
		test_desc["Geometry"]="Square"
		if method == 'FE':
			test_desc["Global_name"]="FE simulation of the 2D Poisson equation"
		elif method == 'FV':
			test_desc["Global_name"]="FV simulation of the 2D Poisson equation"
	elif spaceDim == 3:
		test_desc["Geometry"]="Cube"
		if method == 'FE':
			test_desc["Global_name"]="FE simulation of the 3D Poisson equation"
		elif method == 'FV':
			test_desc["Global_name"]="FV simulation of the 3D Poisson equation"

	# name of result file
	fileName = str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells);

	# computation parameters
	myProblem.setFileName(fileName);

	# Run the computation
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.solveStationaryProblem();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + " is successful !" );
		####################### Postprocessing #########################
		my_ResultField = myProblem.getOutputTemperatureField()
		#The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi
		max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(spaceDim*pi*pi)
		max_sol_num=my_ResultField.max()
		min_sol_num=my_ResultField.min()
		erreur_abs=0
		if method =='FE':
			for i in range(my_mesh.getNumberOfNodes()) :
				if erreur_abs < abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) :
					erreur_abs = abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i])
		else:
			for i in range(my_mesh.getNumberOfCells()) :
				if erreur_abs < abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) :
					erreur_abs = abs(my_RHSfield[i]/(spaceDim*pi*pi) - my_ResultField[i]) 				
		print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
		print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
		print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
			
		assert erreur_abs/max_abs_sol_exacte <1.
		pass
	end = time.time()
	test_desc["Absolute_error"]=erreur_abs
	test_desc["Relative_error"]=erreur_abs/max_abs_sol_exacte
	test_desc["Computational_time_taken_by_run"]=end-start

	print( "------------ !!! End of calculation !!! -----------" );

	myProblem.terminate();
	
	#Postprocessing :
    #================
	if spaceDim == 2 :
		diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,0,0],[1,1,0],resolution)
		if method =='FE':
			PV_routines.Save_PV_data_to_picture_file("StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',"Temperature",'NODES',"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+str(MeshType)+str(BC)+str(nbCells))
		else:
			PV_routines.Save_PV_data_to_picture_file("StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',"Temperature",'CELLS',"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+str(MeshType)+str(BC)+str(nbCells))
	else:
		diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,0,0],[1,1,1],resolution)
		VTK_routines.Clip_VTK_data_to_VTK("StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',[0.5,0.5,0.5], [-0.5,-0.5,-0.5],resolution )
		if method =='FE':
			PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',"Temperature",'NODES',"Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+str(MeshType)+str(BC)+str(nbCells))
		else:
			PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+'_'+str(MeshType)+'_'+str(BC)+str(nbCells)+'_0.vtu',"Temperature",'CELLS',"Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_"+str(spaceDim)+'D'+str(method)+str(MeshType)+str(BC)+str(nbCells))

	with open('test_Poisson'+str(my_mesh.getMeshDimension())+'D_'+str(method)+'_'+str(nbCells)+'Cells'+ '_' +str(BC)+'_'+str(MeshType)+".json", 'w') as outfile:  
		json.dump(test_desc, outfile)
	
	return erreur_abs/max_abs_sol_exacte, nbCells, diag_data, min_sol_num, max_sol_num, end - start 

if __name__ == """__main__""":
	mesh51 = cm.Mesh(0,1,11,0,1,11,0,1,11,1)
	Solver(mesh51,100,"Regular_tetrahedron","Green","FE","Dirichlet",3)
