# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson -\triangle u = f avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage quelconque
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution num"rique avec la solution exacte u=sin(pi*x)*sin(pi*y)
#================================================================================================================================

import cdmath
import time, json
from math import sin, pi
import PV_routines
import VTK_routines

test_desc={}
test_desc["Initial_data"]="None"
test_desc["Boundary_conditions"]="Dirichlet"
test_desc["Global_name"]="FV simulation of the 2D Poisson equation"
test_desc["Global_comment"]="2 points FV diffusion scheme"
test_desc["PDE_model"]="Poisson"
test_desc["PDE_is_stationary"]=True
test_desc["PDE_search_for_stationary_solution"]=False
test_desc["Numerical_method_name"]="VF5"
test_desc["Numerical_method_space_discretization"]="Finite volumes"
test_desc["Numerical_method_time_discretization"]="None"
test_desc["Mesh_is_unstructured"]=True
test_desc["Geometry"]="Square"
test_desc["Part_of_mesh_convergence_analysis"]=True

def solve(my_mesh,filename,resolution, meshType, testColor):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    
    nbCells = my_mesh.getNumberOfCells()
    
    if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
        raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")

    test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
    test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
    test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfCells()
    test_desc["Mesh_cell_type"]=my_mesh.getElementTypesNames()

    print("Mesh groups done")
    print("Number of cells  = ", nbCells)
    
    #Discrétisation du second membre et extraction du nb max de voisins d'une cellule
    #================================================================================
    my_RHSfield = cdmath.Field("RHS_field", cdmath.CELLS, my_mesh, 1)
    maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
    #parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
    for i in range(nbCells): 
        Ci = my_mesh.getCell(i)
        x = Ci.x()
        y = Ci.y()

        my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l edp
        # compute maximum number of neighbours
        maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)
    
    test_desc["Mesh_max_number_of_neighbours"]=maxNbNeighbours

    print("Right hand side discretisation done")
    print("Max nb of neighbours=", maxNbNeighbours)
    
    # Construction de la matrice et du vecteur second membre du système linéaire
    #===========================================================================
    Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours) # warning : third argument is maximum number of non zero coefficients per line of the matrix
    RHS=cdmath.Vector(nbCells)
    #Parcours des cellules du domaine
    for i in range(nbCells):
        RHS[i]=my_RHSfield[i] #la valeur moyenne du second membre f dans la cellule i
        Ci=my_mesh.getCell(i)
        for j in range(Ci.getNumberOfFaces()):# parcours des faces voisinnes
            Fj=my_mesh.getFace(Ci.getFaceId(j))
            if not Fj.isBorder():
                k=Fj.getCellId(0)
                if k==i :
                    k=Fj.getCellId(1)
                Ck=my_mesh.getCell(k)
                distance=Ci.getBarryCenter().distance(Ck.getBarryCenter())
                coeff=Fj.getMeasure()/Ci.getMeasure()/distance
                Rigidite.addValue(i,k,-coeff) # terme extradiagonal
            else:
                coeff=Fj.getMeasure()/Ci.getMeasure()/Ci.getBarryCenter().distance(Fj.getBarryCenter())
                #For the particular case where the mesh boundary does not coincide with the domain boundary
                x=Fj.getBarryCenter().x()
                y=Fj.getBarryCenter().y()
                RHS[i]+=coeff*sin(pi*x)*sin(pi*y)#mettre ici la condition limite du problème de Dirichlet
            Rigidite.addValue(i,i,coeff) # terme diagonal
    
    print("Linear system matrix building done")
    
    # Résolution du système linéaire
    #=================================
    LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"CG","LU")
    LS.setComputeConditionNumber()
    SolSyst=LS.solve()
    
    print "Preconditioner used : ", LS.getNameOfPc()
    print "Number of iterations used : ", LS.getNumberOfIter()
    print("Linear system solved")
    
    test_desc["Linear_solver_algorithm"]=LS.getNameOfMethod()
    test_desc["Linear_solver_preconditioner"]=LS.getNameOfPc()
    test_desc["Linear_solver_precision"]=LS.getTolerance()
    test_desc["Linear_solver_maximum_iterations"]=LS.getNumberMaxOfIter()
    test_desc["Linear_system_max_actual_iterations_number"]=LS.getNumberOfIter()
    test_desc["Linear_system_max_actual_error"]=LS.getResidu()
    test_desc["Linear_system_max_actual_condition number"]=LS.getConditionNumber()

    # Création du champ résultat
    #===========================
    my_ResultField = cdmath.Field("ResultField", cdmath.CELLS, my_mesh, 1)
    for i in range(nbCells):
        my_ResultField[i]=SolSyst[i];
    #sauvegarde sur le disque dur du résultat dans un fichier paraview
    my_ResultField.writeVTK("FiniteVolumes2DPoisson_SQUARE_"+meshType+str(nbCells))
    
    print("Numerical solution of 2D Poisson equation on a square using finite elements done")
    
    end = time.time()

    #Calcul de l'erreur commise par rapport à la solution exacte
    #===========================================================
    #The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi
    max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(2*pi*pi)
    max_sol_num=my_ResultField.max()
    min_sol_num=my_ResultField.min()
    erreur_abs=0
    for i in range(nbCells) :
        if erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) :
            erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i])
    
    print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
    print ("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)

    #Postprocessing :
    #================
	# Extraction of the diagonal data
    diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,1,0],[1,0,0], resolution)
    # save 2D picture
    PV_routines.Save_PV_data_to_picture_file("FiniteVolumes2DPoisson_SQUARE_"+meshType+str(nbCells)+'_0.vtu',"ResultField",'CELLS',"FiniteVolumes2DPoisson_SQUARE_"+meshType+str(nbCells))

    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Absolute_error"]=erreur_abs
    test_desc["Relative_error"]=erreur_abs/max_abs_sol_exacte

    with open('test_Poisson'+str(my_mesh.getMeshDimension())+'D_VF_'+str(nbCells)+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return erreur_abs/max_abs_sol_exacte, nbCells, diag_data, min_sol_num, max_sol_num, end - start


def solve_file( filename,resolution, meshType, testColor):
    my_mesh = cdmath.Mesh(filename+".med")
    return solve(my_mesh, filename,resolution, meshType, testColor)
    
if __name__ == """__main__""":
        mesh51 = cdmath.Mesh(0,1,51,0,1,51)
        solve(mesh51,'51',100,"Regular_squares","Green")
