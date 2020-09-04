# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson -\triangle u = f sur un cube avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage quelconque
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#               Comparaison de la solution num"rique avec la solution exacte u=sin(pi*x)*sin(pi*y)*sin(pi*z)
#================================================================================================================================

import cdmath
import time, json
import PV_routines
import VTK_routines
from math import sin, pi

test_desc={}
test_desc["Initial_data"]="None"
test_desc["Boundary_conditions"]="Dirichlet"
test_desc["Global_name"]="FV simulation of the 3D Poisson equation"
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

def solve(my_mesh, filename,resolution, meshType, testColor):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    
    # Maillage du domaine cubique [0,1]x[0,1]x[0,1], définition des bords
    #====================================================================================
    xmin=0
    xmax=1
    ymin=0
    ymax=1
    zmin=0
    zmax=1
    
    eps=1e-6
    my_mesh.setGroupAtPlan(0,0,eps,"DirichletBorder")#Bord GAUCHE
    my_mesh.setGroupAtPlan(1,0,eps,"DirichletBorder")#Bord DROIT
    my_mesh.setGroupAtPlan(0,1,eps,"DirichletBorder")#Bord BAS
    my_mesh.setGroupAtPlan(1,1,eps,"DirichletBorder")#Bord HAUT
    my_mesh.setGroupAtPlan(0,2,eps,"DirichletBorder")#Bord AVANT
    my_mesh.setGroupAtPlan(1,2,eps,"DirichletBorder")#Bord ARRIERE
    
    if( my_mesh.getSpaceDimension()!=3 or my_mesh.getMeshDimension()!=3) :
        raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")

    nbCells = my_mesh.getNumberOfCells()
    
    test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
    test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
    test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfCells()
    test_desc["Mesh_cell_type"]=my_mesh.getElementTypesNames()

    print("Mesh groups done")
    print("Number of cells =", nbCells)
    
    #Discrétisation du second membre et extraction du nb max de voisins d'une cellule
    #================================================================================
    my_RHSfield = cdmath.Field("RHS_field", cdmath.CELLS, my_mesh, 1)
    maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
    #parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
    for i in range(nbCells): 
        Ci = my_mesh.getCell(i)
        x = Ci.x()
        y = Ci.y()
        z = Ci.z()
        my_RHSfield[i]=3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)#mettre la fonction definie au second membre de l edp
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
            Rigidite.addValue(i,i,coeff) # terme diagonal
    
    print("Linear system matrix building done")
    
    # Résolution du système linéaire
    #=================================
    LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"CG","LU")
    LS.setComputeConditionNumber()
    SolSyst=LS.solve()
    
    print("Preconditioner used : ", LS.getNameOfPc() )
    print("Number of iterations used : ", LS.getNumberOfIter() )
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
    my_ResultField.writeVTK("FiniteVolumes3DPoisson_CUBE_"+meshType+str(nbCells))
    
    print("Numerical solution of 3D Poisson equation on a cube using finite elements done")
    
    end = time.time()

    #Calcul de l'erreur commise par rapport à la solution exacte
    #===========================================================
    #The following formulas use the fact that the exact solution is equal to the right hand side divided by 3*pi*pi
    max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(3*pi*pi)
    max_sol_num=my_ResultField.max()
    min_sol_num=my_ResultField.min()
    erreur_abs=0
    for i in range(nbCells) :
        if erreur_abs < abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i]) :
            erreur_abs = abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i])
    
    print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
    print("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)

    # Postprocessing 
    #==============
    # Extraction of the diagonal data
    diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,0,0],[1,1,1], resolution)
    # save 3D picture
    VTK_routines.Clip_VTK_data_to_VTK("FiniteVolumes3DPoisson_CUBE_"+meshType+str(nbCells)+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "FiniteVolumes3DPoisson_CUBE_"+meshType+str(nbCells)+'_0.vtu',[0.5,0.5,0.5], [-0.5,-0.5,-0.5],resolution )
    PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"FiniteVolumes3DPoisson_CUBE_"+meshType+str(nbCells)+'_0.vtu',"ResultField",'CELLS',"Clip_VTK_data_to_VTK_"+"FiniteVolumes3DPoisson_CUBE_"+meshType+str(nbCells))

    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Absolute_error"]=erreur_abs
    test_desc["Relative_error"]=erreur_abs/max_abs_sol_exacte

    with open('test_Poisson'+str(my_mesh.getMeshDimension())+'D_VF_CUBE_'+meshType+str(nbCells)+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return erreur_abs/max_abs_sol_exacte, nbCells, diag_data, min_sol_num, max_sol_num, end - start


def solve_file( filename,resolution, meshType, testColor):
    my_mesh = cdmath.Mesh(filename+".med")
    return solve(my_mesh, filename,resolution, meshType, testColor)
    
if __name__ == """__main__""":
        mesh51 = cdmath.Mesh(0,1,51,0,1,51,0,1,51)
        solve(mesh51,'51',100,"Regular_cubes","Green")
