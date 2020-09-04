# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson 2D -\triangle u = 0 sur un disque avec conditions aux limites de Dirichlet discontinues
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage triangulaire
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=atan(2*x/(x*x+y*y-1))
#================================================================================================================================

import cdmath
from math import atan, pi
import time, json
import PV_routines
import VTK_routines

test_desc={}
test_desc["Initial_data"]="None"
test_desc["Boundary_conditions"]="Dirichlet"
test_desc["Global_name"]="FV simulation of the 2D Poisson equation"
test_desc["Global_comment"]="Condition limite discontinue, Maillage triangulaire"
test_desc["PDE_model"]="Poisson"
test_desc["PDE_is_stationary"]=True
test_desc["PDE_search_for_stationary_solution"]=False
test_desc["Numerical_method_name"]="VF5"
test_desc["Numerical_method_space_discretization"]="Finite elements"
test_desc["Numerical_method_time_discretization"]="None"
test_desc["Mesh_is_unstructured"]=True
test_desc["Geometry"]="Disk"
test_desc["Part_of_mesh_convergence_analysis"]=True

def solve(my_mesh,filename,resolution, meshName, testColor):
    start = time.time()
    test_desc["Mesh_name"]=meshName
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
    my_ExactSol = cdmath.Field("Exact_field", cdmath.CELLS, my_mesh, 1)
    maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
    eps=1e-6#For coarse meshes
    
    #parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
    for i in range(nbCells): 
        Ci = my_mesh.getCell(i)
        x = Ci.x()
        y = Ci.y()

        #Robust calculation of atan(2x/(x**2+y**2-1)
        if x**2+y**2-1 > eps :
            print("!!! Warning Mesh ",meshName," !!! Cell is not in the unit disk.",", eps=",eps, ", x**2+y**2-1=",x**2+y**2 - 1)
            #raise ValueError("Exact solution computation !!! Domain should be the unit disk.")
        if x**2+y**2-1 < -eps :
            my_ExactSol[i] = atan(2*x/(x**2+y**2-1))
        elif x>0 : #x**2+y**2-1>=0
            my_ExactSol[i] = -pi/2
        elif x<0 : #x**2+y**2-1>=0
            my_ExactSol[i] =  pi/2
        else : #x=0
            my_ExactSol[i] = 0

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
                Rigidite.setValue(i,k,-coeff) # terme extradiagonal
            else:
                coeff=Fj.getMeasure()/Ci.getMeasure()/Ci.getBarryCenter().distance(Fj.getBarryCenter())
                #For the particular case where the mesh boundary does not coincide with the domain boundary
                x=Fj.getBarryCenter().x()
                y=Fj.getBarryCenter().y()
                if x**2+y**2-1 > eps :
                    print("!!! Warning Mesh ", meshName," !!! Face is not in the unit disk.",", eps=",eps, ", x**2+y**2-1=",x**2+y**2 - 1)
                    #raise ValueError("!!! Domain should be the unit disk.")
                if x**2+y**2-1 < -eps :
                    RHS[i]+= coeff*atan(2*x/(x**2+y**2-1))
                elif x>0 : #x**2+y**2-1>=0
                    RHS[i]+= coeff*(-pi/2)
                elif x<0 : #x**2+y**2-1>=0
                    RHS[i]+= coeff*pi/2
                else : #x=0
                    RHS[i]+=  0
            Rigidite.addValue(i,i,coeff) # terme diagonal
    
    print("Linear system matrix building done")
    
    # Résolution du système linéaire
    #=================================
    LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"CG","LU")
    LS.setComputeConditionNumber()
    SolSyst=LS.solve()
    
    print( "Preconditioner used : ", LS.getNameOfPc() )
    print( "Number of iterations used : ", LS.getNumberOfIter() )
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
    my_ResultField.writeVTK("FiniteVolumes2DPoissonStiffBC_DISK_"+meshName+str(nbCells))
    
    print("Numerical solution of 2D Poisson equation on a disk using finite volumes done")
    
    end = time.time()

    #Calcul de l'erreur commise par rapport à la solution exacte
    #===========================================================
    l2_norm_sol_exacte=my_ExactSol.normL2()[0]
    l2_error = (my_ExactSol - my_ResultField).normL2()[0]
    
    print("L2 relative error = norm( exact solution - numerical solution )/norm( exact solution ) = ", l2_error/l2_norm_sol_exacte)
    print("Maximum numerical solution = ", my_ResultField.max(), " Minimum numerical solution = ", my_ResultField.min())
    print("Maximum exact solution = ", my_ExactSol.max(), " Minimum exact solution = ", my_ExactSol.min())

    #Postprocessing :
    #================
	# Extraction of the diagonal data
    diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,-1,0],[0,1,0], resolution)
    # save 2D picture
    PV_routines.Save_PV_data_to_picture_file("FiniteVolumes2DPoissonStiffBC_DISK_"+meshName+str(nbCells)+'_0.vtu',"ResultField",'CELLS',"FiniteVolumes2DPoissonStiffBC_DISK_"+meshName+str(nbCells))

    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Absolute_error"]=l2_error
    test_desc["Relative_error"]=l2_error/l2_norm_sol_exacte

    with open('test_PoissonStiffBC'+str(my_mesh.getMeshDimension())+'D_VF_'+"DISK_"+str(nbCells)+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return l2_error/l2_norm_sol_exacte, nbCells, diag_data, my_ResultField.min(), my_ResultField.max(), end - start


def solve_file( filename,resolution, meshName, testColor):
    my_mesh = cdmath.Mesh(filename+".med")
    return solve(my_mesh, filename,resolution, meshName, testColor)
    
if __name__ == """__main__""":
    solve("diskWithTriangles",100,"Unstructured_triangles","Green")
