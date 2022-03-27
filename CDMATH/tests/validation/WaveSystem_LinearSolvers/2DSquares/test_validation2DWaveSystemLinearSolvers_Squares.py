import WaveSystemLinearSolvers
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
import time, json

    
def test_validation2DWaveSystem_LinearSolver_squares(num_scheme, lin_solver,cfl, mesh_folder):
    start = time.time()
    #### 2D squares meshes
    meshList=['squareWithSquares_1','squareWithSquares_2','squareWithSquares_3','squareWithSquares_4','squareWithSquares_5','squareWithSquares_6']
    meshType="Regular_squares"
    mesh_path='/2DCartesien/'
    mesh_name='squareWithSquares'
    testColor="Green"

    nbMeshes=len(meshList)
    mesh_size_tab=[0]*nbMeshes
    iterations_tab=[0]*nbMeshes
    time_tab=[0]*nbMeshes

    plt.close('all')
    i=0
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        iterations_tab[i], time_tab[i], mesh_size_tab[i] = WaveSystemLinearSolvers.solve_file(mesh_folder+mesh_path+filename, mesh_name, num_scheme, lin_solver,cfl)
        i=i+1
    
    end = time.time()

    # Plot of number of iteration numbers
    plt.close()
    plt.plot(mesh_size_tab, iterations_tab, label='Number of iterations ' + lin_solver)
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Number of iterations')
    plt.title('Number of iterations of finite volumes for the Wave System \n with ' +num_scheme+' scheme, '+lin_solver+' solver, cfl='+str(cfl)+', '+ mesh_name)
    plt.savefig("2DWaveSystemIterations_"+num_scheme+ 'scheme'  +lin_solver+'CFL'+str(cfl)+mesh_name+".png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='cpu time '+ lin_solver)
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('cpu time')
    plt.title('Computational time of finite volumes for the Wave System \n with '+num_scheme+' scheme, '+lin_solver+' solver, cfl='+str(cfl)+', '+ mesh_name)
    plt.savefig("2DWaveSystemComputationalTime_"+num_scheme+ 'scheme'  +lin_solver+'CFL'+str(cfl)+mesh_name+".png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Wave system"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=False
    convergence_synthesis["Numerical_method_name"]=num_scheme			
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]="Implicit"
    convergence_synthesis["Initial_data"]="Constant pressure, divergence free velocity"
    convergence_synthesis["Boundary_conditions"]="Wall"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="Squares"
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

if __name__ == """__main__""":
    if len(sys.argv) >4 :
        cfl = int(sys.argv[1])
        num_scheme = sys.argv[2]
        lin_solver = sys.argv[3]
        mesh_folder = sys.argv[4]
        test_validation2DWaveSystem_LinearSolver_squares(num_scheme, lin_solver,cfl, mesh_folder)
    else :
        raise ValueError("test_validation2DWaveSystem_LinearSolver_squares.py expects a numerical method, a linear solver, a cfl number, and a mesh folder path")

