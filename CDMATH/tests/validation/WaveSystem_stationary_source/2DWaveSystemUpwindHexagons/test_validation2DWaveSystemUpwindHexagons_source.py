import WaveSystemUpwind
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

def test_validation2DWaveSystemSourceUpwindHexagons(bctype,scaling):
    start = time.time()
    #### 2D hexagonal mesh
    meshList=['squareWithHexagons_1','squareWithHexagons_2','squareWithHexagons_3','squareWithHexagons_4','squareWithHexagons_5']
    meshType="Regular hexagons"
    testColor="Green"
    nbMeshes=len(meshList)
    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../../../ressources/2DHexagons/'
    mesh_name='squareWithHexagons'
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    t_final=[0]*nbMeshes
    ndt_final=[0]*nbMeshes
    max_vel=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)

    plt.close('all')
    i=0
    cfl=0.5
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], t_final[i], ndt_final[i], max_vel[i], diag_data_press[i], diag_data_vel[i], time_tab[i] =WaveSystemUpwind.solve_file(mesh_path+filename, mesh_name, resolution,scaling, meshType,testColor,cfl,bctype,True)
        assert max_vel[i]>0.001 and max_vel[i]<0.01
        error_p_tab[i]=log10(error_p_tab[i])
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    # Plot over diagonal line
    for i in range(nbMeshes):
        plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Pressure on diagonal line')
    plt.title('Plot over diagonal line for stationary Wave System with source term \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+'_Pressure_2DWaveSystemSourceUpwind_'+"PlotOverDiagonalLine.png")
    plt.close()

    plt.clf()
    for i in range(nbMeshes):
        plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for the stationary Wave System with source term \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemSourceUpwind_"+"PlotOverDiagonalLine.png")    
    plt.close()

    # Plot of number of time steps
    plt.close()
    plt.plot(mesh_size_tab, ndt_final, label='Number of time step to reach stationary regime')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max time steps for stationary regime')
    plt.title('Number of times steps required for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSourceUpwind_"+"TimeSteps.png")
    
    # Plot of number of stationary time
    plt.close()
    plt.plot(mesh_size_tab, t_final, label='Time where stationary regime is reached')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max time for stationary regime')
    plt.title('Simulated time for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSourceUpwind_"+"TimeFinal.png")
    
    # Plot of number of maximal velocity norm
    plt.close()
    plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max velocity norm')
    plt.title('Maximum velocity norm for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSourceUpwind_"+"MaxVelNorm.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=0.5*log10(mesh_size_tab[i])
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemSourceUpwindHexagonsFV() : Make sure you use distinct meshes and at least two meshes'

    b1p=np.dot(error_p_tab,mesh_size_tab)   
    b2p=np.sum(error_p_tab)
    ap=( a3*b1p-a2*b2p)/det
    bp=(-a2*b1p+a1*b2p)/det
    
    print "FV upwind on 2D hexagonal meshes : scheme order for pressure is ", -ap

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "FV upwind on 2D hexagonal meshes : scheme order for velocity is ", -au
    
    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='|error on stationary pressure|')
    plt.legend()
    plt.xlabel('1/2 log(number of cells)')
    plt.ylabel('log(error p)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_Pressure_2DWaveSystemSourceUpwind_"+"ConvergenceCurve.png")
    
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='|error on stationary velocity|')
    plt.legend()
    plt.xlabel('1/2 log(number of cells)')
    plt.ylabel('log(error u)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemSourceUpwind_"+"ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('1/2 log(number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the stationary Wave System \n with upwind scheme on 2D hexagonal meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSourceUpwind_ComputationalTime.png")
    
    plt.close('all')
    
    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Wave system"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    convergence_synthesis["Numerical_method_name"]="Upwind"
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]="Implicit"
    convergence_synthesis["Initial_data"]="Constant pressure, divergence free velocity"
    convergence_synthesis["Boundary_conditions"]="Periodic"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="Hexagons"
    convergence_synthesis["Numerical_error_velocity"]=error_u_tab
    convergence_synthesis["Numerical_error_pressure"]=error_p_tab
    convergence_synthesis["Max_vel_norm"]=max_vel
    convergence_synthesis["Final_time"]=t_final  
    convergence_synthesis["Final_time_step"]=ndt_final  
    convergence_synthesis["Scheme_order"]=min(-au,-ap)
    convergence_synthesis["Scheme_order_vel"]=-au
    convergence_synthesis["Scheme_order_press"]=-ap
    convergence_synthesis["Scaling_preconditioner"]="None"
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_WaveSystemSource_2DFV_Upwind_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >2 :
        bctype = sys.argv[1]
        scaling = int(sys.argv[2])
        test_validation2DWaveSystemSourceUpwindHexagons(bctype,scaling)
    else :
        raise ValueError("test_validation2DWaveSystemSourceUpwindHexagons.py expects a mesh file name and a scaling parameter")
