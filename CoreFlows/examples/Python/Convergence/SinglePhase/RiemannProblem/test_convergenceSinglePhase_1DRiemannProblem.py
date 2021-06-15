import SinglePhase_1DRiemannProblem
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

    
def test_validationSinglePhase_1DRiemannProblem(cfl,isExplicit,scheme):
    start = time.time()
    #### 1D regular grid
    meshList=[10,20,50,100,200, 400]
    meshType="1D regular grid"
    testColor="Green"
    nbMeshes=len(meshList)
    mesh_size_tab=meshList
    mesh_name='RegularGrid'

    a=0.  ;  b=1.
    x=[0]*nbMeshes
    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    error_T_tab=[0]*nbMeshes
    sol_p=[0]*nbMeshes
    sol_u=[0]*nbMeshes
    sol_T=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes

    plt.close('all')
    i=0

    # Storing of numerical errors, mesh sizes and solution
    for nx in meshList:
        sol_u[i], sol_p[i], sol_T[i], error_u_tab[i], error_p_tab[i], error_T_tab[i], x[i], time_tab[i] = SinglePhase_1DRiemannProblem.solve(a,b,nx,cfl,isExplicit, scheme)
        error_p_tab[i]=log10(error_p_tab[i])
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    if(isExplicit):
        ExplicitOrImplicit="Explicit"
    else:
        ExplicitOrImplicit="Implicit"

    # Plot of results
    for i in range(nbMeshes):
            plt.plot(x[i], sol_p[i], label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('Position (m)')
    plt.ylabel('Pressure (bar)')
    plt.title('Plot of pressure in 1D Euler system \n with '+ExplicitOrImplicit+scheme+' scheme')
    plt.savefig(mesh_name+'_1DEulerSystem'+scheme+'_Pressure.png')
    plt.close()

    plt.clf()
    for i in range(nbMeshes):
            plt.plot(x[i], sol_u[i], label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('Position (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Plot of velocity in 1D Euler system \n with '+ExplicitOrImplicit+scheme+' scheme')
    plt.savefig(mesh_name+'_1DEulerSystem'+scheme+'_Velocity.png')
    plt.close()

    plt.clf()
    for i in range(nbMeshes):
            plt.plot(x[i], sol_T[i], label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('Position (m)')
    plt.ylabel('Temperature (K)')
    plt.title('Plot of temperature in 1D Euler system \n with '+ExplicitOrImplicit+scheme+' scheme')
    plt.savefig(mesh_name+'_1DEulerSystem'+scheme+'_Temperature.png')
    plt.close()

    for i in range(nbMeshes):
        mesh_size_tab[i]=log10(mesh_size_tab[i])
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validationSinglePhase_1DRiemannProblem() : Make sure you use distinct meshes and at least two meshes'

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    a=( a3*b1u-a2*b2u)/det
    b=(-a2*b1u+a1*b2u)/det
    
    print( ExplicitOrImplicit + scheme+" scheme for Euler equation on 1D regular grid : scheme order is ", -a)
    
    assert abs(a+0.26 )<0.01
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='log(|error pressure|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='straight line with slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(|error p|)')
    plt.title('Convergence of finite volumes for the Euler equation \n with '+ExplicitOrImplicit+scheme+' scheme on a 1D regular grid (pressure)')

    plt.savefig(mesh_name+"SinglePhase_1DRiemannProblem_"+scheme+ExplicitOrImplicit+"_CFL"+str(cfl)+"_ConvergenceCurve_pressure.png")
    
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='log(|error velocity|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='straight line with slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(|error u|)')
    plt.title('Convergence of finite volumes for the Euler equation \n with '+ExplicitOrImplicit+scheme+' scheme on a 1D regular grid (velocity)')

    plt.savefig(mesh_name+"SinglePhase_1DRiemannProblem_"+scheme+ExplicitOrImplicit+"_CFL"+str(cfl)+"_ConvergenceCurve_velocity.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the Euler equation \n with '+ExplicitOrImplicit+scheme+' scheme on a 1D regular grid')

    plt.savefig(mesh_name+"SinglePhase_1DRiemannProblem_"+scheme+ExplicitOrImplicit+"_CFL"+str(cfl)+"_ComputationalTime.png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Euler_Equation"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    convergence_synthesis["Numerical_method_name"]=scheme+" scheme"
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]=ExplicitOrImplicit
    convergence_synthesis["Initial_data"]="Riemann problem"
    convergence_synthesis["Boundary_conditions"]="Periodic"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="1D regular grid"
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_SinglePhase_1DRiemannProblem'+ExplicitOrImplicit+'_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >3 :
        cfl = float(sys.argv[1])
        isExplicit = bool(int(sys.argv[2]))
        scheme = sys.argv[3]
        test_validationSinglePhase_1DRiemannProblem(cfl,isExplicit, scheme)
    else :
        test_validationSinglePhase_1DRiemannProblem(0.99,True, "Upwind")
