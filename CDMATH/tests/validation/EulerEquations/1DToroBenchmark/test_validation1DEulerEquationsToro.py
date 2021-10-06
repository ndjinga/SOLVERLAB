import EulerEquations1D
import exact_rs_stiffenedgas

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

    
def test_validation1DEulerEquationsToro(scheme,isImplicit):
    start = time.time()
    #### 1D regular grid
    meshList=[50,100,200,400]
    meshType="1D regular grid"
    testColor="Orange : à debugger"
    nbMeshes=len(meshList)
    mesh_name='RegularGrid'
    simulation_name='Riemann problem'
    a=0.  ;  b=1.
    error_rho_tab=[0]*nbMeshes
    error_v_tab  =[0]*nbMeshes
    error_p_tab  =[0]*nbMeshes
    error_q_tab  =[0]*nbMeshes
    error_h_tab  =[0]*nbMeshes
    error_e_tab  =[0]*nbMeshes
    time_tab     =[0]*nbMeshes

    plt.close('all')

    # Initial parameters for Riemann problem (p in Pa, v in m/s, rho in Kg/m**3)
    
    p_L = 460.894
    p_R = 46.095
    v_L = 19.5975
    v_R = -6.19633
    rho_L = 5.99924
    rho_R = 5.99242    
    
    WL=[rho_L, v_L, p_L]
    WR=[rho_R, v_R, p_R]

    pos_disc = 0.4
    
    precision = 1.e-5
    
    # Stifened gas equation of state
    gamma = 1.4
    p0 = 0.
    q = 0
    
    tmax=0.035
    lam_max= max(abs(v_L)+sqrt(gamma*p_L/rho_L), abs(v_R)+sqrt(gamma*p_R/rho_R) )

    ###### Compute exact solution #####
    RS = exact_rs_stiffenedgas.exact_rs_stiffenedgas(gamma, gamma, p0, p0);# Set the problem EOS
    RS.solve_RP(WL,WR);# Solve the problem structure

    numsamples=max(meshList) #nombre d'échantillons = taille du plus grand maillage
    dx=(b-a)/numsamples
    
    rho_field_exact=[0]*numsamples
    v_field_exact  =[0]*numsamples
    p_field_exact  =[0]*numsamples
    q_field_exact  =[0]*numsamples
    h_field_exact  =[0]*numsamples
    e_field_exact  =[0]*numsamples

    for i in range(numsamples):
    
        soln = RS.sample_solution(WL, WR, ( a+(i+1/2)*dx - pos_disc)/tmax); # Determine the solution ar time t in cell number i
    
        rho_field_exact[i]=soln[0]
        v_field_exact[i]  =soln[1]
        p_field_exact[i]  =soln[2]
        q_field_exact[i]  =rho_field_exact[i]*v_field_exact[i]
        h_field_exact[i]  =exact_rs_stiffenedgas.stiffenedgas_h(rho_field_exact[i], p_field_exact[i], gamma, p0)
        e_field_exact[i]  =exact_rs_stiffenedgas.stiffenedgas_e(rho_field_exact[i], p_field_exact[i], gamma, p0)

    #Set convergence picture parameters    
    max_initial_rho = max(rho_L,rho_R)
    min_initial_rho = min(rho_L,rho_R)
    max_initial_p = max(p_L,p_R)
    min_initial_p = min(p_L,p_R)
    max_initial_v = max(v_L,v_R)
    min_initial_v = min(v_L,v_R)
    max_initial_q = max_initial_rho*max_initial_v
    min_initial_q = min_initial_rho*min_initial_v

    e_L=exact_rs_stiffenedgas.p_to_e_StiffenedGaz(p_L, rho_L, gamma, p0)
    e_R=exact_rs_stiffenedgas.p_to_e_StiffenedGaz(p_R, rho_R, gamma, p0)
    h_L=e_L+p_L/rho_L
    h_R=e_R+p_R/rho_R
    max_initial_e = max(e_L, e_R)
    min_initial_e = min(e_L, e_R)
    min_initial_h = min(h_L,h_R)
    max_initial_h = max(h_L,h_R)

    fig, ([axDensity, axMomentum, axEnthalpie],[axPressure, axVitesse, axEinterne]) = plt.subplots(2, 3,sharex=True, figsize=(14,10))
    plt.gcf().subplots_adjust(wspace = 0.5)

    axDensity.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], rho_field_exact, label='Exact solution')
    axDensity.set(xlabel='x (m)', ylabel='Densité (Kg/m3)')
    axDensity.set_xlim(a,b)
    axDensity.set_ylim(0.9*min_initial_rho , 6*max_initial_rho )
    axDensity.legend()

    axMomentum.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], q_field_exact, label='Exact solution')
    axMomentum.set(xlabel='x (m)', ylabel='Momentum (Kg/m2/s)')
    axMomentum.set_xlim(a,b)
    axMomentum.set_ylim(0.9*min_initial_q , 2.5*max_initial_q )
    axMomentum.legend()
    
    axEnthalpie.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], h_field_exact, label='Exact solution')
    axEnthalpie.set(xlabel='x (m)', ylabel='Enthalpy (J/Kg)')
    axEnthalpie.set_xlim(a,b)
    axEnthalpie.set_ylim(0.9*min_initial_h , 1.75*max_initial_h )
    axEnthalpie.legend()
    
    axPressure.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], p_field_exact, label='Exact solution')
    axPressure.set(xlabel='x (m)', ylabel='Pression (Pa)')
    axPressure.set_xlim(a,b)
    axPressure.set_ylim(0.9*min_initial_p , 4*max_initial_p)
    axPressure.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axPressure.legend()

    axVitesse.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], v_field_exact, label='Exact solution')
    axVitesse.set(xlabel='x (m)', ylabel='Vitesse (m/s)')
    axVitesse.set_xlim(a,b)
    axVitesse.set_ylim(0.9*min_initial_v , 1.1*max_initial_v)
    axVitesse.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axVitesse.legend()

    axEinterne.plot(  [a+0.5*dx + i*dx for i in range(numsamples)], e_field_exact, label='Exact solution')
    axEinterne.set(xlabel='x (m)', ylabel='Energie interne (J/Kg)')
    axEinterne.set_xlim(a,b)
    axEinterne.set_ylim(0.9*min_initial_e , 1.75*max_initial_e)
    axEinterne.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axEinterne.legend()

    ### Main loop to Store of numerical errors, mesh sizes and solutions
    j=0
    for nx in meshList:
        p_field, v_field, e_field, rho_field, q_field, h_field, time_tab[j] = EulerEquations1D.solve(a, b, nx, mesh_name, meshType, isImplicit, scheme, simulation_name, testColor,fig,p_L,v_L,rho_L,p_R,v_R,rho_R, pos_disc, tmax, lam_max, gamma, p0)

        axDensity.plot(  [a+0.5*dx + i*dx for i in range(nx)], rho_field, label=scheme+' scheme on '+str(nx)+' cells') #new picture for video # Returns a tuple of line objects, thus the comma
        axMomentum.plot( [a+0.5*dx + i*dx for i in range(nx)],   q_field, label=scheme+' scheme on '+str(nx)+' cells')    
        axEnthalpie.plot([a+0.5*dx + i*dx for i in range(nx)],   h_field, label=scheme+' scheme on '+str(nx)+' cells')    
        axPressure.plot( [a+0.5*dx + i*dx for i in range(nx)],   p_field, label=scheme+' scheme on '+str(nx)+' cells')
        axVitesse.plot(  [a+0.5*dx + i*dx for i in range(nx)],   v_field, label=scheme+' scheme on '+str(nx)+' cells')
        axEinterne.plot( [a+0.5*dx + i*dx for i in range(nx)],   e_field, label=scheme+' scheme on '+str(nx)+' cells')

        for i in range(nx):
            soln = RS.sample_solution(WL, WR, ( a+(i+1/2)*dx - pos_disc)/tmax); # Determine the solution at time t in cell number i

            error_rho_tab[j]=max(error_rho_tab[j],abs(rho_field[i]-soln[0]))
            error_v_tab[j]  =max(error_v_tab[j],  abs(  v_field[i]-soln[1]))
            error_p_tab[j]  =max(error_p_tab[j],  abs(  p_field[i]-soln[2]))
        
        j=j+1
    
    end = time.time()

    for i in range(nbMeshes):
        error_rho_tab[i]=log10(error_rho_tab[i])
        error_v_tab[i]  =log10(error_v_tab[i])
        error_p_tab[i]  =log10(error_p_tab[i])
        time_tab[i]     =log10(time_tab[i])

    if(isImplicit):
        ImplicitOrExplicit="Implicit"    
    else:
        ImplicitOrExplicit="Explicit"    

    # Plot of the 6 results fields
    plt.suptitle('Euler equations : Convergence of ' + ImplicitOrExplicit + "_" + scheme + ' scheme ')
    plt.savefig(mesh_name+"_1DEulerEquations_Toro4_" + ImplicitOrExplicit + "_" + scheme + ".png")    
    plt.close()

        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(meshList,meshList)
    a2=np.sum(meshList)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation1DEulerEquation'+ImplicitOrExplicit+scheme+' : Make sure you use distinct meshes and at least two meshes'

    b1r=np.dot(error_rho_tab,meshList)   
    b2r=np.sum(error_rho_tab)
    ar=( a3*b1r-a2*b2r)/det
    br=(-a2*b1r+a1*b2r)/det
    
    print(ImplicitOrExplicit + "_" + scheme + " scheme for Euler Equations on 1D regular grid : order for density is ", -ar)
    
    b1u=np.dot(error_v_tab,meshList)   
    b2u=np.sum(error_v_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print(ImplicitOrExplicit + "_" + scheme + " scheme for Euler Equations on 1D regular grid : order for velocity is ", -au)
    
    b1p=np.dot(error_p_tab,meshList)   
    b2p=np.sum(error_p_tab)
    ap=( a3*b1p-a2*b2p)/det
    bp=(-a2*b1p+a1*b2p)/det
    
    print(ImplicitOrExplicit + "_" + scheme + " scheme for Euler Equations on 1D regular grid : order for velocity is ", -ap)
    
    # Plot of convergence curves
    fig_convergence, ([axDensity, axVitesse, axPressure]) = plt.subplots(1, 3,sharex=True, figsize=(14,5))
    plt.gcf().subplots_adjust(wspace = 0.5)

    plt.close()
    axDensity.plot(meshList, error_rho_tab, label='log(|error density|)')
    axDensity.plot(meshList, a*np.array(meshList)+b,label='least square slope : '+'%.3f' % ar)
    axDensity.legend()
    axDensity.set(xlabel='log(Number of cells)', ylabel='log(|error density|)')
    #axDensity.title('Mesh convergence of density')

    axVitesse.plot(meshList, error_v_tab, label='log(|error velocity|)')
    axVitesse.plot(meshList, a*np.array(meshList)+b,label='least square slope : '+'%.3f' % au)
    axVitesse.legend()
    axVitesse.set(xlabel='log(Number of cells)', ylabel='log(|error velocity|)')
    #axVitesse.title('Mesh convergence of Velocity')

    axPressure.plot(meshList, error_p_tab, label='log(|error Pressure|)')
    axPressure.plot(meshList, a*np.array(meshList)+b,label='least square slope : '+'%.3f' % ap)
    axPressure.legend()
    axPressure.set(xlabel='log(Number of cells)', ylabel='log(|error Pressure|)')
    #axPressure.title('Mesh convergence of Pressure')

    plt.suptitle('Convergence of finite volumes for the Euler equations : Toro 4 \n Riemann problem with '+ImplicitOrExplicit+scheme+' scheme on a 1D regular grid')
    plt.savefig(mesh_name+"1DEulerEquation_"+ImplicitOrExplicit+scheme+"_ConvergenceCurves.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(meshList, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the Euler equations : Toro 4 \n Riemann problem with '+ImplicitOrExplicit+scheme+' scheme on a 1D regular grid')
    plt.savefig(mesh_name+"1DEulerEquation_"+ImplicitOrExplicit+scheme+"_ComputationalTime.png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Euler_Equations"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    convergence_synthesis["Numerical_method_name"]="Upwind scheme"
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]=ImplicitOrExplicit
    convergence_synthesis["Initial_data"]="Riemann Problem"
    convergence_synthesis["Boundary_conditions"]="Neumann"
    convergence_synthesis["Space_dimension"]=1
    convergence_synthesis["Mesh_dimension"]=1
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=meshList
    convergence_synthesis["Mesh_cell_type"]="1D regular grid"
    convergence_synthesis["Scheme_order_density"]=-ar
    convergence_synthesis["Scheme_order_velocity"]=-au
    convergence_synthesis["Scheme_order_pressure"]=-ap
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_1DEulerEquations_Toro4_'+ImplicitOrExplicit+scheme+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >2 :
        scheme = sys.argv[1]
        isImplicit = bool(int(sys.argv[2]))
        test_validation1DEulerEquationsToro(scheme,isImplicit)
    else :
        test_validation1DEulerEquationsToro('Roe',False)

