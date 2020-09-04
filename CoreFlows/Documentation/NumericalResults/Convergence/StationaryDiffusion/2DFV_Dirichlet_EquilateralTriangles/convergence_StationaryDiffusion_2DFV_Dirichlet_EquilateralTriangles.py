#!/usr/bin/env python
# -*-coding:utf-8 -*

import validationStationaryDiffusionEquation
import cdmath as cm
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json, os

convergence_synthesis=dict(validationStationaryDiffusionEquation.test_desc)

# !!!!!!!!!! Warning : result is affected by the fact that the mesh if not strictly a partition of the domain. BC should be adapted to find an order 2 of convergence as is the case in CDMATH !!!!!!!!!!!!!

def convergence_StationaryDiffusion_2DFV_Dirichlet_EquilateralTriangles():
    start = time.time() 
    ### 2D FV Equilateral triangle meshes
    method = 'FV'
    BC = 'Dirichlet'
    meshList=['squareWithEquilateralTriangles5','squareWithEquilateralTriangles20','squareWithEquilateralTriangles50','squareWithEquilateralTriangles100','squareWithEquilateralTriangles200']
    mesh_path=os.environ['CDMATH_INSTALL']+'/share/meshes/2DEquilateralTriangles/'
    mesh_name='squareWithEquilateralTriangles'
    meshType="Equilateral_triangles"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    diag_data=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    testColor="Green"
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
		my_mesh=cm.Mesh(mesh_path+filename+".med")
		error_tab[i], mesh_size_tab[i], diag_data[i], min_sol_num, max_sol_num, time_tab[i] =validationStationaryDiffusionEquation.SolveStationaryDiffusionEquation(my_mesh,resolution,meshType,method,BC)
		assert min_sol_num>-1 
		assert max_sol_num<1.2
		plt.plot(curv_abs, diag_data[i], label= str(mesh_size_tab[i]) + ' cells')
		error_tab[i]=log10(error_tab[i])
		time_tab[i]=log10(time_tab[i])
		mesh_size_tab[i] = 0.5*log10(mesh_size_tab[i])
		i=i+1
    end = time.time()
        
   

    # Plot over diagonal line
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Value on diagonal line')
    plt.title('Plot over diagonal line for finite volumes for Poisson problem \n on 2D Equilateral triangles with Dirichlet BC')
    plt.savefig(mesh_name+"_2DFV_StatDiffusion_Dirichlet_PlotOverDiagonalLine.png")

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'convergence_StationaryDiffusion_2DFV_Dirichlet_EquilateralTriangles() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print "FV for diffusion on 2D Equilateral triangle meshes: scheme order is ", -a
    assert abs(a+1.97)<0.01
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.plot(mesh_size_tab, error_tab)
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite elements for Poisson problem \n on 2D equilateral triangles meshes with Dirichlet BC \n Warning : result is affected by the fact that the mesh if not strictly a partition of the domain')
    plt.savefig(mesh_name+"_2DFV_StatDiffusion_Dirichlet_ConvergenceCurve.png")

    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for Poisson problem \n on 2D Equilateral triangles meshes with Dirichlet BC')
    plt.savefig(mesh_name+"_2DFV_StatDiffusion_Dirichlet_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    #convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dim"]=2
    convergence_synthesis["Mesh_dim"]=2
    convergence_synthesis["Mesh_cell_type"]="Triangles"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=round(-a,4)
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["PDE_model"]='Poisson'
    convergence_synthesis["Num_method"]=method
    convergence_synthesis["Bound_cond"]=BC
    convergence_synthesis["Comput_time"]=round(end-start,3)

    with open('Convergence_Poisson_2DFV_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

   
if __name__ == """__main__""":
    convergence_StationaryDiffusion_2DFV_Dirichlet_EquilateralTriangles()
