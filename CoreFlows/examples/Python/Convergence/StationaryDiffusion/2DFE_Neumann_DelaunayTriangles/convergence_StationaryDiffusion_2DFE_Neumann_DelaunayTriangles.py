#!/usr/bin/env python
# -*-coding:utf-8 -*

import validationStationaryDiffusionEquation
import cdmath as cm
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json, os

convergence_synthesis=dict(validationStationaryDiffusionEquation.test_desc)

def convergence_StationaryDiffusion_2DFE_Neumann_DelaunayTriangles():
    start = time.time() 
    ### 2D FE Delaunay triangle meshes
    method = 'FE'
    BC = 'Neumann'
    meshList=['squareWithTriangles_1','squareWithTriangles_2','squareWithTriangles_3','squareWithTriangles_4','squareWithTriangles_5']
    mesh_path=os.environ['CDMATH_INSTALL']+'/share/meshes//2DTriangles/'
    mesh_name='squareWithDelaunayTriangles'
    meshType="Unstructured_triangles"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    diag_data=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    testColor="Orange (not order 2), singular matrix"
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
    plt.title('Plot over diagonal line for finite elements for Poisson problem \n on 2D Delaunay triangles with Neumann BC')
    plt.savefig(mesh_name+"_2DFE_StatDiffusion_Neumann_PlotOverDiagonalLine.png")

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'convergence_StationaryDiffusion_2DFE_Neumann_DelaunayTriangles() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print "FE for diffusion on 2D Delaunay triangle meshes: scheme order is ", -a
    assert abs(a+0.82)<0.02
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.plot(mesh_size_tab, error_tab)
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite elements for Poisson problem \n on 2D Delaunay triangles meshes with Neumann BC')
    plt.savefig(mesh_name+"_2DFE_StatDiffusion_Neumann_ConvergenceCurve.png")

    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite elements for Poisson problem \n on 2D Delaunay triangles meshes with Neumann BC')
    plt.savefig(mesh_name+"_2DFE_StatDiffusion_Neumann_ComputationalTime.png")
    
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

    with open('Convergence_Poisson_2DFE_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

   
if __name__ == """__main__""":
    convergence_StationaryDiffusion_2DFE_Neumann_DelaunayTriangles()
