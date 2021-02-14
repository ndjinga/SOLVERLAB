import FiniteVolumes2DPoisson_SQUARE
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json

convergence_synthesis=dict(FiniteVolumes2DPoisson_SQUARE.test_desc)

def test_validation2DVF_cross_triangles():
    start = time.time()
    ### 2D FV cross triangles mesh
    #meshList=[5,9,15,21,31]
    meshList=['squareWithCrossTriangles_0','squareWithCrossTriangles_1','squareWithCrossTriangles_2','squareWithCrossTriangles_3']#,'squareWithCrossTriangles_4'
    mesh_path='../../../ressources/2DCrossTriangles/'
    meshType="Regular_cross_triangles"
    testColor="Green"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_name='squareWithCrossTriangles'
    diag_data=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
    #for nx in meshList:
        #my_mesh=cdmath.Mesh(0,1,nx,0,1,nx*nx)#Use script provided by Adrien to split each quadrangle in 4 triangles
        #error_tab[i], mesh_size_tab[i], diag_data[i], min_sol_num, max_sol_num, time_tab[i] =FiniteVolumes2DPoisson_SQUARE.solve(my_mesh,str(nx)+'x'+str(nx),resolution,meshType,testColor)
        error_tab[i], mesh_size_tab[i], diag_data[i], min_sol_num, max_sol_num, time_tab[i] =FiniteVolumes2DPoisson_SQUARE.solve_file(mesh_path+filename,resolution,meshType,testColor)
        assert min_sol_num>-0.01 
        assert max_sol_num<1.4
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
    plt.title('Plot over diagonal line for finite volumes \n for Laplace operator on 2D cross triangles meshes')
    plt.savefig(mesh_name+"_2DPoissonFV_PlotOverDiagonalLine.png")

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DVF_cross_triangles() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print( "FV on 2D cross triangles mesh : scheme order is ", -a )
    assert abs(a+0.27)<0.001
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.plot(mesh_size_tab, error_tab)
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite volumes for \n Laplace operator on 2D cross triangles meshes')
    plt.savefig(mesh_name+"_2DPoissonFV_ConvergenceCurve.png")

    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(sqrt(number of cells))')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes \n for Laplace operator on 2D cross triangles meshes')
    plt.savefig(mesh_name+"_2DPoissonFV_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_cell_type"]="cross triangles"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_Poisson_2DVF_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

    import os
    os.system("jupyter-nbconvert --to notebook --execute Convergence_Poisson_FV5_SQUARE_cross_triangles.ipynb")
    os.system("jupyter-nbconvert --to html Convergence_Poisson_FV5_SQUARE_cross_triangles.ipynb")
    os.system("jupyter-nbconvert --to pdf Convergence_Poisson_FV5_SQUARE_cross_triangles.ipynb")

if __name__ == """__main__""":
    test_validation2DVF_cross_triangles()
