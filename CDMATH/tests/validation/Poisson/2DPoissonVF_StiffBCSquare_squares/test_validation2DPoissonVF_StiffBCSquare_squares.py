import FiniteVolumes2DPoissonStiffBC_SQUARE
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json

convergence_synthesis=dict(FiniteVolumes2DPoissonStiffBC_SQUARE.test_desc)

def test_validation2DVF_StiffBC_squares():
    start = time.time()
    #### 2D FV square mesh of a square
    meshList=['squareWithSquares_1','squareWithSquares_2','squareWithSquares_3','squareWithSquares_4','squareWithSquares_5']
    meshType="Regular_squares"
    testColor="Green"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../../../ressources/2DCartesien/'
    mesh_name='squareWithSquares'
    diag_data=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    max_tab=[0]*nbMeshes
    min_tab=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_tab[i], mesh_size_tab[i], diag_data[i], min_tab[i], max_tab[i], time_tab[i] =FiniteVolumes2DPoissonStiffBC_SQUARE.solve_file(mesh_path+filename, resolution, mesh_name, testColor)
        assert min_tab[i]>-1.6
        assert max_tab[i]<1.6
        plt.plot(curv_abs, diag_data[i], label= str(mesh_size_tab[i]) + ' nodes')
        error_tab[i]=log10(error_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    # Plot over diagonal line
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Value on diagonal line')
    plt.title('Plot over diagonal line for finite volumes for Laplace operator \n on 2D square with square meshes')
    plt.savefig(mesh_name+"_2DPoissonVF_StiffBC_PlotOverDiagonalLine.png")

    
    # Plot of min and max curves
    plt.close()
    plt.plot(mesh_size_tab, min_tab, label='Minimum value')
    plt.plot(mesh_size_tab, max_tab, label='Maximum value')
    plt.legend()
    plt.xlabel('Number of nodes')
    plt.ylabel('Value')
    plt.title('Min/Max of Finite volumes for Laplace operator \n on 2D square with square meshes')
    plt.savefig(mesh_name+"_2DPoissonVF_StiffBC_MinMax.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i] = 0.5*log10(mesh_size_tab[i])

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DVF_StiffBC_squares() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print( "FV on 2D square with square mesh : l2 scheme order is ", -a )
    assert abs(a+1.)<0.01

    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('log(sqrt(number of nodes))')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite volumes for Laplace operator \n on 2D square with square meshes')
    plt.savefig(mesh_name+"_2DPoissonVF_StiffBC_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(sqrt(number of nodes))')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for Laplace operator \n on 2D square with square meshes')
    plt.savefig(mesh_name+"_2DPoissonVF_StiffBC_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_cell_type"]="Squares"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_PoissonStiffBC_2DFV_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    test_validation2DVF_StiffBC_squares()
