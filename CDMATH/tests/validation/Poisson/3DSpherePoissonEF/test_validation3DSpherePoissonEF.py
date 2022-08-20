import FiniteElements3DPoissonSphere
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json

convergence_synthesis=dict(FiniteElements3DPoissonSphere.test_desc)
def test_validation3DSphereEF():
    start = time.time()
    #### 3D sphere FE triangle mesh
    meshList=['meshSphere_1','meshSphere_2','meshSphere_3','meshSphere_4','meshSphere_5']
    meshType="Unstructured_3D_triangles"
    testColor="Green"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    iterations_tab=[0]*nbMeshes
    residual_tab=[0]*nbMeshes
    diameter_tab=[0]*nbMeshes
    mesh_path='./'
    mesh_name='SphereWithTriangles'
    diag_data=[0]*nbMeshes
    resolution=100
    plt.close('all')
    i=0
    # Storing of numerical errors and mesh sizes
    for filename in meshList:
        error_tab[i], mesh_size_tab[i], min_sol_num, max_sol_num, time_tab[i], iterations_tab[i], residual_tab[i], diameter_tab[i] =FiniteElements3DPoissonSphere.solve(mesh_path+filename, resolution,meshType,testColor)
        assert min_sol_num>-1.1 
        assert max_sol_num<1.1
        error_tab[i]=log10(error_tab[i])
        time_tab[i]=log10(time_tab[i])
        diameter_tab[i]=log10(diameter_tab[i])
        residual_tab[i]=-log10(residual_tab[i])
        with open('./FiniteElementsOnSpherePoisson_PlotOnSortedLines'+meshType+str(mesh_size_tab[i])+'.csv') as f:
            lines = f.readlines()
            y = [float(line.split(",")[0]) for line in lines[1:]]
            x = [float(line.split(",")[1]) for line in lines[1:]]

        plt.plot(x, y, label= str(mesh_size_tab[i]) + ' nodes')
        i=i+1

    end = time.time()

    # Plot over diagonal line
    plt.legend()
    plt.xlabel('Position on slice circle')
    plt.ylabel('Value on slice circle')
    plt.title('Plot over slice circle for finite elements \n for Laplace operator on 3D sphere meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_Slice.png")
    
    # Plot of iteration number
    plt.close()
    plt.plot(mesh_size_tab, iterations_tab, label='Number of iterations of the linear solver')
    plt.legend()
    plt.xlabel('Number of nodes')
    plt.ylabel('Number of iterations')
    plt.title('Number of CG iterations for finite elements \n for Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_IterationNumber.png")
    
    # Plot of linear solver residual
    plt.close()
    plt.plot(mesh_size_tab, residual_tab, label='Residual of the linear solver')
    plt.legend()
    plt.xlabel('Number of nodes')
    plt.ylabel('-Log(residual)')
    plt.title('CG residual for finite elements \n for Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_LinearSolverResidual.png")
    
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(diameter_tab,diameter_tab)
    a2=np.sum(diameter_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,diameter_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation3DSphereEF() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print( "FE on 3D sphere triangle mesh : scheme order is ", -a)
    assert abs(a-0.775)<0.1

    # Plot of convergence curves
    plt.close()
    plt.plot(diameter_tab, error_tab, label='Log(|numerical-exact|)')
    plt.plot(diameter_tab, a*np.array(diameter_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('Log(h)')
    plt.ylabel('Log(error)')
    plt.title('Convergence of finite elements for \n Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='Log(cpu time))')
    plt.legend()
    plt.xlabel('Log(sqrt(number of nodes)')
    plt.ylabel('Log(cpu time)')
    plt.title('Computational time of finite elements \n for Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dimension"]=3
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_cell_type"]="3DTriangles"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_Poisson_32DFV_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    test_validation3DSphereEF()
