import EulerIsothermal1DSchemeComparison
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json


def test_validationEulerIsothermal_1D_mesh_refinement():
    start = time.time()
    meshList=[50,100,150,200]
    meshType="Regular_1D grid"
    mesh_name='regular1DGrid'
    testColor="Green"
    nbMeshes=len(meshList)
    var_tot_upwind=[0]*nbMeshes
    var_tot_centered=[0]*nbMeshes
    var_tot_staggered=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    #Space Domain of the Simulations is the interval [a,b]
    a=0
    b=2
    cfl=1
    c0=300
    plt.close('all')
    i=0
    # Storing of total variations
    for nx in meshList:
        ntmax=int(nx/4)
        var_tot_upwind[i], var_tot_staggered[i], var_tot_centered[i] =EulerIsothermal1DSchemeComparison.solve( a,b,nx,"1D","RegularGrid",cfl,c0,ntmax)
        i=i+1
    
    end = time.time()

    # Plot of total variation curves
    plt.close()
    plt.figure()
    plt.plot(meshList, var_tot_upwind, label='Implicit upwind')
    plt.plot(meshList, var_tot_staggered, label='Implicit (pseudo)staggered')
    plt.plot(meshList, var_tot_centered, label='Implicit centered')
    plt.xlabel('Mesh size')
    plt.ylabel('Total variation')
    plt.title('Total variation of finite volume schemes')
    plt.legend()
    plt.savefig(mesh_name+"_EulerIsothermalTotalVariationMeshRefinement.png")
    
    
if __name__ == """__main__""":
    test_validationEulerIsothermal_1D_mesh_refinement()
