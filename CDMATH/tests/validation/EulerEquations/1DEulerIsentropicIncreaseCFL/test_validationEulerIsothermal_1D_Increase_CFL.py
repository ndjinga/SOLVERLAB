import EulerIsothermal1DSchemeComparison
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json


def test_validationEulerIsothermal_1D_cfl_increase():
    start = time.time()
    cflList=[ .5, 1, 10, 20]
    meshType="Regular_1D grid"
    mesh_name='regular1DGrid'
    testColor="Green"
    nbSimulations=len(cflList)
    var_tot_upwind   =[0]*nbSimulations
    var_tot_centered =[0]*nbSimulations
    var_tot_staggered=[0]*nbSimulations
    mesh_size_tab    =[0]*nbSimulations
    #Space Domain of the Simulations is the interval [a,b]
    a=0
    b=2
    c0=300
    nx=400
    plt.close('all')
    i=0
    # Storing of total variations
    for cfl in cflList:
        ntmax=int(nx/(4*cfl))
        var_tot_upwind[i], var_tot_staggered[i], var_tot_centered[i] =EulerIsothermal1DSchemeComparison.solve( a,b,nx,"1D","RegularGrid",cfl,c0,ntmax)
        i=i+1
    
    end = time.time()

    # Plot of total variation curves
    plt.close()
    plt.figure()
    plt.plot(cflList, var_tot_upwind,    label='Implicit upwind')
    plt.plot(cflList, var_tot_staggered, label='Implicit staggered')
    plt.plot(cflList, var_tot_centered,  label='Implicit centered')
    plt.xlabel('CFLnumber')
    plt.ylabel('Total variation')
    plt.title('Total variation of finite volume schemes')
    plt.legend()
    plt.savefig(mesh_name+"_EulerIsothermalTotalVariationCFLIncrease.png")
    
    
if __name__ == """__main__""":
    test_validationEulerIsothermal_1D_cfl_increase()
