# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF avec quadrature aux noeuds de l'équation de Poisson 2D -\triangle u = f avec f raide
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2025
# Description : Utilisation de la quadrature standard du premier ordre pour la discrétisation de f sur un triangle ABC :
#                 \int f = 1/3 ( f(A)+f(B)+f(C) ) * |ABC|
#               Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage triangulaire
#               Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte u = tanh(-alpha * (y - 0.5 - 0.25 * sin(2 * pi * x)))+tanh( alpha * (y - x))
#================================================================================================================================

from math import tanh, sin, cos, pi
import CoreFlows as cf
import cdmath

alpha = 20. #Stiffness coefficient for the right hand side function

def equation_terms(x, y):
    A = tanh(-alpha * (y - 0.5 - 0.25 * sin(2 * pi * x) ) )
    B = tanh( alpha * (y - x) )

    dAx  =  alpha/2 * pi * cos(2 * pi * x) * (1 - A**2)
    dAxx =  alpha * pi * (pi * sin(2 * pi * x) * (A**2 - 1) - cos(2 * pi * x) * dAx * A)
    dAy  =  alpha * (A**2 - 1)
    dAyy =2*alpha * dAy * A

    dBx = alpha * (B**2 - 1)
    dBxx = 2*alpha * dBx * B
    dBy = alpha * (1 - B**2)
    dByy = -2*alpha * dBy * B

    fm = dAxx + dAyy + dBxx + dByy

    um = A + B

    return ( -um, -fm )
    
def StationaryDiffusionEquation_2DEF_StiffRHS_NODE_quadrature(mesh_file):
    
    mesh = cdmath.Mesh(str(mesh_file))
    boundaryIds = mesh.getBoundaryNodeIds()
    dirichletValues = {}

    #Remplissage des conditions limites
    for i in boundaryIds:
        Ni = mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()
        (u, f) = equation_terms(x, y)
        dirichletValues[i] = u

    dim = 2
    useFEdiscretisation = True
    problem = cf.StationaryDiffusionEquation(dim, useFEdiscretisation)
    problem.setMesh(mesh)
    problem.setDirichletValues(cf.MapIntDouble(dirichletValues))

    #Remplissage du terme source
    sourceField = cdmath.Field("sourceField", cdmath.NODES, mesh, 1)
    exactSol    = cdmath.Field("exactSol"   , cdmath.NODES, mesh, 1)
    for i in range(mesh.getNumberOfNodes()):
        Ni = mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()

        (u, f) = equation_terms(x, y)
        sourceField[i] = -f
        exactSol[i] = u
    
    problem.setHeatPowerField(sourceField)
    problem.setFileName(mesh_file)
    problem.setLinearSolver(cf.CG, cf.CHOLESKY)
    problem.setPrecision(1e-6)
    problem.initialize()
    problem.solveStationaryProblem()
    
    simSolution = problem.getOutputTemperatureField()
    simSolution.writeVTK("simSolution")
    exactSol.writeVTK("exactSol")
    sourceField.writeVTK("sourceField")
    
    #Calcul de l'erreur commise par rapport à la solution exacte + vérification du principe du maximum
    #==================================================================================================
    max_sol_num=exactSol.max()
    min_sol_num=exactSol.min()
    max_abs_sol_exacte=max(max_sol_num,-min_sol_num)
    erreur_abs= (exactSol - simSolution).normMax()[0]

    print("Absolute error = max(| exact solution - numerical solution |)                         = ",erreur_abs )
    print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
    print("\nSource field : minimum value= ", sourceField.min(),"maximum value= ", sourceField.max() )
    print("Numerical solution : minimum value= ", simSolution.min(),"maximum value= ", simSolution.max() )
    print("Exact solution     : minimum value= ", exactSol.min(),   "maximum value= ", exactSol.max() )

StationaryDiffusionEquation_2DEF_StiffRHS_NODE_quadrature('../resources/squareWithTriangles.med')
