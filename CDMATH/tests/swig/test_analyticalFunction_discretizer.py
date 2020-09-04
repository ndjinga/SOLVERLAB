#Testing the discretization of analytical functions
import cdmath
import analyticalFunction_discretizer as anSol


def Density(x):
    return x * 2


def Pressure(x):
    return x * 3

var_dict = {
    'Density': Density,
    'Pressure': Pressure,
}

xMin=0, xMax=1.0, nx=100
mesh=cdmath.Mesh(xMin, xMax, nx)
solution = anSol.analyticalFunction_discretizer(mesh, output_dir="./tmp",
                                      var_dict=var_dict)


print(solution.var_list)
for varName in solution.var_dict:
    print( varName + "=", solution.var_eval[varName])

solution.save_all_variables()

print( "done" )
