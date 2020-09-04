import cdmath
import os

# Function to use in case of a single field to discretize

def discretizeFunction(func, mesh,var_name):
    myField=cdmath.Field(var, cdmath.CELLS, mesh,1)
    nbCells=mesh.getNumberOfCells()
    spaceDim=mesh.getSpaceDimension()

    if spaceDim==1 :
        for i in range(nbCells):
            x=mesh.getCell(i).x()
            myField[i]=func(x)
    elif spaceDim==2 :
        for i in range(nbCells):
            x=mesh.getCell(i).x()
            y=mesh.getCell(i).y()
            myField[i]=func(x,y)
    elif spaceDim==3 :
        for i in range(nbCells):
            x=mesh.getCell(i).x()
            y=mesh.getCell(i).y()
            z=mesh.getCell(i).z()
            myField[i]=func(x,y,z)
    else :
        raise ValueError("Wrong space dimension : expected dimension 1, 2 or 3")

    return myField

# Class to use in case of many fields to discretize (requires the creation of a dictionary)

class analyticalFunction_discretizer:

    def __init__(self, mesh, var_dict):
        self.mesh = mesh
        self.var_dict = var_dict
        self.var_list = self.var_dict.keys()
        self.fields = {var_name: self.eval(
            var_name) for var_name in self.var_list}

# useless
    # def ShowAllVar(self):
    #     for varName in self.varList:
    #         print varName


    def eval(self, var_name):
        func = self.var_dict[var_name]
        myField=cdmath.Field(var_name, cdmath.CELLS, self.mesh,1)
        nbCells=self.mesh.getNumberOfCells()
        spaceDim=self.mesh.getSpaceDimension()

        if spaceDim==1 :
            for i in range(nbCells):
                x=self.mesh.getCell(i).x()
                myField[i]=func(x)
        elif spaceDim==2 :
            for i in range(nbCells):
                x=self.mesh.getCell(i).x()
                y=self.mesh.getCell(i).y()
                myField[i]=func(x,y)
        elif spaceDim==3 :
            for i in range(nbCells):
                x=self.mesh.getCell(i).x()
                y=self.mesh.getCell(i).y()
                z=self.mesh.getCell(i).z()
                myField[i]=func(x,y,z)
        else :
            raise ValueError("Wrong space dimension : expected dimension 1, 2 or 3")

        return myField

    def get_variable(self, var):
        return self.fields[var]

    def save_all_variablesVTK(self):
        for var in self.var_list:
            self.get_variable(var).writeVTK(var)

    def save_all_variablesMED(self):
        for var in self.var_list:
            self.get_variable(var).writeMED(var)

    def save_all_variablesCSV(self):
        for var in self.var_list:
            self.get_variable(var).writeCSV(var)
