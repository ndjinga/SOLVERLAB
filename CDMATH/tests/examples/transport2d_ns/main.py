#!/usr/bin/env python
# -*-coding:utf-8 -*

import math

import cdmath


def initial_conditions(my_mesh):
    rayon = 0.15
    xcentre = 0.25
    ycentre = 0.25
    y_field = cdmath.Field("Y field", cdmath.CELLS, my_mesh, 1)
    nbCells = my_mesh.getNumberOfCells()
    for j in range(nbCells):
        x = my_mesh.getCell(j).x()
        y = my_mesh.getCell(j).y()
        valX = (x - xcentre) * (x - xcentre)
        valY = (y - ycentre) * (y - ycentre)
        val = math.sqrt(valX + valY)
        if val < rayon:
            y_field[j] = 1.0
            pass
        else:
            y_field[j] = 0.0
            pass
        pass
    return y_field


def sigma_flux(VitesseX, VitesseY, cfl, y_field, indexFacesPerio):
    # Calculation of fluxes #
    SumFlux = cdmath.Field("Fluxes", cdmath.CELLS, y_field.getMesh(), 1)
    my_mesh = y_field.getMesh()
    nbCells = my_mesh.getNumberOfCells()
    normU = math.sqrt(VitesseX * VitesseX + VitesseY * VitesseY)
    for j in range(nbCells):
        Cj = my_mesh.getCell(j)
        nbFace = Cj.getNumberOfFaces()
        SumF = 0.0
        minlengthFk = 1.E30
        for k in range(nbFace):
            indexFace = Cj.getFacesId()[k]
            Fk = my_mesh.getFace(indexFace)
            NormalX = Cj.getNormalVector(k, 0)
            NormalY = Cj.getNormalVector(k, 1)
            LengthFk = Fk.getMeasure()
            UN = VitesseX * NormalX + VitesseY * NormalY
            minlengthFk = min(minlengthFk, LengthFk / abs(UN))
            minlengthFk = min(minlengthFk, LengthFk / abs(VitesseX))
            minlengthFk = min(minlengthFk, LengthFk / abs(VitesseY))
            conc = 0.0
            cellCourante = j
            cellAutre = -1
            if (not Fk.isBorder()):
                indexC1 = Fk.getCellsId()[0]
                indexC2 = Fk.getCellsId()[1]
                # hypothesis: the cell of index indexC1 is the current cell of index j #
                if (indexC1 == j):
                    # hypothese is verified #
                    cellCourante = indexC1
                    cellAutre = indexC2
                    pass
                elif (indexC2 == j):
                    # hypothesis is not verified #
                    cellCourante = indexC2
                    cellAutre = indexC1
                    pass
                # define left and right cell with the product of velocity * outward normal vector
                # if u*n>0: nothing to do, else invert left and right
                if (UN > 1.E-15):
                    conc = y_field[cellCourante]
                    pass
                else:
                    conc = y_field[cellAutre]
                    pass
                pass
            else:
                # homogeneous Neumann boundary conditions #
                if (Fk.getGroupName() == "GAUCHE" or Fk.getGroupName() == "DROITE"):
                    if (UN > 1.E-15):
                        conc = y_field[cellCourante]
                        pass
                    else:
                        conc = 0.0
                        pass
                    pass
                # periodical boundary conditions #
                if (Fk.getGroupName() == "BAS" or Fk.getGroupName() == "HAUT"):
                    indexFP = indexFacesPerio[indexFace]
                    # another way to get the index of the periodical face #
                    # int indexFP=my_mesh.getIndexFacePeriodic(indexFace);
                    Fp = my_mesh.getFace(indexFP)
                    indexCp = Fp.getCellsId()[0]
                    if (UN > 1.E-15):
                        conc = y_field[cellCourante]
                        pass
                    else:
                        conc = y_field[indexCp]
                        pass
                    pass
                pass
            SumF = SumF + UN * LengthFk * conc
            pass
        dt = cfl * minlengthFk / normU
        SumFlux[j] = dt * SumF / Cj.getMeasure()
        pass
    return dt, SumFlux


def EquationTransport2D(tmax, VitesseX, VitesseY, cfl, freqSortie, my_mesh, output_filename):

    # Initial conditions #
    print("Construction of the initial condition …")
    y_field = initial_conditions(my_mesh)
    #
    # MED output of the initial condition at t=0 and iter = 0
    #

    it = 0
    time = 0.
    print("Saving the solution at T=" + str(time) + "…")
    y_field.setTime(time, it)
    y_field.writeMED(output_filename)
    y_field.writeVTK(output_filename)
    y_field.writeCSV(output_filename)

    # Time loop #
    print("Resolution of the transport equation with an UPWIND scheme…")
    ntmax = 3
    indexFacesPerio = my_mesh.getIndexFacePeriodic()
    dt = 0.
    while (it < ntmax and time <= tmax):
        dt, SumFlux = sigma_flux(VitesseX, VitesseY, cfl, y_field, indexFacesPerio)
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))

        # Advancing one time step #
        y_field -= SumFlux
        time += dt
        it += 1
        # Output every freq times
        if (it % freqSortie == 0):
            y_field.setTime(time, it)
            y_field.writeMED(output_filename, False)
            y_field.writeVTK(output_filename, False)
            y_field.writeCSV(output_filename)
            pass
        pass
    return


def main():
    print("Resolution of the 2D transport equation:")
    print("- DOMAIN: SQUARE")
    print("- MESH: TRIANGULAR, GENERATED WITH SALOME")
    print("- PERIODICAL BC UP AND DOWN")
    print("- HOMOGENEOUS NEUMANN BC LEFT AND RIGHT")

    # Problem data
    cfl = 0.4
    VitesseX = 1.0
    VitesseY = 1.0
    tmax = 1.
    freqSortie = 10

    print("Loading triangular mesh …")
    my_mesh = cdmath.Mesh("../../tests/ressources/meshSquare.med")
    output_filename = "Exercice2PyTest"
    EquationTransport2D(tmax, VitesseX, VitesseY, cfl, freqSortie, my_mesh, output_filename)
    print("CDMATH calculation done.")
    return


if __name__ == """__main__""":
    main()
