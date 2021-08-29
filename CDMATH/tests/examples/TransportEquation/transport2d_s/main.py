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
    # Fluxes calculation #
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
                # hypothese: La cellule d'index indexC1 est la cellule courante index j #
                if (indexC1 == j):
                    # hypothese verifie #
                    cellCourante = indexC1
                    cellAutre = indexC2
                    pass
                elif (indexC2 == j):
                    # hypothese non verifie #
                    cellCourante = indexC2
                    cellAutre = indexC1
                    pass
                # definir la cellule gauche et droite par le prduit vitesse * normale sortante
                # si u*n>0 : rien a faire sinon inverser la gauche et la droite
                if (UN > 1.E-15):
                    conc = y_field[cellCourante]
                    pass
                else:
                    conc = y_field[cellAutre]
                    pass
                pass
            else:
                # conditions aux limites neumann homogene #
                if (Fk.getGroupName() == "LeftEdge" or Fk.getGroupName() == "RightEdge"):
                    if (UN > 1.E-15):
                        conc = y_field[cellCourante]
                        pass
                    else:
                        conc = 0.0
                        pass
                    pass
                # conditions aux limites periodiques #
                if (Fk.getGroupName() == "BottomEdge" or Fk.getGroupName() == "TopEdge"):
                    indexFP = indexFacesPerio[indexFace]
                    # une autre manière de recuperer l'index de la face periodique #
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


def EquationTransport2D(tmax, VitesseX, VitesseY, cfl, output_freq, my_mesh, file):

    # Initial conditions #
    print("Construction of the initial condition …")
    y_field = initial_conditions(my_mesh)
    #
    #  MED output of the initial conditions at t=0 and iteration = 0
    #

    iteration = 0
    time = 0.
    print("Saving the solution at T=" + str(time) + "…")
    y_field.setTime(time, iteration)
    y_field.writeMED(file)
    y_field.writeVTK(file)
    y_field.writeCSV(file)

    # Time loop #
    print("Resolution of the transport equation with an UPWIND scheme …")
    ntmax = 3
    indexFacesPerio = my_mesh.getIndexFacePeriodic()
    dt = 0.
    while (iteration < ntmax and time <= tmax):
        dt, SumFlux = sigma_flux(VitesseX, VitesseY, cfl, y_field, indexFacesPerio)
        print("-- Iter: " + str(iteration) + ", Time: " + str(time) + ", dt: " + str(dt))

        # Advancing time step #
        y_field -= SumFlux
        time += dt
        iteration += 1
        # Output every output_freq iterations
        if (iteration % output_freq == 0):
            y_field.setTime(time, iteration)
            y_field.writeMED(file, False)
            y_field.writeVTK(file, False)
            y_field.writeCSV(file)
            pass
        pass
    return


def main():
    print("RESOLUTION OF THE 2D TRANSPORT EQUATION:")
    print("- DOMAIN: SQUARE [0,1]x[0,1]")
    print("- MESH: CARTESIAN, INTERNAL GENERATION WITH CDMATH")
    print("- PERIODIC BC ON TOP AND BOTTOM")
    print("- HOMOGENEOUS NEUMANN BC ON LEFT AND RIGHT")

    # Problem data
    cfl = 0.4
    VitesseX = 1.0
    VitesseY = 1.0
    tmax = 1.
    output_freq = 10

    print("Construction of a cartesian mesh …")
    xinf = 0.0
    xsup = 1.0
    yinf = 0.0
    ysup = 1.0
    nx = 100
    ny = 100
    my_mesh = cdmath.Mesh(xinf, xsup, nx, yinf, ysup, ny)
    eps = 1.E-10
    my_mesh.setGroupAtPlan(xsup, 0, eps, "RightEdge")
    my_mesh.setGroupAtPlan(xinf, 0, eps, "LeftEdge")
    my_mesh.setGroupAtPlan(yinf, 1, eps, "BottomEdge")
    my_mesh.setGroupAtPlan(ysup, 1, eps, "TopEdge")
    fileOutPutCart = "Exercice1PyTest"
    EquationTransport2D(tmax, VitesseX, VitesseY, cfl, output_freq, my_mesh, fileOutPutCart)
    print("CDMATH calculation done.")
    return


if __name__ == """__main__""":
    main()
