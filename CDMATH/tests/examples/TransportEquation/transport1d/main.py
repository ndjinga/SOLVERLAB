#!/usr/bin/env python
# -*-coding:utf-8 -*

import math

import cdmath


def main():
    a = -5.0
    b = 5.0
    nx = 1000
    ntmax = 1000
    dx = (b - a) / nx
    pi = 3.1415927
    # Transport velocity
    cfl = 0.5
    u = 3.
    dt = cfl * dx / u

    my_mesh = cdmath.Mesh(a, b, nx)
    conc = cdmath.Field("Concentration", cdmath.CELLS, my_mesh, 1)

    # Initial conditions
    sigma = math.sqrt(0.2)
    for i in range(my_mesh.getNumberOfCells()):
        x = my_mesh.getCell(i).x()
        conc[i] = 0.5 / (sigma * math.sqrt(2 * pi)) * math.exp(-0.5 * math.pow((x / sigma), 2))
        pass

    time = 0.
    tmax = 3.0
    it = 0

    print("MED post-treatment of the solution at T=" + str(time) + "…")
    output_filename = "EqTr1D"
    conc.setTime(time, it)
    conc.writeMED(output_filename)
    conc.writeVTK(output_filename)
    conc.writeCSV(output_filename)
    output_freq = 10

    # Time loop
    while (it < ntmax and time <= tmax):
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
        conc[0] = conc[0] - u * dt / dx * (conc[0] - conc[my_mesh.getNumberOfCells() - 1])
        for j in range(1, my_mesh.getNumberOfCells()):
            conc[j] = conc[j] - u * dt / dx * (conc[j] - conc[j - 1])
            pass
        time += dt
        it += 1
        if (it % output_freq == 0):
            conc.setTime(time, it)
            conc.writeMED(output_filename, False)
            conc.writeVTK(output_filename, False)
            conc.writeCSV(output_filename)
            pass
        pass
    print("CDMATH calculation done.")
    return


if __name__ == """__main__""":
    main()
