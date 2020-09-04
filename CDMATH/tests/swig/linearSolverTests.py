#!/usr/bin/env python
# -*-coding:utf-8 -*-

import unittest

from cdmath import *


class TestsLinearSolverSwig(unittest.TestCase):

    def testClassLinearSolver(self):
        A = Matrix(2, 2)
        A[0, 0] = 3.
        A[0, 1] = -2.
        A[1, 0] = -2.
        A[1, 1] = 4.

        A *= A.transpose()

        self.assertEqual(A[0,0], 13)
        self.assertEqual(A[0,1],-14)
        self.assertEqual(A[1,0],-14)
        self.assertEqual(A[1,1], 20)

        Xana = Vector(2)
        Xana[0] = 1.
        Xana[1] = 2.

        B = A * Xana

        self.assertEqual(B[0],-15)
        self.assertEqual(B[1], 26)

        LS=LinearSolver(A,B,500,1.E-10,"GMRES","LU");
        X=LS.solve();
        self.assertTrue(abs(X[0]-Xana[0])<1.E-10);
        self.assertTrue(abs(X[1]-Xana[1])<1.E-10);
        self.assertEqual(LS.getStatus(),True);

        self.assertEqual(LS.getNumberMaxOfIter(),500);
        self.assertEqual(LS.getTolerance(),1.E-10);
        self.assertEqual(LS.getNameOfMethod(),"GMRES");
        self.assertEqual(LS.getNumberOfIter(),1);
        self.assertEqual(LS.isMatrixSingular(),False);
        self.assertEqual(LS.getNameOfPc(),"LU");

        LS2 = LinearSolver(A, B, 500, 1.E-10, "CG")
        A1 = SparseMatrixPetsc(2, 2,4)
        A1[0, 0] =  1.
        A1[0, 1] = -1.
        A1[1, 0] = -1.
        A1[1, 1] =  1.

        B1 = Vector(2)
        B1[0] = 2.
        B1[1] =-2.

        LS2.setMatrix(A1 * -1.)
        # LS2.setSndMember(B1 * -1)
        # LS2.setTolerance(1.E-10)
        # LS2.setNumberMaxOfIter(10)
        # LS2.setMatrixIsSingular(True)
        # X2 = LS2.solve()
        # self.assertTrue(abs(X2[0] -   1) < 1.E-10)
        # self.assertTrue(abs(X2[1] - (-1)) < 1.E-10)
        # self.assertEqual(LS2.getStatus(), True)
        # self.assertEqual(LS2.getNumberOfIter(), 1)
        # self.assertEqual(LS2.isMatrixSingular(), True)
        # self.assertEqual(LS2.getNameOfMethod(), "CG")
# #
        # LS3 = LinearSolver(A, B, 500, 1.E-10, "BICG")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "BICG")
# #
        # LS3 = LinearSolver(A, B, 500, 1.E-10, "CR")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "CR")
# #
        # LS3 = LinearSolver(A, B, 500, 1.E-10, "CGS")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "CGS")

        # LS3 = LinearSolver(A, B, 500, 1.E-10, "GMRES")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "GMRES")

        # LS3=LinearSolver(A,B,500,1.E-10,"BICG","LU");
        # X3=LS3.solve();
        # self.assertTrue(abs(X3[0]-Xana[0])<1.E-10);
        # self.assertTrue(abs(X3[1]-Xana[1])<1.E-10);
        # self.assertEqual(LS3.getStatus(),True);
        # self.assertEqual(LS3.getNameOfMethod(),"BICG");
        # self.assertEqual(LS3.getNameOfPc(),"LU");

        # LS3 = LinearSolver(A, B, 500, 1.E-10, "BICG")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "BICG")
        # self.assertEqual(LS3.getNameOfPc(),"");

        # LS3 = LinearSolver(A, B, 500, 1.E-10, "GCR")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "GCR")

        # LS3 = LinearSolver(A, B, 500, 1.E-10, "LSQR")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "LSQR")
# #
        # LS3 = LinearSolver(A, B, 500, 1.E-10, "CHOLESKY")
        # X3 = LS3.solve()
        # self.assertTrue(abs(X3[0] - Xana[0]) < 1.E-10)
        # self.assertTrue(abs(X3[1] - Xana[1]) < 1.E-10)
        # self.assertEqual(LS3.getStatus(), True)
        # self.assertEqual(LS3.getNameOfMethod(), "CHOLESKY")

        # LS3=LinearSolver(A,B,500,1.E-10,"LU");
        # X3=LS3.solve();
        # self.assertTrue(abs(X3[0]-Xana[0])<1.E-10);
        # self.assertTrue(abs(X3[1]-Xana[1])<1.E-10);
        # self.assertEqual(LS3.getStatus(),True);
        # self.assertEqual(LS3.getNameOfMethod(),"LU");
        # self.assertEqual(LS3.getNameOfPc(),"");

        A2 = SparseMatrixPetsc(6, 6, 16)
        A2[0, 0] = 2.
        A2[0, 1] = -1.

        A2[1, 0] = -1.
        A2[1, 1] = 2.
        A2[1, 2] = -1.

        A2[2, 1] = -1.
        A2[2, 2] = 2.
        A2[2, 3] = -1.

        A2[3, 2] = -1.
        A2[3, 3] = 2.
        A2[3, 4] = -1.

        A2[4, 3] = -1.
        A2[4, 4] = 2.
        A2[4, 5] = -1.

        A2[5, 4] = -1.
        A2[5, 5] = 2.

        Xana2 = Vector(6)
        Xana2[0] = 1.
        Xana2[1] = 2.
        Xana2[2] = 3.
        Xana2[3] = 4.
        Xana2[4] = 5.
        Xana2[5] = 6.

        B2 = A2 * Xana2

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "GMRES", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "CG", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "LGMRES", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "BICG", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "CR", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "CGS", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "BICG", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "GCR", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 1)

        LS11 = LinearSolver(A2, B2, 500, 1.E-10, "LSQR", "ILU")
        X11 = LS11.solve()
        for i in range(X11.getNumberOfRows()):
            self.assertTrue(abs(X11[i] - Xana2[i]) < 1.E-10)
            pass

        self.assertEqual(LS11.getStatus(), True)

        self.assertEqual(LS11.getNumberMaxOfIter(), 500)
        self.assertEqual(LS11.getTolerance(), 1.E-10)
        self.assertEqual(LS11.getNumberOfIter(), 6)


if __name__ == """__main__""":
    unittest.main()
