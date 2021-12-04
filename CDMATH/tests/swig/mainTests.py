# coding: latin-1

import math
import unittest

from cdmath import *


class TestsCDMATHSwig(unittest.TestCase):

    def testClassPoint(self):
        P1 = Point(1., 2., 3.)
        self.assertTrue(P1.x() == 1.)
        self.assertTrue(P1.y() == 2.)
        self.assertTrue(P1.z() == 3.)
        self.assertTrue(P1[0] == 1.)
        self.assertTrue(P1[1] == 2.)
        self.assertTrue(P1[2] == 3.)
        P2 = Point(1., 2., 3.)
        self.assertTrue(14., P1.dot(P2))

        P3 = P1 + P2
        self.assertTrue(2., P3.x())
        self.assertTrue(4., P3.y())
        self.assertTrue(6., P3.z())

        P5 = Point(1., 2., 3.)
        P6 = Point(3., 5., 0.)
        P4 = P5 - P6
        self.assertTrue(-2., P4.x())
        self.assertTrue(-3., P4.y())
        self.assertTrue(3., P4.z())

        P5 += P6
        self.assertTrue(4., P5.x())
        self.assertTrue(7., P5.y())
        self.assertTrue(3., P5.z())

        P7 = Point()
        P7[0] = 1.0
        P7[1] = 2.0
        P7[2] = 3.0
        P8 = Point(3., 5., 0.)
        P7 -= P8
        self.assertTrue(-2., P7.x())
        self.assertTrue(-3., P7.y())
        self.assertTrue(3., P7.z())

        P9 = Point()
        P9 = P1 * 3.0
        self.assertTrue(3., P9.x())
        self.assertTrue(6., P9.y())
        self.assertTrue(9., P9.z())

        P10 = Point(1., 2., 3.)
        P10 *= 3.0
        self.assertTrue(3., P10.x())
        self.assertTrue(6., P10.y())
        self.assertTrue(9., P10.z())

        norm = P1.norm()
        self.assertTrue(math.sqrt(14.), norm)

        P11 = Point(1., 2., 3.)
        P12 = Point(4., 5., 6.)
        dx = P12.x() - P11.x()
        dy = P12.y() - P11.y()
        dz = P12.z() - P11.z()
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
        self.assertTrue(distance, P11.distance(P12))

        P13 = Point(3., 6., 9.)
        P14 = Point()
        P14 = P13 / 3.0
        self.assertTrue(1., P14.x())
        self.assertTrue(2., P14.y())
        self.assertTrue(3., P14.z())

        P15 = Point(3., 6., 9.)
        P15 /= 3.0
        self.assertTrue(1., P15.x())
        self.assertTrue(2., P15.y())
        self.assertTrue(3., P15.z())
        return

    def testClassField(self):
        M = Mesh(0.0, 1.0, 10, 0., 1., 5)

        conc1 = Field("CONCENTRATION", CELLS, M, 2, 1.2)
        self.assertTrue(1.2 == conc1.getTime())
        for j in range(conc1.getNumberOfComponents()):
            for i in range(conc1.getNumberOfElements()):
                conc1[i, j] = i + j

        conc1n = Field("CONCENTRATION", NODES, M, 2, 1.2)
        self.assertTrue(1.2 == conc1n.getTime())
        for j in range(conc1n.getNumberOfComponents()):
            for i in range(conc1n.getNumberOfElements()):
                conc1n[i, j] = i

        conc3 = Field("CONCENTRATION", CELLS, M, 2, 1.2)
        for j in range(conc3.getNumberOfComponents()):
            for i in range(conc3.getNumberOfElements()):
                conc3[i, j] = -(i + j)

        v1 = conc3.getValuesOnComponent(1)
        v2 = conc3.getValuesOnAllComponents(4)

        for i in range(conc3.getNumberOfElements()):
            self.assertTrue(-(i + 1) == v1[i])

        for j in range(conc3.getNumberOfComponents()):
            self.assertTrue(-(4 + j) == v2[j])

        fileNameVTK = "champc"
        conc1.writeVTK(fileNameVTK)

        fileNameMED = "champc"
        conc1.writeMED(fileNameMED)
        conc1.setTime(2.3, 1)
        conc1.writeMED(fileNameMED, False)
        for i in range(conc1.getNumberOfElements()):
            self.assertTrue(1.0 * i == conc1[i])
            pass

        self.assertTrue(2 == conc1.getNumberOfComponents())
        self.assertTrue(50 == conc1.getNumberOfElements())
        self.assertTrue(2.3 == conc1.getTime())

        fileNameVTK = "champn"
        conc1n.writeVTK(fileNameVTK)

        fileNameMED = "champn"
        conc1n.writeMED(fileNameMED)
        conc1n.setTime(2.3, 1)
        conc1n.writeMED(fileNameMED, False)
        for i in range(conc1n.getNumberOfElements()):
            self.assertTrue(1.0 * i == conc1n[i])
            pass

        self.assertTrue(2 == conc1n.getNumberOfComponents())
        self.assertTrue(66 == conc1n.getNumberOfElements())
        self.assertTrue(2.3 == conc1n.getTime())

        conc6 = Field("CONCENTRATION", NODES, M, 2)
        for i in range(conc6.getNumberOfComponents()):
            for j in range(conc6.getNumberOfElements()):
                conc6[j, i] = i * 1.0 + 2. * j
                pass
            pass
        for i in range(conc6.getNumberOfComponents()):
            for j in range(conc6.getNumberOfElements()):
                self.assertTrue(1.0 * i + 2. * j == conc6[j, i])
                self.assertTrue(1.0 * i + 2. * j == conc6.getValues()
                                [i + j * conc6.getNumberOfComponents()])
                pass
            pass

        conc6 = Field("CONCENTRATION", CELLS, M, 2)
        conc6.setInfoOnComponent(0, "compo1")
        conc6.setInfoOnComponent(1, "compo2")
        self.assertTrue(conc6.getInfoOnComponent(0) == "compo1")
        self.assertTrue(conc6.getInfoOnComponent(1) == "compo2")
        for i in range(conc6.getNumberOfComponents()):
            for j in range(conc6.getNumberOfElements()):
                conc6[j, i] = i * 1.0 + 2. * j
                pass
            pass
        for i in range(conc6.getNumberOfComponents()):
            for j in range(conc6.getNumberOfElements()):
                self.assertTrue(1.0 * i + 2. * j == conc6[j, i])
                self.assertTrue(1.0 * i + 2. * j == conc6.getValues()
                                [i + j * conc6.getNumberOfComponents()])
                pass
            pass

        self.assertTrue(2 == conc1.getNumberOfComponents())
        self.assertTrue(50 == conc1.getNumberOfElements())

        conc3 = conc1
        for i in range(conc3.getNumberOfElements()):
            conc3[i, 0] = i * 1.0
            pass

        x = conc3[2]
        self.assertTrue(x == 2.0)

        for i in range(conc3.getNumberOfElements()):
            self.assertTrue(1.0 * i == conc3[i])
            pass
        self.assertTrue(2 == conc3.getNumberOfComponents())
        self.assertTrue(50 == conc3.getNumberOfElements())

        conc6 = conc3 + conc1
        for i in range(conc6.getNumberOfElements()):
            self.assertTrue(2.0 * i == conc6[i])
            pass
        self.assertTrue(2 == conc6.getNumberOfComponents())
        self.assertTrue(50 == conc6.getNumberOfElements())

        conc6 = conc3 - conc1
        for i in range(conc6.getNumberOfElements()):
            self.assertTrue(0.0 == conc6[i])
            pass
        self.assertTrue(2 == conc6.getNumberOfComponents())
        self.assertTrue(50 == conc6.getNumberOfElements())

        conc6 = conc1
        conc6 += conc1
        for i in range(conc6.getNumberOfElements()):
            self.assertTrue(2.0 * i == conc6[i, 0])
            pass
        self.assertTrue(2 == conc6.getNumberOfComponents())
        self.assertTrue(50 == conc6.getNumberOfElements())

        for i in range(conc6.getNumberOfElements()):
            conc6[i, 0] = i * 1.0
            pass
        conc6 *= 2.0
        for i in range(conc6.getNumberOfElements()):
            self.assertTrue(2.0 * i == conc6[i])
        self.assertTrue(2 == conc6.getNumberOfComponents())
        self.assertTrue(50 == conc6.getNumberOfElements())

#       conc7=Field("CONCENTRATION",NODES,M,2) ;
#        conc7.setFieldByMEDCouplingFieldDouble(conc1n.getField());
#        conc7.setName("CONC")
#        self.assertTrue( conc7.getName() == "CONC" );
#        for i in range(conc7.getNumberOfElements()):
#            self.assertTrue( conc1n[i] == conc7[i] );
#            pass
#        self.assertTrue( 2 == conc7.getNumberOfComponents() );
#        self.assertTrue( 66 == conc7.getNumberOfElements() );

#        conc7=Field("CONCENTRATION",CELLS,M,2) ;
#        conc7.setFieldByMEDCouplingFieldDouble(conc1.getField());
#        conc7.setName("CONC")
#        self.assertTrue( conc7.getName() == "CONC" );
#        for i in range(conc7.getNumberOfElements()):
#            self.assertTrue( conc1[i] == conc7[i] );
#            pass
#        self.assertTrue( 2 == conc7.getNumberOfComponents() );
#        self.assertTrue( 50 == conc7.getNumberOfElements() );

        conc8 = Field("CONCENTRATION", CELLS, M)
        for i in range(conc8.getNumberOfElements()):
            conc8[i] = i * 1.0
            pass
        for i in range(conc8.getNumberOfElements()):
            self.assertTrue(1.0 * i == conc8[i])
            pass
        self.assertTrue(1 == conc8.getNumberOfComponents())
        self.assertTrue(50 == conc8.getNumberOfElements())

        conc8 = Field("CONCENTRATION", NODES, M)
        for i in range(conc8.getNumberOfElements()):
            conc8[i] = i * 1.0
            pass
        for i in range(conc8.getNumberOfElements()):
            self.assertTrue(1.0 * i == conc8[i])
            pass
        self.assertTrue(1 == conc8.getNumberOfComponents())
        self.assertTrue(66 == conc8.getNumberOfElements())

        conc9 = Field("CONCENTRATION", CELLS, M)
        for i in range(conc9.getNumberOfElements()):
            conc9[i] = i * 1.0
            pass
        conc9 /= 2.0
        for i in range(conc9.getNumberOfElements()):
            self.assertTrue(1.0 * i / 2. == conc9[i])
            pass

        conc10 = conc8
        for i in range(conc10.getNumberOfElements()):
            conc10[i] = i * 1.0
            pass
        conc10 -= 2.0
        for i in range(conc10.getNumberOfElements()):
            self.assertTrue(1.0 * i - 2.0 == conc10[i])

        conc11 = conc8
        for i in range(conc11.getNumberOfElements()):
            conc11[i] = i * 1.0
            pass
        conc11 += 2.0
        for i in range(conc11.getNumberOfElements()):
            self.assertTrue(1.0 * i + 2. == conc11[i])

        conc12 = conc8
        for i in range(conc12.getNumberOfElements()):
            conc12[i] = i * 1.0
            pass
        conc12 += conc8
        for i in range(conc12.getNumberOfElements()):
            self.assertTrue(2.0 * i == conc12[i])

        conc13 = conc8
        for i in range(conc13.getNumberOfElements()):
            conc13[i] = i * 1.0
            pass
        conc13 -= conc8
        for i in range(conc13.getNumberOfElements()):
            self.assertTrue(0.0 == conc13[i])

        conc15 = conc1 * 2.
        conc16 = conc1 * 0.3 

        for i in range(conc15.getNumberOfElements()):
            self.assertTrue(conc1[i] * 2. == conc15[i])
            self.assertTrue(conc1[i] * 0.3 == conc16[i])

        MF = Mesh(0.0, 1.0, 3, 0., 1., 3)
        concF1 = Field("CONCENTRATION", FACES, MF)
        for j in range(concF1.getNumberOfComponents()):
            for i in range(concF1.getNumberOfElements()):
                concF1[i, j] = i + j

        for j in range(concF1.getNumberOfComponents()):
            for i in range(concF1.getNumberOfElements()):
                self.assertTrue(i + j == concF1[i, j])

        self.assertTrue(1 == concF1.getNumberOfComponents())
        self.assertTrue(0.0 == concF1.getTime())
        self.assertTrue(24 == concF1.getNumberOfElements())
        return

    def testClassCell(self):
        P = Point(0.5, 0.5, 0.0)
        c1 = Cell(4, 4, 1.0, P)
        c1.addNormalVector(0, 0.2, 0.3, 0.0)
        self.assertTrue(0.2 == c1.getNormalVector(0, 0))
        self.assertTrue(0.3 == c1.getNormalVector(0, 1))
        self.assertTrue(0.2 == c1.getNormalVectors()[0])
        self.assertTrue(0.3 == c1.getNormalVectors()[1])
        c = Cell()
        c = c1
        self.assertTrue(1.0 == c.getMeasure())
        self.assertTrue(4 == c.getNumberOfNodes())
        self.assertTrue(4 == c.getNumberOfFaces())
        self.assertTrue(0.5 == c.getBarryCenter().x())
        self.assertTrue(0.5 == c.getBarryCenter().y())
        self.assertTrue(0.0 == c.getBarryCenter().z())
        self.assertTrue(0.5 == c.x())
        self.assertTrue(0.5 == c.y())
        self.assertTrue(0.0 == c.z())
        c2 = c1
        c2.addNormalVector(1, 0.4, 0.6, 0.0)
        self.assertTrue(0.2 == c2.getNormalVector(0, 0))
        self.assertTrue(0.3 == c2.getNormalVector(0, 1))
        self.assertTrue(0.4 == c2.getNormalVector(1, 0))
        self.assertTrue(0.6 == c2.getNormalVector(1, 1))

        c2 = c1
        c2.addFaceId(0, 10)
        c2.addFaceId(1, 11)
        c2.addFaceId(2, 12)
        c2.addFaceId(3, 13)
        c2.addNodeId(0, 20)
        c2.addNodeId(1, 21)
        c2.addNodeId(2, 22)
        c2.addNodeId(3, 23)

        self.assertTrue(10 == c2.getFacesId()[0])
        self.assertTrue(11 == c2.getFacesId()[1])
        self.assertTrue(12 == c2.getFacesId()[2])
        self.assertTrue(13 == c2.getFacesId()[3])
        self.assertTrue(20 == c2.getNodesId()[0])
        self.assertTrue(21 == c2.getNodesId()[1])
        self.assertTrue(22 == c2.getNodesId()[2])
        self.assertTrue(23 == c2.getNodesId()[3])
        return

    def testClassNode(self):
        P = Point(0.5, 0.5, 0.0)
        n1 = Node(4, 4, 3, P)
        n = Node()
        n = n1
        self.assertTrue(4 == n.getNumberOfCells())
        self.assertTrue(4, n.getNumberOfFaces())
        self.assertTrue(3, n.getNumberOfEdges())
        self.assertTrue(0.5 == n.getPoint().x())
        self.assertTrue(0.5 == n.getPoint().y())
        self.assertTrue(0.0 == n.getPoint().z())
        self.assertTrue(0.5 == n.x())
        self.assertTrue(0.5 == n.y())
        self.assertTrue(0.0 == n.z())
        n2 = n1
        n2.addFaceId(0, 10)
        n2.addFaceId(1, 11)
        n2.addFaceId(2, 12)
        n2.addFaceId(3, 13)
        n2.addCellId(0, 20)
        n2.addCellId(1, 21)
        n2.addCellId(2, 22)
        n2.addCellId(3, 23)

        self.assertTrue(10 == n2.getFacesId()[0])
        self.assertTrue(11 == n2.getFacesId()[1])
        self.assertTrue(12 == n2.getFacesId()[2])
        self.assertTrue(13 == n2.getFacesId()[3])
        self.assertTrue(20 == n2.getCellsId()[0])
        self.assertTrue(21 == n2.getCellsId()[1])
        self.assertTrue(22 == n2.getCellsId()[2])
        self.assertTrue(23 == n2.getCellsId()[3])
        n2 = n
        self.assertTrue(0. == n.distance(n2))
        return

    def testClassFace(self):
        p = Point(0, 1, 2)
        f1 = Face(2, 2, 1.0, p, 1., 2., 3.)
        f = Face()
        f = f1
        self.assertTrue(1.0 == f.getMeasure())
        self.assertTrue(2 == f.getNumberOfNodes())
        self.assertTrue(2 == f.getNumberOfCells())
        self.assertTrue(p.x() == f.getBarryCenter().x())
        self.assertTrue(p.y() == f.getBarryCenter().y())
        self.assertTrue(p.z() == f.getBarryCenter().z())
        self.assertTrue(p.x() == f.x())
        self.assertTrue(p.y() == f.y())
        self.assertTrue(p.z() == f.z())
        self.assertTrue(-1 == f.getRegion())
        self.assertTrue(False == f.isBorder())
        f.setGroupName("Bord1")
        self.assertTrue(0 == f.getRegion())
        self.assertTrue(True == f.isBorder())
        self.assertTrue("Bord1" == f.getGroupName())

        f2 = f1
        f2.addCellId(0, 10)
        f2.addCellId(1, 11)
        f2.addNodeId(0, 20)
        f2.addNodeId(1, 21)

        self.assertTrue(10 == f2.getCellsId()[0])
        self.assertTrue(11 == f2.getCellsId()[1])
        self.assertTrue(20 == f2.getNodesId()[0])
        self.assertTrue(21 == f2.getNodesId()[1])

        f2 = f
        self.assertTrue(1.0 == f2.getMeasure())
        self.assertTrue(2 == f2.getNumberOfNodes())
        self.assertTrue(2 == f2.getNumberOfCells())
        self.assertTrue(0 == f2.getRegion())
        self.assertTrue(True == f2.isBorder())
        self.assertTrue("Bord1" == f2.getGroupName())
        return

    def testClassMesh(self):
        import os
        M1 = Mesh(0.0, 4.0, 4)
        self.assertTrue(1 == M1.getSpaceDimension())
        self.assertTrue(5 == M1.getNumberOfNodes())
        self.assertTrue(4 == M1.getNumberOfCells())
        self.assertTrue(5 == M1.getNumberOfFaces())
        self.assertTrue(0. == M1.getFace(0).x())
        self.assertTrue(0. == M1.getNode(0).x())
        self.assertTrue(1. == M1.getFace(1).x())
        self.assertTrue(1. == M1.getNode(1).x())
        self.assertTrue(2. == M1.getFace(2).x())
        self.assertTrue(2. == M1.getNode(2).x())
        self.assertTrue(3. == M1.getFace(3).x())
        self.assertTrue(3. == M1.getNode(3).x())
        self.assertTrue(4. == M1.getFace(4).x())
        self.assertTrue(4. == M1.getNode(4).x())
        x11 = M1.getCell(1).x()
        y11 = M1.getCell(1).y()
        self.assertTrue(x11 == 1.5)
        self.assertTrue(y11 == 0.0)
        M1.setGroupAtFaceByCoords(0., 0., 0., 1.E-14, "LeftEdge")
        M1.setGroupAtFaceByCoords(4., 0., 0., 1.E-14, "RightEdge")
        M1.setGroupAtNodeByCoords(0., 0., 0., 1.E-14, "LeftEdge")
        M1.setGroupAtNodeByCoords(4., 0., 0., 1.E-14, "RightEdge")
        self.assertTrue(M1.getFace(0).isBorder() == True)
        self.assertTrue(M1.getFace(1).isBorder() == False)
        self.assertTrue(M1.getFace(2).isBorder() == False)
        self.assertTrue(M1.getFace(3).isBorder() == False)
        self.assertTrue(M1.getFace(4).isBorder() == True)
        self.assertTrue(len(M1.getNameOfFaceGroups()) ==3)
        self.assertTrue(M1.getNameOfFaceGroups()[2] == "RightEdge")
        self.assertTrue(M1.getNameOfFaceGroups()[1] == "LeftEdge")
        self.assertTrue(M1.getNameOfFaceGroups()[0] == "Boundary")
        self.assertTrue(M1.getNameOfNodeGroups()[2] == "RightEdge")
        self.assertTrue(M1.getNameOfNodeGroups()[1] == "LeftEdge")
        self.assertTrue(M1.getNameOfNodeGroups()[0] == "Boundary")

        M1 = Mesh(0.0, 1.0, 4)
        self.assertTrue(1 == M1.getSpaceDimension())
        self.assertTrue(5 == M1.getNumberOfNodes())
        self.assertTrue(4 == M1.getNumberOfCells())
        self.assertTrue(5 == M1.getNumberOfFaces())

        xinf = 0.0
        xsup = 4.0
        yinf = 0.0
        ysup = 4.0
        M2 = Mesh(xinf, xsup, 4, yinf, ysup, 4)

        self.assertEqual(4, M2.getNy())
        self.assertEqual(4, M2.getNx())
        self.assertEqual(2, M2.getSpaceDimension())
        self.assertEqual(25, M2.getNumberOfNodes())
        self.assertEqual(16, M2.getNumberOfCells())
        self.assertEqual(40, M2.getNumberOfFaces())
#        x1=M2.getCells()[4].x();
#        y1=M2.getCells()[4].y();
#        self.assertTrue( x1==0.5 );
#        self.assertTrue( y1==1.5 );

#        x2=M2.getNodes()[24].x();
#        y2=M2.getNodes()[24].y();
#        self.assertTrue( x2==4. );
#        self.assertTrue( y2==4. );
        eps = 1.E-10
        M2.setGroupAtPlan(xsup, 0, eps, "RightEdge")
        M2.setGroupAtPlan(xinf, 0, eps, "LeftEdge")
        M2.setGroupAtPlan(yinf, 1, eps, "BottomEdge")
        M2.setGroupAtPlan(ysup, 1, eps, "TopEdge")
        self.assertTrue(len(M2.getNameOfFaceGroups()) == 5)
        self.assertTrue(M2.getNameOfFaceGroups()[0] == "Boundary")
        self.assertTrue(M2.getNameOfFaceGroups()[1] == "RightEdge")
        self.assertTrue(M2.getNameOfFaceGroups()[2] == "LeftEdge")
        self.assertTrue(M2.getNameOfFaceGroups()[3] == "BottomEdge")
        self.assertTrue(M2.getNameOfFaceGroups()[4] == "TopEdge")
        nbFaces = M2.getNumberOfFaces()
        M2.setPeriodicFaces()
        indexFaces = M2.getIndexFacePeriodic()
        for i in range(nbFaces):
            #            x=M2.getFaces()[i].x();
            #            y=M2.getFaces()[i].y();
            x = M2.getFace(i).x()
            y = M2.getFace(i).y()
            if (abs(y) < 1.E-10 and abs(x - 0.5) < 1.E-10):
                indexFace = M2.getIndexFacePeriodic(i)
                xi = M2.getFace(indexFace).x()
                yi = M2.getFace(indexFace).y()
                self.assertTrue(xi == x)
                self.assertTrue(yi == ysup)
                self.assertTrue(True == M2.getFace(indexFace).isBorder())
                pass
                pass
            pass

        M2.writeMED("TestMesh")
        M22 = Mesh("TestMesh.med")
        self.assertTrue(2 == M22.getSpaceDimension())
        self.assertTrue(25 == M22.getNumberOfNodes())
        self.assertTrue(16 == M22.getNumberOfCells())
        self.assertTrue(40 == M22.getNumberOfFaces())

        M23 = Mesh("meshSquare.med")
        self.assertTrue(len(M23.getNameOfFaceGroups()) == 5)
        print( M23.getNameOfFaceGroups() )
        self.assertTrue(M23.getNameOfFaceGroups()[1] == "Bottom")
        self.assertTrue(M23.getNameOfFaceGroups()[2] == "Left")
        self.assertTrue(M23.getNameOfFaceGroups()[3] == "Right")
        self.assertTrue(M23.getNameOfFaceGroups()[4] == "Top")
        self.assertTrue(M23.getNameOfFaceGroups()[0] == "Boundary")

        M3 = M1
        self.assertTrue(1 == M3.getSpaceDimension())
        self.assertTrue(5 == M3.getNumberOfNodes())
        self.assertTrue(4 == M3.getNumberOfCells())
        self.assertTrue(5 == M3.getNumberOfFaces())

        M3 = M2
        self.assertTrue(2 == M3.getSpaceDimension())
        self.assertTrue(25 == M3.getNumberOfNodes())
        self.assertTrue(16 == M3.getNumberOfCells())
        self.assertTrue(40 == M3.getNumberOfFaces())

        M4 = M3
        self.assertTrue(2 == M4.getSpaceDimension())
        self.assertTrue(25 == M4.getNumberOfNodes())
        self.assertTrue(16 == M4.getNumberOfCells())
        self.assertTrue(40 == M4.getNumberOfFaces())

        M5 = Mesh(0.0, 1.0, 4, 0.0, 1.0, 4, 0.0, 1.0, 4)
        self.assertTrue(3 == M5.getSpaceDimension())
        fileNameVTK = "TestMesh"
        M4.writeVTK(fileNameVTK)
        fileNameMED = "TestMesh"
        M4.writeMED(fileNameMED)
        M6 = Mesh(fileNameMED + ".med")
        self.assertTrue(2 == M6.getSpaceDimension())
        self.assertTrue(25 == M6.getNumberOfNodes())
        self.assertTrue(16 == M6.getNumberOfCells())
        self.assertTrue(40 == M6.getNumberOfFaces())
        return

    def testClassMatrix(self):
        A = Matrix(2, 2)
        A[0, 0] = 1.
        A[0, 1] = 2.
        A[1, 0] = 3.
        A[1, 1] = 4.

        self.assertTrue(1.0, A[0, 0])
        self.assertTrue(2.0, A[0, 1])
        self.assertTrue(3.0, A[1, 0])
        self.assertTrue(4.0, A[1, 1])

        self.assertTrue(-2., A.determinant())

        A1 = Matrix(2, 2)
        A1 = A
        self.assertTrue(1.0, A1[0, 0])
        self.assertTrue(2.0, A1[0, 1])
        self.assertTrue(3.0, A1[1, 0])
        self.assertTrue(4.0, A1[1, 1])

        A11 = Matrix(2, 2)

        A11 = A1 - A

        self.assertTrue(0.0 == A11[0, 0])
        self.assertTrue(0.0 == A11[0, 1])
        self.assertTrue(0.0 == A11[1, 0])
        self.assertTrue(0.0 == A11[1, 1])

        A11 = A1 + A

        self.assertTrue(2.0 == A11[0, 0])
        self.assertTrue(4.0 == A11[0, 1])
        self.assertTrue(6.0 == A11[1, 0])
        self.assertTrue(8.0 == A11[1, 1])

        A22 = Matrix(2, 2)
        A22 = A * 2

        self.assertTrue(2.0 == A22[0, 0])
        self.assertTrue(4.0 == A22[0, 1])
        self.assertTrue(6.0 == A22[1, 0])
        self.assertTrue(8.0 == A22[1, 1])

        A22 = A * 3
        self.assertTrue(3.0 == A22[0, 0])
        self.assertTrue(6.0 == A22[0, 1])
        self.assertTrue(9.0 == A22[1, 0])
        self.assertTrue(12.0 == A22[1, 1])

        A22 = A *0.5
        self.assertTrue(0.5 == A22[0, 0])
        self.assertTrue(1.0 == A22[0, 1])
        self.assertTrue(1.5 == A22[1, 0])
        self.assertTrue(2.0 == A22[1, 1])

        A2 = Matrix(A1)
        A2 *= 2
        self.assertTrue(2.0 == A2[0, 0])
        self.assertTrue(4.0 == A2[0, 1])
        self.assertTrue(6.0 == A2[1, 0])
        self.assertTrue(8.0 == A2[1, 1])

        A2 /= 2
        self.assertTrue(1.0 == A2[0, 0])
        self.assertTrue(2.0 == A2[0, 1])
        self.assertTrue(3.0 == A2[1, 0])
        self.assertTrue(4.0 == A2[1, 1])

        A2 -= A
        self.assertTrue(0.0 == A2[0, 0])
        self.assertTrue(0.0 == A2[0, 1])
        self.assertTrue(0.0 == A2[1, 0])
        self.assertTrue(0.0 == A2[1, 1])

        A2 += A
        self.assertTrue(1.0 == A2[0, 0])
        self.assertTrue(2.0 == A2[0, 1])
        self.assertTrue(3.0 == A2[1, 0])
        self.assertTrue(4.0 == A2[1, 1])

        X = Vector(2)
        X[0] = 1.
        X[1] = 2.

        X1 = A * X
        self.assertTrue(5. == X1[0])
        self.assertTrue(11.0 == X1[1])

        self.assertTrue(True == A2.isSquare())
        self.assertTrue(False == A2.isSymmetric())

        A3 = Matrix(2, 3)
        A3[0, 0] = 1.
        A3[0, 1] = 2.
        A3[0, 2] = 2.
        A3[1, 0] = 3.
        A3[1, 1] = 4.
        A3[1, 2] = 4.
        self.assertTrue(False == A3.isSquare())

        A[0, 0] = 1.
        A[0, 1] = -2.
        A[1, 0] = -2.
        A[1, 1] = 4.
        self.assertTrue(True == A.isSymmetric())

        A4 = Matrix(4, 4)
        A4[0, 0] = 1.
        A4[0, 1] = 2.
        A4[0, 2] = 3.
        A4[0, 3] = 4.
        A4[1, 0] = 5.
        A4[1, 1] = 6.
        A4[1, 2] = 7.
        A4[1, 3] = 8.
        A4[2, 0] = 9.
        A4[2, 1] = 10.
        A4[2, 2] = 11.
        A4[2, 3] = 12.
        A4[3, 0] = 13.
        A4[3, 1] = 14.
        A4[3, 2] = 15.
        A4[3, 3] = 16.
        A5 = Matrix(A4.transpose())
        self.assertTrue(1. == A5[0, 0])
        self.assertTrue(5. == A5[0, 1])
        self.assertTrue(9. == A5[0, 2])
        self.assertTrue(13. == A5[0, 3])
        self.assertTrue(2. == A5[1, 0])
        self.assertTrue(6. == A5[1, 1])
        self.assertTrue(10. == A5[1, 2])
        self.assertTrue(14. == A5[1, 3])
        self.assertTrue(3. == A5[2, 0])
        self.assertTrue(7. == A5[2, 1])
        self.assertTrue(11. == A5[2, 2])
        self.assertTrue(15. == A5[2, 3])
        self.assertTrue(4. == A5[3, 0])
        self.assertTrue(8. == A5[3, 1])
        self.assertTrue(12. == A5[3, 2])
        self.assertTrue(16. == A5[3, 3])

        self.assertTrue(0. == A5.determinant())

    def testClassVector(self):
        A = Vector(2)
        A[0] = 1.
        A[1] = 2.
        self.assertTrue(1.0 == A[0])
        self.assertTrue(2.0 == A[1])
        self.assertTrue(math.sqrt(5.) == A.norm())

        B = A
        self.assertTrue(1.0 == B[0])
        self.assertTrue(2.0 == B[1])

        C = A + B
        self.assertTrue(2.0 == C[0])
        self.assertTrue(4.0 == C[1])

        val = A * C
        self.assertTrue(10.0 == val)

        D = A - B
        self.assertTrue(0.0 == D[0])
        self.assertTrue(0.0 == D[1])

        E = A * 2
        self.assertTrue(2.0 == E[0])
        self.assertTrue(4.0 == E[1])

        E *= 0.5
        self.assertTrue(1.0 == E[0])
        self.assertTrue(2.0 == E[1])

        E = A * 2
        self.assertTrue(2.0, E[0])
        self.assertTrue(4.0, E[1])

        F = A * 0.5
        self.assertTrue(A[0] / 2 == F[0])
        self.assertTrue(A[1] / 2 == F[1])

        F *= 2
        self.assertTrue(A[0] == F[0])
        self.assertTrue(A[1] == F[1])
        a = A[0]
        self.assertTrue(A[0] == a)

#        v3=Vector(4);
#        v3[0]=1;
#        v3[1]=2;
#        v3[2]=3;
#        v3[3]=4;

#        v4=Vector(3);
#        v4[0]=1.;
#        v4[1]=2.;
#        v4[2]=3.;

#        v5=v3^v4;
#
#         self.assertTrue( 1.==v5[0,0] );
#         self.assertTrue( 2.==v5[0,1] );
#         self.assertTrue( 3.==v5[0,2] );
#         self.assertTrue( 2.==v5[1,0] );
#         self.assertTrue( 4.==v5[1,1] );
#         self.assertTrue( 6.==v5[1,2] );
#         self.assertTrue( 3.==v5[2,0] );
#         self.assertTrue( 6.==v5[2,1] );
#         self.assertTrue( 9.==v5[2,2] );
#         self.assertTrue( 4.==v5[3,0] );
#         self.assertTrue( 8.==v5[3,1] );
#         self.assertTrue( 12.==v5[3,2] );

if __name__ == '__main__':
    unittest.main()
