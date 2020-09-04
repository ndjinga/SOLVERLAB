/*
 * nodetests.cxx
 *
 *  Created on: 02 fevrier. 2012
 *      Authors: CDMAT
 */

#include "NodeTests.hxx"

//----------------------------------------------------------------------
void
NodeTests::testClassNode( void )
//----------------------------------------------------------------------
{
    Point P(0.5,0.5,0.0);
	Node n1(4,4,3,P);
	Node n(n1);
	CPPUNIT_ASSERT_EQUAL( 4, n.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( 4, n.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 3, n.getNumberOfEdges() );
	CPPUNIT_ASSERT_EQUAL( 0.5, n.getPoint().x() );
	CPPUNIT_ASSERT_EQUAL( 0.5, n.getPoint().y() );
	CPPUNIT_ASSERT_EQUAL( 0.0, n.getPoint().z() );
	CPPUNIT_ASSERT_EQUAL( 0.5, n.x() );
	CPPUNIT_ASSERT_EQUAL( 0.5, n.y() );
	CPPUNIT_ASSERT_EQUAL( 0.0, n.z() );
	Node n2;
    n2=n1;
    n2.addFaceId(0,10);
    n2.addFaceId(1,11);
    n2.addFaceId(2,12);
    n2.addFaceId(3,13);

    CPPUNIT_ASSERT_EQUAL( 10, n2.getFacesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 11, n2.getFacesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 12, n2.getFacesId()[2] );
	CPPUNIT_ASSERT_EQUAL( 13, n2.getFacesId()[3] );

    n2.addCellId(0,20);
    n2.addCellId(1,21);
    n2.addCellId(2,22);
    n2.addCellId(3,23);

    CPPUNIT_ASSERT_EQUAL( 20, n2.getCellsId()[0] );
	CPPUNIT_ASSERT_EQUAL( 21, n2.getCellsId()[1] );
	CPPUNIT_ASSERT_EQUAL( 22, n2.getCellsId()[2] );
	CPPUNIT_ASSERT_EQUAL( 23, n2.getCellsId()[3] );

    n2.addNeighbourNodeId(0,5);
    n2.addNeighbourNodeId(1,6);
    n2.addNeighbourNodeId(2,7);

    CPPUNIT_ASSERT_EQUAL( 5, n2.getNeighbourNodesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 6, n2.getNeighbourNodesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 7, n2.getNeighbourNodesId()[2] );

    n2=n;
    CPPUNIT_ASSERT_EQUAL( 0., n.distance(n2) );
}
