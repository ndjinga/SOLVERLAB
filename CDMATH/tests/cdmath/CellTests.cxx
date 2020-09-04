/*
 * celltests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#include "CellTests.hxx"

//----------------------------------------------------------------------
void
CellTests::testClassCell( void )
//----------------------------------------------------------------------
{
    Point P(0.5,0.5,0.0);
	Cell c1(4,4,1.0,P);
    c1.addNormalVector(0,0.2,0.3,0.0);
	CPPUNIT_ASSERT_EQUAL( 0.2, c1.getNormalVector(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.3, c1.getNormalVector(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, c1.getNormalVector(0,2) );
	Cell c(c1);
	CPPUNIT_ASSERT_EQUAL( 1.0, c.getMeasure() );
	CPPUNIT_ASSERT_EQUAL( 4, c.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 4, c.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 0.5, c.getBarryCenter().x() );
	CPPUNIT_ASSERT_EQUAL( 0.5, c.getBarryCenter().y() );
	CPPUNIT_ASSERT_EQUAL( 0.0, c.getBarryCenter().z() );
	CPPUNIT_ASSERT_EQUAL( 0.5, c.x() );
	CPPUNIT_ASSERT_EQUAL( 0.5, c.y() );
	CPPUNIT_ASSERT_EQUAL( 0.0, c.z() );
	Cell c2;
    c2=c1;
    c2.addNormalVector(1,0.4,0.6,0.0);
	CPPUNIT_ASSERT_EQUAL( 0.2, c2.getNormalVector(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.3, c2.getNormalVector(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, c2.getNormalVector(0,2) );
	CPPUNIT_ASSERT_EQUAL( 0.4, c2.getNormalVector(1,0) );
	CPPUNIT_ASSERT_EQUAL( 0.6, c2.getNormalVector(1,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, c2.getNormalVector(1,2) );

	c2=c1;
    c2.addFaceId(0,10);
    c2.addFaceId(1,11);
    c2.addFaceId(2,12);
    c2.addFaceId(3,13);
    c2.addNodeId(0,20);
    c2.addNodeId(1,21);
    c2.addNodeId(2,22);
    c2.addNodeId(3,23);

    CPPUNIT_ASSERT_EQUAL( 10, c2.getFacesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 11, c2.getFacesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 12, c2.getFacesId()[2] );
	CPPUNIT_ASSERT_EQUAL( 13, c2.getFacesId()[3] );
    CPPUNIT_ASSERT_EQUAL( 20, c2.getNodesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 21, c2.getNodesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 22, c2.getNodesId()[2] );
	CPPUNIT_ASSERT_EQUAL( 23, c2.getNodesId()[3] );

    CPPUNIT_ASSERT_EQUAL( 10, c2.getFacesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 11, c2.getFacesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 12, c2.getFacesId()[2] );
	CPPUNIT_ASSERT_EQUAL( 13, c2.getFacesId()[3] );
    CPPUNIT_ASSERT_EQUAL( 20, c2.getNodesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 21, c2.getNodesId()[1] );
	CPPUNIT_ASSERT_EQUAL( 22, c2.getNodesId()[2] );
	CPPUNIT_ASSERT_EQUAL( 23, c2.getNodesId()[3] );
}
