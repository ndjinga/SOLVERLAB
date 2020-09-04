/*
 * facetests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMATH
 */

#include "FaceTests.hxx"

//----------------------------------------------------------------------
void
FaceTests::testClassFace( void )
//----------------------------------------------------------------------
{
	Point p(0,1,2);
	Face f1(2,2,1.0,p,1.,2.,3.);
	Face f(f1);
	CPPUNIT_ASSERT_EQUAL( 1.0, f.getMeasure() );
	CPPUNIT_ASSERT_EQUAL( 2, f.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 2, f.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( p.x(), f.getBarryCenter().x() );
	CPPUNIT_ASSERT_EQUAL( p.y(), f.getBarryCenter().y() );
	CPPUNIT_ASSERT_EQUAL( p.z(), f.getBarryCenter().z() );
	CPPUNIT_ASSERT_EQUAL( p.x(), f.x() );
	CPPUNIT_ASSERT_EQUAL( p.y(), f.y() );
	CPPUNIT_ASSERT_EQUAL( p.z(), f.z() );
	CPPUNIT_ASSERT_EQUAL( 1., f.getXN() );
	CPPUNIT_ASSERT_EQUAL( 2., f.getYN() );
	CPPUNIT_ASSERT_EQUAL( 3., f.getZN() );

	Face f2;
	f2=f1;
    f2.addCellId(0,10);
    f2.addCellId(1,11);
    f2.addNodeId(0,20);
    f2.addNodeId(1,21);

    CPPUNIT_ASSERT_EQUAL( 10, f2.getCellsId()[0] );
	CPPUNIT_ASSERT_EQUAL( 11, f2.getCellsId()[1] );
    CPPUNIT_ASSERT_EQUAL( 20, f2.getNodesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 21, f2.getNodesId()[1] );

    CPPUNIT_ASSERT_EQUAL( 10, f2.getCellsId()[0] );
	CPPUNIT_ASSERT_EQUAL( 11, f2.getCellsId()[1] );
    CPPUNIT_ASSERT_EQUAL( 20, f2.getNodesId()[0] );
	CPPUNIT_ASSERT_EQUAL( 21, f2.getNodesId()[1] );

	f2=f;
	CPPUNIT_ASSERT_EQUAL( 1.0, f2.getMeasure() );
	CPPUNIT_ASSERT_EQUAL( 2, f2.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 2, f2.getNumberOfCells() );
}
