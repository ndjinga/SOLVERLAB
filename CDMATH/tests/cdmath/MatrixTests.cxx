/*
 * DoubleTabTests.cxx
 *
 *  Created on: 18 mars 2013
 *      Author: mekkas
 */



/*
 * celltests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMATH
 */

#include "MatrixTests.hxx"

#include "Vector.hxx"

using namespace std;

//----------------------------------------------------------------------
void
MatrixTests::testClassMatrix( void )
//----------------------------------------------------------------------
{
    Matrix A(2,2);
    A(0,0)=1.;
    A(0,1)=2.;
    A(1,0)=3.;
    A(1,1)=4.;

	CPPUNIT_ASSERT_EQUAL( 1.0, A(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A(1,1) );

	CPPUNIT_ASSERT_EQUAL( -2., A.determinant() );

    Matrix A1(2,2);
    A1=A;
	CPPUNIT_ASSERT_EQUAL( 1.0, A1(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A1(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A1(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A1(1,1) );

	Matrix A11(2,2);

	A11 = A1-A;

	CPPUNIT_ASSERT_EQUAL( 0.0, A11(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(1,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(1,1) );

	A11 = A1+A;

	CPPUNIT_ASSERT_EQUAL( 2.0, A11(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A11(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A11(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A11(1,1) );

	Matrix A22(2,2);
	A22 = 2*A;

	CPPUNIT_ASSERT_EQUAL( 2.0, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A22(1,1) );

	A22 = A*3;
	CPPUNIT_ASSERT_EQUAL( 3.0, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 9.0, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 12.0, A22(1,1) );

	A22 = A/2;
	CPPUNIT_ASSERT_EQUAL( 0.5, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 1.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 1.5, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A22(1,1) );

    Matrix A2(A1);
    A2*=2;
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A2(1,1) );

    A2/=2;
	CPPUNIT_ASSERT_EQUAL( 1.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(1,1) );

    A2-=A;
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(1,1) );

    A2+=A;
	CPPUNIT_ASSERT_EQUAL( 1.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(1,1) );

	Vector X(2);
    X(0)=1.;
    X(1)=2.;

    Vector X1(X);
    X1=A*X;
	CPPUNIT_ASSERT_EQUAL( 5., X1(0) );
	CPPUNIT_ASSERT_EQUAL( 11.0, X1(1) );

	CPPUNIT_ASSERT_EQUAL( true, A2.isSquare() );
	CPPUNIT_ASSERT_EQUAL( false, A2.isSymmetric() );

	Matrix A3(2,3);
    A3(0,0)=1.;
    A3(0,1)=2.;
    A3(0,2)=2.;
    A3(1,0)=3.;
    A3(1,1)=4.;
    A3(1,2)=4.;
	CPPUNIT_ASSERT_EQUAL( false, A3.isSquare() );

    A(0,0)=1.;
    A(0,1)=-2.;
    A(1,0)=-2.;
    A(1,1)=4.;
	CPPUNIT_ASSERT_EQUAL( true, A.isSymmetric() );

	Matrix A4(4,4);
    A4(0,0)=1.;
    A4(0,1)=2.;
    A4(0,2)=3.;
    A4(0,3)=4.;
    A4(1,0)=5.;
    A4(1,1)=6.;
    A4(1,2)=7.;
    A4(1,3)=8.;
    A4(2,0)=9.;
    A4(2,1)=10.;
    A4(2,2)=11.;
    A4(2,3)=12.;
    A4(3,0)=13.;
    A4(3,1)=14.;
    A4(3,2)=15.;
    A4(3,3)=16.;
    Matrix A5(A4.transpose());
	CPPUNIT_ASSERT_EQUAL( 1., A5(0,0) );
	CPPUNIT_ASSERT_EQUAL( 5., A5(0,1) );
	CPPUNIT_ASSERT_EQUAL( 9., A5(0,2) );
	CPPUNIT_ASSERT_EQUAL( 13., A5(0,3) );
	CPPUNIT_ASSERT_EQUAL( 2., A5(1,0) );
	CPPUNIT_ASSERT_EQUAL( 6., A5(1,1) );
	CPPUNIT_ASSERT_EQUAL( 10., A5(1,2) );
	CPPUNIT_ASSERT_EQUAL( 14., A5(1,3) );
	CPPUNIT_ASSERT_EQUAL( 3., A5(2,0) );
	CPPUNIT_ASSERT_EQUAL( 7., A5(2,1) );
	CPPUNIT_ASSERT_EQUAL( 11., A5(2,2) );
	CPPUNIT_ASSERT_EQUAL( 15., A5(2,3) );
	CPPUNIT_ASSERT_EQUAL( 4., A5(3,0) );
	CPPUNIT_ASSERT_EQUAL( 8., A5(3,1) );
	CPPUNIT_ASSERT_EQUAL( 12., A5(3,2) );
	CPPUNIT_ASSERT_EQUAL( 16., A5(3,3) );

	CPPUNIT_ASSERT_EQUAL( 0., A5.determinant() );
}

