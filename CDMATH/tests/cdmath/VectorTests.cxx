/*
 * VectorTests.cxx
 *
 *  Created on: 18 mars 2013
 *      Author: mekkas
 */


#include "VectorTests.hxx"
#include <cmath>

using namespace std;

//----------------------------------------------------------------------
void
VectorTests::testClassVector( void )
//----------------------------------------------------------------------
{
    Vector A(2);
    A(0)=1.;
    A(1)=2.;
	CPPUNIT_ASSERT_EQUAL( 1.0, A(0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A(1) );
	CPPUNIT_ASSERT_EQUAL( sqrt(5.), A.norm() );

	Vector B;
	B=A;
	CPPUNIT_ASSERT_EQUAL( 1.0, B(0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, B(1) );

	Vector C(A+B);
	CPPUNIT_ASSERT_EQUAL( 2.0, C(0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, C(1) );

	double val=A*C;
	CPPUNIT_ASSERT_EQUAL( 10.0, val );

	Vector D(A-B);
	CPPUNIT_ASSERT_EQUAL( 0.0, D(0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, D(1) );

	Vector E;
	E=2*A;
	CPPUNIT_ASSERT_EQUAL( 2.0, E(0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, E(1) );

	E/=2;
	CPPUNIT_ASSERT_EQUAL( 1.0, E(0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, E(1) );

	E=A*2;
	CPPUNIT_ASSERT_EQUAL( 2.0, E(0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, E(1) );

	Vector F;
	F=A/2;
	CPPUNIT_ASSERT_EQUAL( A(0)/2, F(0) );
	CPPUNIT_ASSERT_EQUAL( A(1)/2, F(1) );

	F*=2;
	CPPUNIT_ASSERT_EQUAL( A(0), F(0) );
	CPPUNIT_ASSERT_EQUAL( A(1), F(1) );
	double a=A(0);
	CPPUNIT_ASSERT_EQUAL( A(0), a );

	Vector v3(4);
	v3(0)=1;
	v3(1)=2;
	v3(2)=3;
	v3(3)=4;

	Vector v4(3);
	v4(0)=1.;
	v4(1)=2.;
	v4(2)=3.;

	Matrix v5=v3^v4;

	CPPUNIT_ASSERT_EQUAL( 1., v5(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2., v5(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3., v5(0,2) );
	CPPUNIT_ASSERT_EQUAL( 2., v5(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4., v5(1,1) );
	CPPUNIT_ASSERT_EQUAL( 6., v5(1,2) );
	CPPUNIT_ASSERT_EQUAL( 3., v5(2,0) );
	CPPUNIT_ASSERT_EQUAL( 6., v5(2,1) );
	CPPUNIT_ASSERT_EQUAL( 9., v5(2,2) );
	CPPUNIT_ASSERT_EQUAL( 4., v5(3,0) );
	CPPUNIT_ASSERT_EQUAL( 8., v5(3,1) );
	CPPUNIT_ASSERT_EQUAL( 12., v5(3,2) );
}

