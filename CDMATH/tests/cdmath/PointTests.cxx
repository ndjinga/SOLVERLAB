/*
 * pointtests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#include "PointTests.hxx"
#include <cmath>

//----------------------------------------------------------------------
void
PointTests::testClassPoint( void )
//----------------------------------------------------------------------
{
   Point P1(1., 2., 3.);
   CPPUNIT_ASSERT_EQUAL( 1., P1.x() );
   CPPUNIT_ASSERT_EQUAL( 1., P1[0] );
   CPPUNIT_ASSERT_EQUAL( 2., P1.y() );
   CPPUNIT_ASSERT_EQUAL( 2., P1[1] );
   CPPUNIT_ASSERT_EQUAL( 3., P1.z() );
   CPPUNIT_ASSERT_EQUAL( 3., P1[2] );

   Point P2(1., 2., 3.);
   CPPUNIT_ASSERT_EQUAL( 14., P1.dot(P2) );

   Point P3;
   P3=P1+P2;
   CPPUNIT_ASSERT_EQUAL( 2., P3.x());
   CPPUNIT_ASSERT_EQUAL( 4., P3.y());
   CPPUNIT_ASSERT_EQUAL( 6., P3.z());

   Point P5(1., 2., 3.);
   Point P6(3., 5., 0.);
   Point P4=P5-P6;
   CPPUNIT_ASSERT_EQUAL( -2., P4.x());
   CPPUNIT_ASSERT_EQUAL( -3., P4.y());
   CPPUNIT_ASSERT_EQUAL( 3., P4.z());

   P5+=P6;
   CPPUNIT_ASSERT_EQUAL( 4., P5.x());
   CPPUNIT_ASSERT_EQUAL( 7., P5.y());
   CPPUNIT_ASSERT_EQUAL( 3., P5.z());

   Point P7;
   P7[0]=1.0;
   P7[1]=2.0;
   P7[2]=3.0;
   Point P8(3., 5., 0.);
   P7-=P8;
   CPPUNIT_ASSERT_EQUAL( -2., P7.x());
   CPPUNIT_ASSERT_EQUAL( -3., P7.y());
   CPPUNIT_ASSERT_EQUAL( 3., P7.z());

   Point P9;
   P9=P1*3.;
   CPPUNIT_ASSERT_EQUAL( 3., P9.x());
   CPPUNIT_ASSERT_EQUAL( 6., P9.y());
   CPPUNIT_ASSERT_EQUAL( 9., P9.z());

   Point P10(1., 2., 3.);
   P10*=3.0;
   CPPUNIT_ASSERT_EQUAL( 3., P10.x());
   CPPUNIT_ASSERT_EQUAL( 6., P10.y());
   CPPUNIT_ASSERT_EQUAL( 9., P10.z());

   double norm=P1.norm();
   CPPUNIT_ASSERT_EQUAL( sqrt(14.), norm);

   Point P11(1., 2., 3.);
   Point P12(4., 5., 6.);
   double dx=P12.x()-P11.x();
   double dy=P12.y()-P11.y();
   double dz=P12.z()-P11.z();
   double distance=sqrt(dx*dx+dy*dy+dz*dz);
   CPPUNIT_ASSERT_EQUAL( distance, P11.distance(P12));

   Point P13(3., 6., 9.);
   Point P14;
   P14=P13/3.0;
   CPPUNIT_ASSERT_EQUAL( 1., P14.x());
   CPPUNIT_ASSERT_EQUAL( 2., P14.y());
   CPPUNIT_ASSERT_EQUAL( 3., P14.z());

   Point P15(3., 6., 9.);
   P15/=3.0;
   CPPUNIT_ASSERT_EQUAL( 1., P15.x());
   CPPUNIT_ASSERT_EQUAL( 2., P15.y());
   CPPUNIT_ASSERT_EQUAL( 3., P15.z());

   Point P16(3., 6., 9.);
   Point P17(P16);
   CPPUNIT_ASSERT_EQUAL( P16.x(), P17.x());
   CPPUNIT_ASSERT_EQUAL( P16.y(), P17.y());
   CPPUNIT_ASSERT_EQUAL( P16.z(), P17.z());
}
