/*
 * PointTests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#ifndef POINTTESTS_HXX_
#define POINTTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Point.hxx"

class PointTests : public CppUnit::TestCase 
{
    public: //----------------------------------------------------------------
	void testClassPoint( void );

	static CppUnit::Test *suite()
	{
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Point Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<PointTests>("testClassPoint",  &PointTests::testClassPoint));
        return suiteOfTests;
    }
};

#endif /* POINTTESTS_HXX_ */
