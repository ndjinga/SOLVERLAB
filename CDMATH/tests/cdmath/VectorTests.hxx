/*
 * DoubleTabTests.hxx
 *
 *  Created on: 18 mars 2013
 *      Author: mekkas
 */

#ifndef VECTORTESTS_HXX_
#define VECTORTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Vector.hxx"

class VectorTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassVector( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class DoubleTab Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<VectorTests>("testClassVector", &VectorTests::testClassVector));
        return suiteOfTests;
    }
};

#endif /* VECTORTESTS_HXX_ */
