/*
 * DoubleTabTests.hxx
 *
 *  Created on: 18 mars 2013
 *      Author: mekkas
 */

#ifndef MATRIXTESTS_HXX_
#define MATRIXTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Matrix.hxx"

class MatrixTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassMatrix( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class DoubleTab Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<MatrixTests>("testClassMatrix", &MatrixTests::testClassMatrix));
        return suiteOfTests;
    }
};

#endif /* MATRIXTESTS_HXX_ */
