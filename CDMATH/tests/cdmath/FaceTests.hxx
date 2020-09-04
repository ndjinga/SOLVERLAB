/*
 * facetests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMATH
 */

#ifndef FACETESTS_HXX_
#define FACETESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Point.hxx"
#include "Face.hxx"

class FaceTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassFace( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Face Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<FaceTests>("testClassFace", &FaceTests::testClassFace));
        return suiteOfTests;
    }
};
#endif /* FACETESTS_HXX_ */
