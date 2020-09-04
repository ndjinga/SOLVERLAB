/*
 * celltests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#ifndef CELLTESTS_HXX_
#define CELLTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Point.hxx"
#include "Cell.hxx"

class CellTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassCell( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Cell Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<CellTests>("testClassCell", &CellTests::testClassCell));
        return suiteOfTests;
    }
};
#endif /* CELLTESTS_HXX_ */
