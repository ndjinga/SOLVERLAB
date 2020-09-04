/*
 * fieldtests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#ifndef FIELDTESTS_HXX_
#define FIELDTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Mesh.hxx"
#include "Field.hxx"

class FieldTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassField( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Field Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<FieldTests>("testClassField", &FieldTests::testClassField));
        return suiteOfTests;
    }
};
#endif /* FIELDTESTS_HXX_ */
