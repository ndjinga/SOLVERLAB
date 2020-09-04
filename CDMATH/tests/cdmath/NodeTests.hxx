/*
 * vertextests.hxx
 *
 *  Created on: 02 fevrier. 2012
 *      Authors: CDMAT
 */

#ifndef NODETESTS_HXX_
#define NODETESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Node.hxx"

class NodeTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassNode( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Node Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<NodeTests>("testClassNode",   &NodeTests::testClassNode));
        return suiteOfTests;
    }
};
#endif /* NODETESTS_HXX_ */
