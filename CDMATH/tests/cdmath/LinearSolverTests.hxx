/*
 * celltests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#ifndef LINEARSOLVERTESTS_HXX_
#define LINEARSOLVERTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "LinearSolver.hxx"

class LinearSolverTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassLinearSolver( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class LinearSolver Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<LinearSolverTests>("testClassLinearSolver", &LinearSolverTests::testClassLinearSolver));
        return suiteOfTests;
    }
};
#endif /* LINEARSOLVERTESTS_HXX_ */
