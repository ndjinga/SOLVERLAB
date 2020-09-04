/*
 * DoubleTabTests.hxx
 *
 *  Created on: 18 mars 2013
 *      Author: mekkas
 */

#ifndef SparseMatrixPetscTESTS_HXX_
#define SparseMatrixPetscTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "SparseMatrixPetsc.hxx"

class SparseMatrixPetscTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testClassSparseMatrixPetsc( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class DoubleTab Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<SparseMatrixPetscTests>("testClassSparseMatrixPetsc", &SparseMatrixPetscTests::testClassSparseMatrixPetsc));
        return suiteOfTests;
    }
};

#endif /* SparseMatrixPetscTESTS_HXX_ */
