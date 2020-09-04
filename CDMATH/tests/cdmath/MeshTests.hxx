/*
 * meshtests.hxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#ifndef MESHTESTS_HXX_
#define MESHTESTS_HXX_

#include <string>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>

#include "Mesh.hxx"

class MeshTests : public CppUnit::TestCase
{
    public: //----------------------------------------------------------------
      void testNormals( Mesh mesh );
      void testClassMesh( void );
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Sample Unit Class Mesh Tests");
        suiteOfTests->addTest(new CppUnit::TestCaller<MeshTests>("testClassMesh", &MeshTests::testClassMesh));
        return suiteOfTests;
    }
};
#endif /* MESHTESTS_HXX_ */
