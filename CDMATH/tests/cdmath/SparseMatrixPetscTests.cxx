/*
 * SparseMatrixPetsctests.cxx
 *
 *  Created on: Dec. 2017
 *      Authors: CDMATH
 */

#include "SparseMatrixPetscTests.hxx"

#include "Vector.hxx"

using namespace std;

//----------------------------------------------------------------------
void
SparseMatrixPetscTests::testClassSparseMatrixPetsc( void )
//----------------------------------------------------------------------
{
    SparseMatrixPetsc A(2,2);
    A.setValue(0,0,1.);
    A.setValue(0,1,2.);
    A.setValue(1,0,3.);
    A.setValue(1,1,4.);

	CPPUNIT_ASSERT_EQUAL( 1.0, A(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A(1,1) );

	CPPUNIT_ASSERT_EQUAL( false, A.isSymmetric(1.e-5) );

    SparseMatrixPetsc A1(2,2);
    A1=A;
	CPPUNIT_ASSERT_EQUAL( 1.0, A1(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A1(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A1(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A1(1,1) );

	SparseMatrixPetsc A11(2,2);

	A11 = A1-A;

	CPPUNIT_ASSERT_EQUAL( 0.0, A11(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(1,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A11(1,1) );

	A11 = A1+A;

	CPPUNIT_ASSERT_EQUAL( 2.0, A11(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A11(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A11(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A11(1,1) );

	SparseMatrixPetsc A22(2,2);

	A22 = 2*A;
	CPPUNIT_ASSERT_EQUAL( 2.0, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A22(1,1) );

	A22 = A*3;
	CPPUNIT_ASSERT_EQUAL( 3.0, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 9.0, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 12.0, A22(1,1) );

	A22 = A/2;
	CPPUNIT_ASSERT_EQUAL( 0.5, A22(0,0) );
	CPPUNIT_ASSERT_EQUAL( 1.0, A22(0,1) );
	CPPUNIT_ASSERT_EQUAL( 1.5, A22(1,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A22(1,1) );

	SparseMatrixPetsc A2(A1);
    A2*=2;
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A2(1,1) );

    A2/=2;
	CPPUNIT_ASSERT_EQUAL( 1.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(1,1) );

    A2-=A;
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 0.0, A2(1,1) );

    A2+=A;
	CPPUNIT_ASSERT_EQUAL( 1.0, A2(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A2(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A2(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A2(1,1) );

	Vector X(2);
    X(0)=1.;
    X(1)=2.;

    Vector X1(X);
    X1=A*X;
	CPPUNIT_ASSERT_EQUAL( 5., X1(0) );
	CPPUNIT_ASSERT_EQUAL( 11.0, X1(1) );

	CPPUNIT_ASSERT_EQUAL( true, A2.isSquare() );
	CPPUNIT_ASSERT_EQUAL( false, A2.isSymmetric() );

	SparseMatrixPetsc A3(2,3);
    A3.setValue(0,0,1.);
    A3.setValue(0,1,2.);
    A3.setValue(0,2,2.);
    A3.setValue(1,0,3.);
    A3.setValue(1,1,4.);
    A3.setValue(1,2,4.);
	CPPUNIT_ASSERT_EQUAL( false, A3.isSquare() );

    A.setValue(0,0,1.);
    A.setValue(0,1,-1.);
    A.setValue(1,0,-1.);
    A.setValue(1,1,1.);
	CPPUNIT_ASSERT_EQUAL( true, A.isSymmetric(1.e-10) );

	/* Eigenvalues and eigenvectors */
	std::vector< double > vp = A.getEigenvalues(2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, vp[0],1.e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 2, vp[1],1.e-5);

	std::vector< Vector > Vp = A.getEigenvectors(2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(Vp[0][0]) - abs(Vp[0][1]),1.e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(Vp[1][0]) - abs(Vp[1][1]),1.e-5);

	MEDCoupling::DataArrayDouble * VpArrayDouble = A.getEigenvectorsDataArrayDouble(2);
	const double *values=VpArrayDouble->getConstPointer();
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(values[0]) - abs(values[1]),1.e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(values[2]) - abs(values[3]),1.e-5);

    A.setValue(0,0,-1.);
    A.setValue(0,1, 1.);
    A.setValue(1,0, 1.);
    A.setValue(1,1,-1.);

	std::vector< double > sigma = A.getSingularValues(2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, sigma[0],1.e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 2, sigma[1],1.e-5);

	std::vector< Vector > Vs = A.getSingularVectors(2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(Vs[0][0]) - abs(Vs[0][1]),1.e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, abs(Vs[1][0]) - abs(Vs[1][1]),1.e-5);

	/* Condition number of a symmetric non singular matrix */
    A.setValue(0,0, 0.);
    A.setValue(0,1, 1.);
    A.setValue(1,0, 1.);
    A.setValue(1,1, 0.);

	CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, A.getConditionNumber(1.e-10),1.e-5);

	/* Condition number of a symmetric singular matrix */
	SparseMatrixPetsc A33(3,3,2);
    A33.setValue(1,2, 1.);
    A33.setValue(2,1, 1.);

	CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, A33.getConditionNumber(true,1.e-10),1.e-5);

	/* Condition number of a non symmetric singular matrix */
    A33.setValue(1,2, 1.);
    A33.setValue(2,1, 4.);

	CPPUNIT_ASSERT_DOUBLES_EQUAL( 4, A33.getConditionNumber(true,1.e-10),1.e-5);

	SparseMatrixPetsc A4(4,4);
    A4.setValue(0,0,1.);
    A4.setValue(0,1,2.);
    A4.setValue(0,2,3.);
    A4.setValue(0,3,4.);
    A4.setValue(1,0,5.);
    A4.setValue(1,1,6.);
    A4.setValue(1,2,7.);
    A4.setValue(1,3,8.);
    A4.setValue(2,0,9.);
    A4.setValue(2,1,10.);
    A4.setValue(2,2,11.);
    A4.setValue(2,3,12.);
    A4.setValue(3,0,13.);
    A4.setValue(3,1,14.);
    A4.setValue(3,2,15.);
    A4.setValue(3,3,16.);

	A4.viewMatrix();//Display matrix coefficients on the screen
	A4.viewMatrix(true,0.05, "A4");//Open an x windows displaying the matrix nonzero strcture
    //The following line would pause the x window until the user presses right mouse : left mouse->zoom in, middle mouse->zoom out, right mouse->continue with the simulation
	//A4.viewMatrix(true,-1);//This pauses the x window until the user presses right mouse
	A4.getEigenvalues(    4, EPS_SMALLEST_MAGNITUDE, 1e-6, EPSKRYLOVSCHUR, true, 0.05, "A4");//Plot eigenvalues in a X-Windows and write the image in a file
	A4.getSingularValues( 4, SVD_SMALLEST          , 1e-6, SVDCYCLIC     , true, 0.05, "A4");//Plot eigenvalues in a X-Windows and write the image in a file
	
    SparseMatrixPetsc A5(A4.transpose());
	CPPUNIT_ASSERT_EQUAL( 1., A5(0,0) );
	CPPUNIT_ASSERT_EQUAL( 5., A5(0,1) );
	CPPUNIT_ASSERT_EQUAL( 9., A5(0,2) );
	CPPUNIT_ASSERT_EQUAL( 13., A5(0,3) );
	CPPUNIT_ASSERT_EQUAL( 2., A5(1,0) );
	CPPUNIT_ASSERT_EQUAL( 6., A5(1,1) );
	CPPUNIT_ASSERT_EQUAL( 10., A5(1,2) );
	CPPUNIT_ASSERT_EQUAL( 14., A5(1,3) );
	CPPUNIT_ASSERT_EQUAL( 3., A5(2,0) );
	CPPUNIT_ASSERT_EQUAL( 7., A5(2,1) );
	CPPUNIT_ASSERT_EQUAL( 11., A5(2,2) );
	CPPUNIT_ASSERT_EQUAL( 15., A5(2,3) );
	CPPUNIT_ASSERT_EQUAL( 4., A5(3,0) );
	CPPUNIT_ASSERT_EQUAL( 8., A5(3,1) );
	CPPUNIT_ASSERT_EQUAL( 12., A5(3,2) );
	CPPUNIT_ASSERT_EQUAL( 16., A5(3,3) );

	SparseMatrixPetsc AA(4,4);

	AA.setValue(1,1,1.);
	AA.setValue(3,3,3.);
	AA.setValue(2,2,2.);

	CPPUNIT_ASSERT_EQUAL( 1., AA(1,1) );
	CPPUNIT_ASSERT_EQUAL( 2., AA(2,2) );
	CPPUNIT_ASSERT_EQUAL( 3., AA(3,3) );

	SparseMatrixPetsc BB(4,4,3);

	BB.setValue(1,1,1.);
	BB.setValue(3,3,3.);
	BB.setValue(2,2,2.);

	CPPUNIT_ASSERT_EQUAL( 1., BB(1,1) );
	CPPUNIT_ASSERT_EQUAL( 2., BB(2,2) );
	CPPUNIT_ASSERT_EQUAL( 3., BB(3,3) );

    Matrix A6(2,2);
    A6(0,0)=1.;
    A6(0,1)=2.;
    A6(1,0)=3.;
    A6(1,1)=4.;

    SparseMatrixPetsc A7(2,2);
    A7.setValue(0,0,A6);

	CPPUNIT_ASSERT_EQUAL( 1.0, A7(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A7(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A7(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A7(1,1) );

    A7.addValue(0,0,A6);

	CPPUNIT_ASSERT_EQUAL( 2.0, A7(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A7(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A7(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A7(1,1) );

    SparseMatrixPetsc A8(2,2,2,4);
    A8.setValuesBlocked(0,0,A6);

	CPPUNIT_ASSERT_EQUAL( 1.0, A8(0,0) );
	CPPUNIT_ASSERT_EQUAL( 2.0, A8(0,1) );
	CPPUNIT_ASSERT_EQUAL( 3.0, A8(1,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A8(1,1) );

    A8.addValuesBlocked(0,0,A6);

	CPPUNIT_ASSERT_EQUAL( 2.0, A8(0,0) );
	CPPUNIT_ASSERT_EQUAL( 4.0, A8(0,1) );
	CPPUNIT_ASSERT_EQUAL( 6.0, A8(1,0) );
	CPPUNIT_ASSERT_EQUAL( 8.0, A8(1,1) );
	
	A8.viewMatrix();//Display matrix coefficients on the screen
	A8.viewMatrix(true,0.05, "A8");//Open an x windows displaying the matrix nonzero strcture
    //The following line would pause the x window until the user presses right mouse : left mouse->zoom in, middle mouse->zoom out, right mouse->continue with the simulation
	//A8.viewMatrix(true,-1);//This pauses the x window until the user presses right mouse
	A8.getEigenvalues(    4, EPS_SMALLEST_MAGNITUDE, 1e-6, EPSKRYLOVSCHUR, true, 0.05, "A8");//Plot eigenvalues in a X-Windows and write the image in a file
	A8.getSingularValues( 4, SVD_SMALLEST          , 1e-6, SVDCYCLIC     , true, 0.05, "A8");//Plot eigenvalues in a X-Windows and write the image in a file
}
