/*
 * SparseMatrixPetsc.hxx
 *
 *  Created on: 03/11/2017. 2017
 *      Author: mndjinga
 */

#ifndef SOURCE_DIRECTORY__BASE_INC_SparseMatrixPetsc_HXX_
#define SOURCE_DIRECTORY__BASE_INC_SparseMatrixPetsc_HXX_

#include <iostream>
#include <vector>
#include <cstring>

#include "MEDCouplingUMesh.hxx"
#include "Vector.hxx"
#include <petsc.h>

#include <slepceps.h>
#include <slepcsvd.h>

class SparseMatrixPetsc: public GenericMatrix {

public:
	SparseMatrixPetsc();
	virtual ~SparseMatrixPetsc();

	/**
	 * constructor for a block sparse matrix
	 * @param numberOfRows : The number of rows
	 * @param numberOfColumns : The number of columns
	 */
	SparseMatrixPetsc( int numberOfRows, int numberOfColumns ) ;

	/**
	 * constructor for a block sparse matrix with number of non zero coefficients given
	 * @param numberOfRows : The number of rows
	 * @param numberOfColumns : The number of columns
	 * @param nnz : The maximum number of nonzeros coefficients per line (or an upper bound)
	 */
	SparseMatrixPetsc( int numberOfRows, int numberOfColumns, int nnz ) ;

	/**
	 * constructor for a sparse matrix with block structure
	 * @param blockSize : The block size
	 * @param numberOfRows : The number of rows
	 * @param numberOfColumns : The number of columns
	 * @param nnz : The maximum number of nonzeros coefficients per line (or an upper bound)
     * @comment blockSize should always divide numberOfRows and numberOfColumns
	 */
    SparseMatrixPetsc( int blockSize, int numberOfRows, int numberOfColumns, int nnz );

	/**
	 * constructor by copy
	 * @param SparseMatrixPetsc : The SparseMatrixPetsc object to be copied
	 */
	SparseMatrixPetsc ( const SparseMatrixPetsc& sparseMatrixPetsc ) ;

	/**
	 * constructor with data
	 * @param Petsc matris : mat
	 */
	SparseMatrixPetsc(Mat mat);

	SparseMatrixPetsc transpose() const ;

	double operator()( int i, int j ) const ;

	void setValue( int i, int j, double value ) ;
	void addValue( int i, int j, double value ) ;

	void setValue( int i, int j, Matrix M ) ;
	void addValue( int i, int j, Matrix M ) ;

	void setValuesBlocked( int i, int j, Matrix M ) ;
	void addValuesBlocked( int i, int j, Matrix M ) ;

	SparseMatrixPetsc& operator+= (const SparseMatrixPetsc& SparseMatrixPetsc) ;

	SparseMatrixPetsc& operator-= (const SparseMatrixPetsc& SparseMatrixPetsc) ;

	SparseMatrixPetsc& operator*= (double value) ;

	SparseMatrixPetsc& operator/= (double value) ;

	SparseMatrixPetsc& operator*= (const SparseMatrixPetsc& matrix) ;

	Vector operator* (const Vector& vector) const ;

	const SparseMatrixPetsc& operator= ( const SparseMatrixPetsc& SparseMatrixPetsc ) ;

	friend SparseMatrixPetsc operator+ (const SparseMatrixPetsc& SparseMatrixPetsc1, const SparseMatrixPetsc& SparseMatrixPetsc2);

	friend SparseMatrixPetsc operator- (const SparseMatrixPetsc& SparseMatrixPetsc1, const SparseMatrixPetsc& SparseMatrixPetsc2);

	friend SparseMatrixPetsc operator* (double value , const SparseMatrixPetsc& SparseMatrixPetsc ) ;

	friend SparseMatrixPetsc operator* (const SparseMatrixPetsc& SparseMatrixPetsc, double value ) ;

	friend SparseMatrixPetsc operator/ (const SparseMatrixPetsc& SparseMatrixPetsc, double value) ;

	friend SparseMatrixPetsc operator*(const SparseMatrixPetsc& M, const SparseMatrixPetsc& N) ;

	void viewMatrix(bool useXWindow=false, double pause_lenght=0, std::string matrixName="") const ;
	//Save matrix coefficients into a file in ascii or binary mode
	void saveMatrix(std::string filename, bool binaryMode=false) const ;
	double getMatrixCoeff(int i, int j) const;    

	bool containsPetscMatrix() const;
	Mat getPetscMatrix() const;
	//returns the array of matrix coefficients
	std::vector< double > getArray();
    
    void diagonalShift(double lambda);
    void zeroEntries();//sets the matrix coefficients to zero
    
    std::vector< double > getEigenvalues( int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6, EPSType type = EPSKRYLOVSCHUR, bool viewEigenvaluesInXWindows=false, double pause_lenght=0, std::string matrixName="") const;
    std::vector< Vector > getEigenvectors(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6, EPSType type = EPSKRYLOVSCHUR) const;
    void plotEigenvalues(std::string matrixName="", int nev=-1, double pause_lenght=0, double tol=1e-6, EPSWhich which=EPS_SMALLEST_MAGNITUDE, EPSType type = EPSKRYLOVSCHUR) const;
    MEDCoupling::DataArrayDouble * getEigenvectorsDataArrayDouble(int nev, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6, EPSType type = EPSKRYLOVSCHUR) const;
    std::vector< double > getSingularValues( int nsv, SVDWhich which=SVD_SMALLEST, double tol=1e-6, SVDType type = SVDCYCLIC, bool viewSingularValuesInXWindows=false, double pause_lenght=0, std::string matrixName="") const;
    std::vector< Vector > getSingularVectors(int nsv, SVDWhich which=SVD_SMALLEST, double tol=1e-6, SVDType type = SVDCYCLIC) const;
    double getConditionNumber(bool isSingular=false, double tol=1e-6) const;
        
    bool isSymmetric(double tol=1.e-6) const ;
    
    void  leftDiagonalScale(Vector v);
    void rightDiagonalScale(Vector v);
    
private:
	Mat _mat;

	int _numberOfNonZeros ;//The maximum number of nonzeros coefficients per line (or an upper bound)
	
	int computeSpectrum(int nev, double ** valP, double ***vecP, EPSWhich which=EPS_SMALLEST_MAGNITUDE, double tol=1e-6, EPSType type = EPSKRYLOVSCHUR, bool viewEigenvaluesInXWindows=false, double pause_lenght=0, std::string matrixName="") const;
	int computeSVD     (int nsv, double ** valS, double ***vecS, SVDWhich which=SVD_SMALLEST          , double tol=1e-6, SVDType type = SVDCYCLIC, bool viewSingularValuesInXWindows=false, double pause_lenght=0, std::string matrixName="") const;

	Vector vecToVector(const Vec& vec) const ;
	Vec vectorToVec( const Vector& myVector ) const ;
};


#endif /* SOURCE_DIRECTORY__BASE_INC_SparseMatrixPetsc_HXX_ */
