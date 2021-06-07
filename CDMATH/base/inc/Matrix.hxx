/*
 * Matrix.hxx
 *
 *  Created on: 13 April. 2013
 *      Authors: CDMAT
 */

#ifndef MATRIX_HXX_
#define MATRIX_HXX_


/**
 * Matrix class is defined by
 * - number of rows
 * - number of columns
 * - values array
 */

#include <iostream>

#include "GenericMatrix.hxx"

class Vector ;

class Matrix: public GenericMatrix
{
    public: //----------------------------------------------------------------
    /**
     * default constructor
     */
    Matrix ( void ) ;

    /**
     * constructor with data
     * @param dim : The number of rows and columns
     */
    Matrix ( int dim ) ;

    /**
     * constructor with data
     * @param numberOfRows : The number of rows
     * @param numberOfColumns : The number of columns
     */
    Matrix ( int numberOfRows, int numberOfColumns ) ;

    /**
     * constructor by copy
     * @param matrix : The Matrix object to be copied
     */
    Matrix ( const Matrix& matrix ) ;

    /**
     * deep copy of a matrix (values are copied)
     * @param matrix : The Matrix object to be copied
     */
    Matrix deepCopy(  ) const;

    /**
     * destructor
     */
    virtual ~Matrix ( void ) ;

    const DoubleTab& getValues( void ) const ;

	DoubleTab getValues( void ) ;

	void setValues(const DoubleTab& values) ;

	//returns the array of matrix coefficients
	std::vector< double > getArray();

	//returns the maximum coefficient
    double max() const;

	//returns the minimum coefficient
    double min() const;

    bool isSparseMatrix( void ) const ;

    double& operator () ( int i, int j ) ;

    double operator () ( int i, int j ) const ;

    Matrix& operator+= (const Matrix& matrix) ;

    Matrix& operator-= (const Matrix& matrix) ;

    Matrix& operator*= (double value) ;

    Matrix& operator*= (const Matrix& matrix) ;

    Matrix& operator/= (double value) ;

    Vector operator* (const Vector& vector) const ;

    Matrix transpose() const ;

    Matrix partMatrix(int row, int column) const ;

    double determinant() const ;

    const Matrix& operator= ( const Matrix& matrix ) ;

    friend Matrix operator+ (const Matrix& matrix1, const Matrix& matrix2);

    friend Matrix operator- (const Matrix& matrix1, const Matrix& matrix2);

    friend Matrix operator* (double value , const Matrix& matrix ) ;

    friend Matrix operator* (const Matrix& matrix, double value ) ;

    friend Matrix operator/ (const Matrix& matrix, double value) ;

    friend Matrix operator*(const Matrix& M, const Matrix& N) ;

    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix ) ;


    protected: //----------------------------------------------------------------

    DoubleTab _values ;

	//This function is used in the computation of the determinant
    int coefficient(int index) const ;
};

#endif /* MATRIX_HXX_ */
