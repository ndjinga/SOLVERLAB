/*
 * GenericMatrix.hxx
 *
 *  Created on: 13 April. 2013
 *      Authors: CDMATH
 */

#ifndef GENERICMATRIX_HXX_
#define GENERICMATRIX_HXX_


/**
 * GenericMatrix class is defined by
 * - number of rows
 * - number of columns
 * - values array
 */

#include <iostream>
#include <vector>
#include "DoubleTab.hxx"

class GenericMatrix
{
    public: //----------------------------------------------------------------
    /**
     * default constructor
     */
	GenericMatrix ( void ) ;

    /**
     * destructor
     */
    virtual ~GenericMatrix ( void ) ;

    /**
     * return number of rows in this matrix
     * @return _numberOfRows
     */
    int getNumberOfRows ( void ) const ;

    /**
     * return number of columns in this matrix
     * @return _numberOfColumns
     */
    int getNumberOfColumns ( void ) const ;

    virtual double operator ()( int i, int j ) const = 0;

    virtual bool isSymmetric(double tol=1e-6) const ;

    bool isSquare() const ;

    bool isSparseMatrix( void ) const ;

    virtual  bool containsPetscMatrix() const { return false; };

    void view() const ;

    protected: //----------------------------------------------------------------

    /*
     * The number of rows.
     */
    int _numberOfRows ;

    /*
     * The number of columns.
     */
    int _numberOfColumns ;

    bool _isSparseMatrix ;

};

#endif /* GENERICMATRIX_HXX_ */
