/*
 * Matrix.cxx
 *
 *  Created on: 13 April. 2013
 *      Authors: CDMATH
 */

#include <math.h> 

#include "GenericMatrix.hxx"
#include "CdmathException.hxx"

using namespace std;

GenericMatrix::GenericMatrix()
{
	_numberOfRows = 0;
	_numberOfColumns = 0;
	_isSparseMatrix=false;
}

GenericMatrix::~GenericMatrix()
{
}

bool
GenericMatrix::isSparseMatrix( void ) const
{
	return _isSparseMatrix;
}


int
GenericMatrix::getNumberOfRows() const
{
	return _numberOfRows ;
}

int
GenericMatrix::getNumberOfColumns() const
{
	return _numberOfColumns ;
}

const DoubleTab&
GenericMatrix::getValues( void ) const
{
	return _values;
}

//----------------------------------------------------------------------
DoubleTab
GenericMatrix::getValues()
//----------------------------------------------------------------------
{
	return _values;
}

void
GenericMatrix::setValues(const DoubleTab& values)
{
	_values=values;
}

bool
GenericMatrix::isSymmetric(double tol) const
{
	if( ! isSquare() )
		throw "isSymmetric::Matrix is not square!!!";

	bool res = true;

	int dim = _numberOfRows;

	for(int i=0; i<dim-1; i++)
		for(int j=i+1; j<dim; j++)
			if(fabs((*this)(i,j) - (*this)(j,i))> tol )
			{
				res = false;
				break;
			}
	return res;
}

bool
GenericMatrix::isSquare() const
{
	if(_numberOfRows == _numberOfColumns)
		return true;
	return false;
}

int
GenericMatrix::coefficient(int index) const
{
	if(! (index % 2) )
		return (1);
	return (-1);
}


double 
GenericMatrix::max() const
{
    return _values.max();
}


double 
GenericMatrix::min() const
{
    return _values.max();	
}


void
GenericMatrix::view() const
{
	for (int i=0; i<_numberOfRows;i++)
	{
		for (int j=0;j<_numberOfColumns; j++)
		{
			cout.width(6);
			cout.precision(6);
			cout<<(*this)(i,j);
		}
		cout<<endl;
	}
}
