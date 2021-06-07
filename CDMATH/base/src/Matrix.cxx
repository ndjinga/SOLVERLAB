/*
 * Matrix.cxx
 *
 *  Created on: 13 April. 2013
 *      Authors: CDMATH
 */
#include <iostream>
#include <cstring>

#include "Matrix.hxx"
#include "Vector.hxx"

#include "CdmathException.hxx"

using namespace std;

Matrix::Matrix()
{
	_isSparseMatrix=false;
}

Matrix::Matrix(int dim)
{
	_numberOfRows = dim;
	_numberOfColumns = dim;
	_isSparseMatrix=false;
	_values=DoubleTab(_numberOfRows*_numberOfColumns,0.);
}

Matrix::Matrix(int numberOfRows, int numberOfColumns)
{
	_numberOfRows = numberOfRows;
	_numberOfColumns = numberOfColumns;
	_isSparseMatrix=false;
	_values=DoubleTab(_numberOfRows*_numberOfColumns,0.);
}

Matrix::~Matrix()
{
}

Matrix::Matrix(const Matrix& matrix)
{
	_numberOfRows = matrix.getNumberOfRows();
	_numberOfColumns = matrix.getNumberOfColumns();
	_isSparseMatrix=matrix.isSparseMatrix();
	_values=DoubleTab (_numberOfRows*_numberOfColumns,matrix.getValues().getValues());
}

const DoubleTab&
Matrix::getValues( void ) const
{
	return _values;
}

//----------------------------------------------------------------------
DoubleTab
Matrix::getValues()
//----------------------------------------------------------------------
{
	return _values;
}

void
Matrix::setValues(const DoubleTab& values)
{
	_values=values;
}

double
Matrix::max() const
{
    return _values.max();
}

double 
Matrix::min() const
{
    return _values.max();	
}

std::vector< double > 
Matrix::getArray() 
{
	int numberOfRows  =getNumberOfRows();
	int numberOfColums=getNumberOfColumns();
	int size=_numberOfRows*numberOfColums;
	
	vector< double >  result(size);	
	double* values = result.data();
	
    memcpy(values,_values.getPointer(),size*sizeof(double)) ;

	return result;
}

bool
Matrix::isSparseMatrix( void ) const
{
	return _isSparseMatrix;
}

double&
Matrix::operator()(int i, int j)
{
	if(i>_numberOfRows || j>_numberOfColumns || i<0 || j<0)
		throw CdmathException("double& Matrix::operator()(int i, int j) : i>number of rows or j>number of columns !");

	return _values.getPointer()[i+_numberOfRows*j];
}

double
Matrix::operator()(int i, int j) const
{
	if(i>_numberOfRows || j>_numberOfColumns || i<0 || j<0)
		throw CdmathException("double& Matrix::operator()(int i, int j) : i>number of rows or j>number of columns !");
	return _values[i+_numberOfRows*j];
}

Vector
Matrix::operator* (const Vector& vector) const
{
	Vector res(_numberOfRows);
	for(int i=0; i<_numberOfRows; i++)
	{
		double sum=0.;
		for(int j=0; j<_numberOfColumns; j++)
			sum=sum+(*this)(i,j)*vector(j);
		res(i) = sum;
	}
	return res;
}

Matrix&
Matrix::operator*= (const Matrix& matrix)
{
	int numberOfRows2 = matrix.getNumberOfRows();
	int numberOfColumns2 = matrix.getNumberOfColumns();

	if(_numberOfColumns!=numberOfRows2)
	{
		string msg="Matrix Matrix::operator()*(const Matrix& matrix1, const Matrix& matrix2) : dimensions of the matrices is incompatible!";
		throw CdmathException(msg);
	}
	Matrix res(_numberOfRows, numberOfColumns2);
	for(int i=0;i<_numberOfRows;i++)
	{
		for(int j=0;j<numberOfColumns2;j++)
		{
			double som=0.;
			for(int k=0;k<_numberOfColumns;k++)
				som+=(*this)(i,k)*matrix(k,j);
			res(i,j)=som;
		}
	}
	(*this)=res;
	return (*this);
}

Matrix
Matrix::transpose() const
{
	Matrix res(_numberOfColumns, _numberOfRows);
	for(int i=0; i<_numberOfRows; i++)
		for(int j=0; j<_numberOfColumns; j++)
			res(i,j) = (*this)(j,i);
	return res;
}

Matrix
Matrix::partMatrix(int row, int column) const
{
	int r = 0;
	int c = 0;
	Matrix res(_numberOfRows-1, _numberOfColumns-1);

	for (int i=0; i<_numberOfRows; i++)
	{
		c = 0;
		if(i != row)
		{
			for(int j=0; j<_numberOfColumns; j++)
				if(j != column)
					res(r,c++) = (*this)(i,j);
			r++;
		}
	}
	return res;
}

int
Matrix::coefficient(int index) const
{
	if(! (index % 2) )
		return (1);
	return (-1);
}

double
Matrix::determinant() const
{
	if( ! isSquare() )
		throw "isSymmetric::Matrix is not square!!!";
	else
	{
		double res = 0.0;
		int dim = _numberOfRows;
		if(dim==1)
			return (*this)(0,0);

		for(int i=0; i<dim; i++)
		{
			Matrix matrix = this->partMatrix(i,0);
			res += ( coefficient(i)*(*this)(i,0)*(matrix.determinant() ) );
		}
		return res;
	}
}

Matrix
operator+ (const Matrix& matrix1, const Matrix& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="Matrix Matrix::operator()+(const Matrix& matrix1, const Matrix& matrix2) : number of rows or columns of the matrices is diffrerent!";
		throw CdmathException(msg);
	}
	Matrix res(numberOfRows, numberOfColumns);
	for(int i=0;i<numberOfRows;i++)
		for(int j=0;j<numberOfColumns;j++)
			res(i,j)=matrix1(i,j)+matrix2(i,j);
	return res;
}

Matrix
operator- (const Matrix& matrix1, const Matrix& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="Matrix Matrix::operator()+(const Matrix& matrix1, const Matrix& matrix2) : number of rows or columns of the matrices is diffrerent!";
		throw CdmathException(msg);
	}
	Matrix res(numberOfRows, numberOfColumns);
	for(int i=0;i<numberOfRows;i++)
		for(int j=0;j<numberOfColumns;j++)
			res(i,j)=matrix1(i,j)-matrix2(i,j);
	return res;
}

Matrix
operator*(const Matrix& matrix1, const Matrix& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfColumns!=numberOfRows2)
	{
		string msg="Matrix Matrix::operator()*(const Matrix& matrix1, const Matrix& matrix2) : dimensions of the matrices is incompatible!";
		throw CdmathException(msg);
	}
	Matrix res(numberOfRows, numberOfColumns2);
	for(int i=0;i<numberOfRows;i++)
	{
		for(int j=0;j<numberOfColumns2;j++)
		{
			double som=0.;
			for(int k=0;k<numberOfColumns;k++)
				som+=matrix1(i,k)*matrix2(k,j);
			res(i,j)=som;
		}
	}
	return res;
}

Matrix
operator* (double value , const Matrix& matrix )
{
	Matrix res(matrix);
	DoubleTab t1=res.getValues();
	t1*=value;
	res.setValues(t1);
	return res;
}

Matrix
operator* (const Matrix& matrix, double value )
{
	Matrix res(matrix);
	DoubleTab t1=res.getValues();
	t1*=value;
	res.setValues(t1);
	return res;
}

Matrix
operator/ (const Matrix& matrix, double value)
{
	Matrix res(matrix);
	DoubleTab t1=res.getValues();
	t1/=value;
	res.setValues(t1);
	return res;
}

Matrix&
Matrix::operator+= (const Matrix& matrix)
{
	_values+=matrix.getValues();
	return (*this);
}

Matrix&
Matrix::operator-= (const Matrix& matrix)
{
	_values-=matrix.getValues();
	return (*this);
}

Matrix&
Matrix::operator*= (double value)
{
	_values*=value;
	return (*this);
}

Matrix&
Matrix::operator/= (double value)
{
	_values/=value;
	return (*this);
}

//----------------------------------------------------------------------
const Matrix&
Matrix::operator= ( const Matrix& matrix )
//----------------------------------------------------------------------
{
	_numberOfRows=matrix.getNumberOfRows();
	_numberOfColumns=matrix.getNumberOfColumns();
	_isSparseMatrix=matrix.isSparseMatrix();
	_values=DoubleTab (_numberOfRows*_numberOfColumns,matrix.getValues().getValues());
	
    return *this;
}

Matrix Matrix::deepCopy(  ) const
{
    Matrix A(getNumberOfRows(), getNumberOfColumns()) ;
    DoubleTab values(_numberOfRows*_numberOfColumns,getValues().getValues());
    A.setValues(values);
    
    return A;
}

ostream&
operator<<(ostream& out, const Matrix& matrix)
{
	for (int i=0; i<matrix.getNumberOfRows();i++)
	{
		for (int j=0;j<matrix.getNumberOfColumns(); j++)
		{
			out.width(6);
			out.precision(6);
			out<<"\t"<<matrix(i,j);
		}
		out<<endl;
	}
	return out;
}
