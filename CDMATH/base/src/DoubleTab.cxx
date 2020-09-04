/*
 * DoubleTab.cxx
 *
 *  Created on: 22 févr. 2013
 *      Author: mekkas
 */



#include "DoubleTab.hxx"
#include <string.h>

using namespace std;

DoubleTab::~DoubleTab()
{
	delete [] _values;
}

DoubleTab::DoubleTab()
{
	_numberOfElements=0;
	_values=NULL;
}

DoubleTab::DoubleTab(int size)
{
	_numberOfElements=size;
	_values = new double [size];
}

DoubleTab::DoubleTab(const int size, double initialValue)
{
	_values = new double [size];
	_numberOfElements=size;
	for (int i=0;i<size;i++)
		_values[i] = initialValue ;
}

DoubleTab::DoubleTab(const int size, const double* value)
{
	_values = new double [size];
	_numberOfElements=size;
	memcpy(_values,value,_numberOfElements*sizeof(double)) ;
}

void
DoubleTab::resize(const int size)
{
	double* oldvalues=new double [_numberOfElements] ;
	int oldsize=_numberOfElements;
	copy(_values,_values+oldsize,oldvalues);
	delete [] _values;

	_numberOfElements+=size;
	_values = new double [_numberOfElements] ;
	for (int i=0;i<oldsize;i++)
		_values[i] = oldvalues[i] ;
	for (int i=oldsize;i<_numberOfElements;i++)
		_values[i] = 0. ;
	delete [] oldvalues;
}

DoubleTab::DoubleTab(const DoubleTab& dt)
{
	_numberOfElements=dt.size();
	_values = new double [_numberOfElements];
	memcpy(_values,dt.getValues(),_numberOfElements*sizeof(double)) ;
}

DoubleTab&
DoubleTab::operator=(const DoubleTab & dt)
{
	_numberOfElements=dt.size();
    if (_values)
        delete [] _values ;
	_values = new double [_numberOfElements];
	memcpy(_values,dt.getValues(),_numberOfElements*sizeof(double)) ;
	return *this;
}

DoubleTab&
DoubleTab::operator=(double value)
{
	for (int i=0;i<_numberOfElements;i++)
		_values[i] = value ;
	return *this;
}

double&
DoubleTab::operator[](int i)
{
	return _values[i];
}

const double&
DoubleTab::operator[](int i) const
{
	return _values[i];
}

double&
DoubleTab::operator()(int i)
{
	return _values[i];
}

const double&
DoubleTab::operator()(int i) const
{
	return _values[i];
}

int
DoubleTab::size() const
{
	return _numberOfElements;
}

double*
DoubleTab::getValues(void) const
{
	return _values;
}

double
DoubleTab::max() const
{
	double res=_values[0];
	for (int i=0;i<_numberOfElements;i++)
		if (_values[i]>res)
			res=_values[i];
	return res;
}

double
DoubleTab::min() const
{
	double res=_values[0];
	for (int i=0;i<_numberOfElements;i++)
		if (_values[i]<res)
			res=_values[i];
	return res;
}

double*
DoubleTab::getPointer(void)
{
	return _values;
}

DoubleTab&
DoubleTab::operator+=(const DoubleTab& dt)
{
    if(size() != dt.size())
    {
        cout<<"Warning : adding DoubleTab of different sizes"<<endl;
        cout<<"First DoubleTab size "<<size()<<", second DoubleTab size "<<dt.size()<<endl;
    }
	for (int i=0;i<size();i++)
		_values[i] += dt[i] ;
	return *this;
}

DoubleTab&
DoubleTab::operator+=(double value)
{
	for (int i=0;i<_numberOfElements;i++)
		_values[i] += value ;
	return *this;
}

DoubleTab&
DoubleTab::operator*=(double value)
{
	for (int i=0;i<_numberOfElements;i++)
		_values[i] *= value ;
	return *this;
}

DoubleTab&
DoubleTab::operator/=(double value)
{
	for (int i=0;i<_numberOfElements;i++)
		_values[i] /= value ;
	return *this;
}

DoubleTab&
DoubleTab::operator-=(const DoubleTab& dt)
{
    if(size() != dt.size())
    {
        cout<<"Warning : subtracting DoubleTab of different sizes"<<endl;
        cout<<"First DoubleTab size "<<size()<<", second DoubleTab size "<<dt.size()<<endl;
    }
	for (int i=0;i<size();i++)
		_values[i] -= dt[i] ;
	return *this;
}

DoubleTab&
DoubleTab::operator-=(double value)
{
	for (int i=0;i<_numberOfElements;i++)
		_values[i] -= value ;
	return *this;
}

DoubleTab
operator+(const DoubleTab& U,const DoubleTab& V)
{
	int size = U.size();
	DoubleTab res(size);
	for (int i=0; i<size; i++)
         res[i] = U[i] + V[i];
	return res;
}

DoubleTab
operator-(const DoubleTab& U,const DoubleTab& V)
{
	int size = U.size();
	DoubleTab res(size);
	for (int i=0; i<size; i++)
         res[i] = U[i] - V[i];
	return res;
}

double
operator*(const DoubleTab& U,const DoubleTab& V)
{

	int size = U.size();
	double res = 0.0;
	for(int i=0; i<size; i++)
		res += U[i] * V[i];
   return res;
}

DoubleTab
operator*(double value, const DoubleTab& U)
{

	int size = U.size();
	DoubleTab res(size);
	for(int i=0; i<size; i++)
		res[i] = U[i] * value;
   return res;
}

DoubleTab
operator*(const DoubleTab& U, double value)
{

	int size = U.size();
	DoubleTab res(size);
	for(int i=0; i<size; i++)
		res[i] = U[i] * value;
   return res;
}

DoubleTab
operator/(const DoubleTab& U, double value)
{

	int size = U.size();
	DoubleTab res(size);
	for(int i=0; i<size; i++)
		res[i] = U[i] / value;
   return res;
}

ostream&
operator<<(ostream& out, const DoubleTab& U)
{
	for (int i=0; i<U.size();i++)
	{
		out.width(6);
		out.precision(6);
		out<<U[i];
		out<<endl;
	}
	return out;
}
