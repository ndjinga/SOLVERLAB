/*
 * DoubleTab.hxx
 *
 *  Created on: 22 févr. 2013
 *      Author: mekkas
 */

#ifndef DOUBLETAB_HXX_
#define DOUBLETAB_HXX_

#include <iostream>

class DoubleTab
{
	public://----------------------------------------------------------------

	DoubleTab();

	DoubleTab(const int size);

	DoubleTab(const int size, double initialValue);

	DoubleTab(const DoubleTab& dt);

	DoubleTab(const int size, const double* value) ;

	void resize(const int size) ;

	~DoubleTab();

	int size() const;

	double* getPointer(void) ;

	double* getValues(void) const ;

	double max() const ;

	double min() const ;

	DoubleTab & operator=(const DoubleTab & dt);

	DoubleTab & operator=(double value);

	double & operator[](int i);

	const double & operator[](int i) const;

	double & operator()(int i);

	const double & operator()(int i) const;

	DoubleTab& operator+=(const DoubleTab& dt);

	DoubleTab& operator+=(double value);

	DoubleTab& operator-=(const DoubleTab& dt);

	DoubleTab& operator-=(double value);

	DoubleTab& operator*= (double value) ;

	DoubleTab& operator/= (double value) ;

	friend DoubleTab operator+ (const DoubleTab& U, const DoubleTab& V) ;

	friend DoubleTab operator- (const DoubleTab& U, const DoubleTab& V) ;

    friend DoubleTab operator*(double value , const DoubleTab& V) ;

    friend DoubleTab operator*(const DoubleTab& U, double value ) ;

	friend DoubleTab operator/(const DoubleTab& U, double value ) ;

	friend double  operator* (const DoubleTab& U, const DoubleTab& V) ;

	friend std::ostream& operator<<(std::ostream& out, const DoubleTab& U ) ;

	private://----------------------------------------------------------------
	/*
	 * Values.
	 */
	double* _values;

	/*
	 * number of elements.
	 */
	int _numberOfElements;
};


#endif /* DOUBLETAB_HXX_ */
