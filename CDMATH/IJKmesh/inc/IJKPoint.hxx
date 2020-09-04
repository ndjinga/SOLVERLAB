/*
 * IJKPoint.hxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#ifndef IJKPOINT_HXX_
#define IJKPOINT_HXX_

#include <iostream>

/**
 * Point class is defined by
 * - The x-index
 * - The y-index
 * - The z-index
 */

class Point
{
    public: //----------------------------------------------------------------
      /**
       * default constructor
       */
      Point ( void ) ;

      /**
       * constructor with data
       * Create a point with integer indices at (I, J, K)
       * @param i : The x-index
       * @param j : The y-index
       * @param k : The z-index
       */
      Point( const int I, const int J, const int K ) ;

      /**
       * constructor by copy
       * @param p : The Point object to be copied
      */
      Point ( const Point & p ) ;

      /**
       * destructor
       */
      ~Point ( void ) ;

     /**
      * return address of coordinate in direction i
      * @param i : Direction
      * @return The address of coordinate in the given direction
      */
      int& operator[] ( int i ) ;

     /**
      * return x-index
      * @return x-index
      */
      int I () const ;

     /**
      * return y-index
      * @return y-index
      */
      double J () const ;

     /**
      * return z-index
      * @return z-index
      */
      double K () const ;

     /**
      * return Compute sum of two points
      * @param p : Point
      * @return The Compute sum of two points
      */
      Point operator+ ( const Point& p ) const ;

     /**
      * return Compute difference of two points
      * @param p : Point
      * @return The Compute difference of two points
      */
      Point operator- ( const Point& p ) const ;

     /**
      * return Add given point
      * @param p : Point
      * @return Add given point
      */
      const Point& operator+= ( const Point& p ) ;

     /**
      * return Subtract given point
      * @param p : Point
      * @return Subtract given point
      */
      const Point& operator-= ( const Point& p ) ;

     /**
      * return Multiplication with scalar
      * @param s : Scalar
      * @return Multiplication with scalar
      */
      Point operator* ( double s ) const ;

     /**
      * return Incremental multiplication with scalar
      * @param s : Scalar
      * @return Incremental multiplication with scalar
      */
      const Point& operator*= ( double s ) ;

     /**
      * return Division with scalar
      * @param s : Scalar
      * @return Division with scalar
      */
      Point operator/ ( double s ) const ;

     /**
      * return Incremental Division with scalar
      * @param s : Scalar
      * @return Incremental Division with scalar
      */
      const Point& operator/= ( double s ) ;

     /**
      * return Assignment operator
      * @param p : Point
      * @return Assignment operator
      */
      const Point& operator= ( const Point& p ) ;

     /**
      * return Compute distance to given point
      * @param p : Point
      * @return The distance
      */
      double distance( const Point& p ) const ;


    private: //----------------------------------------------------------------
      int _x[3] ;
};

#endif /* POINT_HXX_ */
