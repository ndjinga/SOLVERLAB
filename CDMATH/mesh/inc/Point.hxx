/*
 * point.hxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#ifndef POINT_HXX_
#define POINT_HXX_

#include <iostream>

/**
 * Point class is defined by
 * - The x-coordinate
 * - The y-coordinate
 * - The z-coordinate
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
       * Create a point at (x, y, z)
       * @param x : The x-coordinate
       * @param y : The y-coordinate
       * @param z : The z-coordinate
       */
      Point( const double x, const double y, const double z ) ;

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
      double& operator[] ( int i ) ;

     /**
      * return x-coordinate
      * @return x-coordinate
      */
      double x () const ;

     /**
      * return y-coordinate
      * @return y-coordinate
      */
      double y () const ;

     /**
      * return z-coordinate
      * @return z-coordinate
      */
      double z () const ;

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

     /**
      * return Compute norm of point representing a vector from the origin
      * @return The norm
      */
      double norm( void ) const ;

     /**
      * return Compute dot product with given vector
      * @param p : Point
      * @return The dot product with given vector
      */
      double dot(const Point& p) const;

    private: //----------------------------------------------------------------
      double _x[3] ;
};

#endif /* POINT_HXX_ */
