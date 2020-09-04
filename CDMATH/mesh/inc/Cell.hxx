/*
 * cell.hxx
 *
 *  Created on: 23 January. 2012
 *      Authors: CDMATH
 */

#ifndef CELL_HXX_
#define CELL_HXX_

/**
 * Cell class is defined by
 * - indices of nodes surounding this cell
 * - indices of faces surounding this cell
 * - measure of this cell
 * - barycenter of this cell
 */

#include <vector>
#include "Point.hxx"
#include "Vector.hxx"

class Cell
{
    public: //----------------------------------------------------------------
	/**
	 * default constructor
	 */
	Cell ( void ) ;

	/**
	 * constructor with data
	 * @param numberOfNodes : The number of nodes allocated for this cell
	 * @param numberOfFaces : The number of faces allocated for this cell
	 * @param measure : The measure of this cell
	 * @param p : The barycenter of this cell
	 */
	Cell ( int numberOfNodes, int numberOfFaces, double measure, const Point p ) ;

	/**
	 * constructor by copy
	 * @param cell : The cell object to be copied
	 */
	Cell ( const Cell& cell ) ;

	/**
	 * destructor
	 */
	~Cell ( void ) ;

	/**
	 * return number of nodes in this cell
	 * @return _numberOfNodes
	 */
	int getNumberOfNodes ( void ) const ;

	/**
	 * return number of faces in this cell
	 * @return _countFace
	 */
	int getNumberOfFaces ( void ) const ;

	/**
	 * return nodes ID in this cell
	 * @return _nodes
	 */
	std::vector< int > getNodesId ( void ) const ;

	/**
	 * return faces ID in this cell
	 * @return _faces
	 */
	std::vector< int > getFacesId ( void ) const ;

	/**
	 * return cordinate numComposant of the normal vector in this cell
	 */
	double getNormalVector ( const int numNormalVector, const int numComposant ) const ;

	/**
	 * return normal vectors in this cell
	 */
	Vector getNormalVectors (void) const ;

	/**
	 * return the measure of this cell (length in 1D, surface in 2D or
	 * or volume in 3D
	 * @return _measure
	 */
	double getMeasure ( void ) const ;

	/**
	 * return barrycenter in this cell
	 * @return _point
	 */
	Point getBarryCenter ( void ) const ;

	/**
	 * return cordinate x of the barycenter in this cell
	 */
	double x ( void ) const ;

	/**
	 * return cordinate y of the barycenter in this cell
	 */
	double y ( void ) const ;

	/**
	 * return cordinate z of the barycenter in this cell
	 */
	double z ( void ) const ;

	/**
	 * add a face faceId in this cell
	 */
	void addFaceId ( int numFace, int faceId ) ;

	/**
	 * add a node nodeId in this cell
	 */
	void addNodeId ( int numNode, int nodeId ) ;

	/**
	 * surcharge opertor =
	 * @param cell : The cell object to be copied
	 */
	const Cell& operator= ( const Cell& cell ) ;

	/**
	 * add a normal vector numNormalVector in this cell
	 */
	void addNormalVector ( int numNormalVector, double x, double y, double z ) ;

	int getFaceId(int localId) const ;

	int getNodeId(int localId) const ;

   private: //----------------------------------------------------------------

	/*
	 * The nodes ID in this cell.
	 */
	std::vector< int > _nodesId ;

	/*
	 * The number of nodes in this cell.
	 */
	int _numberOfNodes ;

	/*
	 * The faces ID in this cell.
	 */
	std::vector< int > _facesId ;

	/*
	 * The number of faces in this cell.
	 */
	int _numberOfFaces ;

	/*
	 * The length or surface or volume of the cell.
	 */
	double _measure ;

	/*
	 * The coordinate of barycenter the cell.
	 */
	Point _point ;

	/*
	 * The coordinate of normal vector the cell.
	 */
	Vector _normalVectors ;
};

#endif /* CELL_HXX_ */
