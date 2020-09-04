/*
 * IJKCell.hxx
 *
 *  Created on: 24 March. 2019
 *      Authors: CDMATH
 */

#ifndef IJKCELL_HXX_
#define IJKCELL_HXX_

/**
 * IJKCell class is defined by
 * - number of nodes surrounding every cell of the mesh : known constant
 * - number of faces surrounding every cell of the mesh : known constant
 * - number of cells surrounding this cell : two cells are neighbours if they have a common face (neighbour_policy 0), a comon edge (neighbour_policy 1) or a common node (neighbour_policy 2)
 * - measure of this cell : known constant
 * - barycenter of this cell
 * 
 * - methods to determine neighbouring faces and nodes indices from mesh structure
 */

#include "Point.hxx"
#include "IntTab.hxx"
#include "Vector.hxx"

class IJKCell
{
    public: //----------------------------------------------------------------
	/**
	 * default constructor
	 */
	IJKCell ( void ) ;

	/**
	 * constructor with data
	 * @param numberOfNodes : The number of nodes allocated for this cell
	 * @param numberOfFaces : The number of faces allocated for this cell
	 * @param measure : The measure of this cell
	 * @param p : The barycenter of this cell
	 */
	IJKCell ( const Point p, std::vector< int > IJKCoords, int cellIndex, std::vector< int > surroundingNodes, std::vector< int > surroundingFaces, double measure, std::vector< Vector > normalVectorsAroundCells ) ;

	/**
	 * constructor by copy
	 * @param cell : The cell object to be copied
	 */
	IJKCell ( const Cell& cell ) ;

	/**
	 * destructor
	 */
	~IJKCell ( void ) ;

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
	 * return nodes ID surrounding this cell
	 * @return _nodes
	 */
	std::vector< int >  getNodesId ( void ) const ;

	/**
	 * return faces ID surrounding this cell
	 * @return _faces
	 */
	std::vector< int >  getFacesId ( void ) const ;

	/**
	 * return cells ID around this cell depending on neighbouring policy
	 * @param neighbouring_policy : two cells are neighbours if they have a common face (value 0), a comon edge (value 1) or a common node (value 2)
	 */
	std::vector< int >  getSurroundingCellsId ( int neighbouring_policy=0 ) const ;

	/**
	 * return cordinate numComposant of the normal vector in this cell
	 */
	double getNormalVector ( const int numNormalVector, const int numComposant ) const ;

	/**
	 * return normal vectors surrounding this cell
	 */
	std::vector< Vector > getNormalVectors (void) const ;

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
	 * surcharge opertor =
	 * @param cell : The cell object to be copied
	 */
	const IJKCell& operator= ( const IJKCell& cell ) ;

   private: //----------------------------------------------------------------

	/*
	 * The (i,j,k) coordinate of this cell.
	 */
    std::vector< int > _IJKCoords;
	/*
	 * The index of this cell.
	 */
    int _cellIndex;// necessaire ???
    
	/*
	 * The coordinate of barycenter the cell.
	 */
	Point _point ;// necessaire ??? se déduit des cordonnées IJK

    /*
     * The number of cells surrounding this cell.
     */
    int _numberOfCells ;

    /*
     * The indices of cells surrounding this Node.
     */
    const std::vector< int >  _surroundingCells ;

	/*
	 * The normal vectors surrounding each cell of the mesh.
	 */
	std::vector< Vector > _normalVectorsAroundCells ;

	/*
	 * The mesh MEDCouplingIMesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingIMesh> _mesh;

};

#endif /* IJKCELL_HXX_ */
