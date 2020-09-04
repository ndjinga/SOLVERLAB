/*
 * IJKFace.hxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#ifndef IJKFACE_HXX_
#define IJKFACE_HXX_

/**
 * IJKFace class is defined by
 * - indices of nodes surrounding this face
 * - indices of cells surrounding this face
 * - normal vector to this face
 * - measure of this face
 * - barycenter of this face
 */

#include "IntTab.hxx"
#include "Point.hxx"

#include <vector>
#include <string>

class IJKFace
{
public: //----------------------------------------------------------------
    /**
     * default constructor
     */
    IJKFace ( void ) ;

    /**
     * constructor by copy
     * @param face : The Face object to be copied
     */
    IJKFace ( const Face& face ) ;

    /**
     * constructor with data
     * @param numberOfNodes : The number of nodes allocated for this face
     * @param numberOfCells : The number of cells allocated for this face
     * @param measure : The measure of this cell
     * @param p : The barycenter of this cell
     * @param xN : x coordinate N
     * @param yN : y coordinate N
     * @param zN : z coordinate N
     */
    IJKFace( const Point p, std::vector< int > IJKCoords, int faceIndex, const std::vector< int > surroundingNodes, const std::vector< int > surroundingCells, const double measure, Vector normalVector, std::vector<std::string> groupNames ) ;

    /**
     * destructor
     */
    ~IJKFace ( void ) ;

    /**
     * The cells ID that this face belongs to
     * @return _cellsId
     */
    const std::vector< int >  getCellsId ( void ) const ;

    /**
     * The nodes ID that this face belongs to
     * @return _nodesId
     */
    const std::vector< int >  getNodesId ( void ) const ;

    /**
     * return the measure of this face (length in 2D or
     * or surface in 3D
     * @return _measure
     */
    double getMeasure( void ) const ;

    /**
     * return number of cells in this face
     * @return _numberOfCells
     */
    int getNumberOfCells ( void ) const ;

    /**
     * return number of nodes in this face
     * @return _numberOfNodes
     */
    int getNumberOfNodes ( void ) const ;

    /**
     * return barrycenter of this face
     * @return _point
     */
    Point getBarryCenter( void ) const ;

    /**
     * return cordinate x of the barycenter in this face
     */
    double x ( void ) const ;

    /**
     * return cordinate y of the barycenter in this face
     */
    double y ( void ) const ;

    /**
     * return cordinate z of the barycenter in this face
     */
    double z ( void ) const ;

    /**
     * return the list of group names of this face
     */
    std::vector<std::string> getGroupNames(void) const;

    /**
     * return a groupe name of this face
     */
    std::string getGroupName(int igroup=0) const;

    /**
     * return 0 if the face is on the border of domain
     * else -1
     */
    int getRegion(void) const ;

    /**
     * return True if the face is on the border of domain
     * else False
     */
    bool isBorder(void) ;

    /**
     * surcharge operator =
     * @param face : The Face object to be copied
     */
    const Face& operator= ( const Face& face ) ;

	/**
	 * return normal vectors of this face
	 */
     
	Vector getNormalVectors (void) const ;

    double getXN(void) const ;

    double getYN(void) const ;

    double getZN(void) const ;

    int getNodeId(int localId) const ;
    
    int getCellId(int localId) const ;


private: //----------------------------------------------------------------

    double _xN;

    double _yN;

    double _zN;

	/*
	 * The (igrid, i,j,k) coordinate of this face.
	 */
    std::vector< int > _IJKCoords;
    
	/*
	 * The index of this face.
	 */
    int _faceIndex;
    
    /*
     * The number of cells surrounding this face.
     */
    int _numberOfCells ;

    /*
     * The length of this face.
     */
    double _measure ;

    /*
     * The coordinate of barycenter the cell.
     */
    Point _point ;

    /*
     * The region of this face. -1 internal or number of edge that this face belongs to
     */
    int _region ;

    /*
     * The group names of the face.
     */
    std::vector<std::string> _groupNames ;

	/*
	 * The mesh MEDCouplingIMesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingIMesh> _mesh;

    /*
     * The indices of cells surrounding this face.
     */
    const std::vector< int >  _surroundingCells ;

    /*
     * The indices of nodes surrounding this face.
     */
    const std::vector< int >  _surroundingNodes ;
};

#endif /* IJKFACE_HXX_ */
