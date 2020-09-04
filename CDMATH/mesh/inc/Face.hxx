/*
 * face.hxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#ifndef FACE_HXX_
#define FACE_HXX_

/**
 * face class is defined by
 * - number of nodes allocated for this face
 * - number of cells allocated for this face
 * - measure of this face
 * - barycenter of this face
 */

#include "Point.hxx"

#include <vector>
#include <string>

class Face
{
public: //----------------------------------------------------------------
    /**
     * default constructor
     */
    Face ( void ) ;

    /**
     * constructor by copy
     * @param face : The Face object to be copied
     */
    Face ( const Face& face ) ;

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
    Face( const int numberOfNodes, const int numberOfCells, const double measure, const Point p, double xN, double yN, double zN ) ;

    /**
     * destructor
     */
    ~Face ( void ) ;

    /**
     * The cells ID that this face belongs to
     * @return _cellsId
     */
    std::vector< int > getCellsId ( void ) const ;

    /**
     * The nodes ID that this face belongs to
     * @return _nodesId
     */
    std::vector< int > getNodesId ( void ) const ;

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
     * @param groupName : set a groupe name for this face
     */
    void setGroupName(const std::string groupName);

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
     * @param numCell : index local of cell to add in this face
     * @param cellId : index global of cell to add in this face
     */
    void addCellId (const int numCell, const int cellId ) ;

    /**
     * @param numNode : index local of node to add in this face
     * @param nodeId : index global of node to add in this face
     */
    void addNodeId (const int numNode, const int nodeId ) ;

    /**
     * surcharge operator =
     * @param face : The Face object to be copied
     */
    const Face& operator= ( const Face& face ) ;

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
     * The cell id that this face belongs to.
     */
    std::vector< int > _cellsId ;

    /*
     * The vertex id that this face belongs to.
     */
    std::vector< int > _nodesId ;

    /*
     * The number of cells allocated for this face.
     */
    int _numberOfCells ;

    /*
     * The number of nodes allocated for this face.
     */
    int _numberOfNodes ;

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
};

#endif /* FACE_HXX_ */
