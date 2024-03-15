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
 * - normal vector to the face
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
     * The cells ID that this face belongs to (real geometric neighbours, excluding periodic boundary pairing)
     * @return _cellsId
     */
    std::vector< int > getCellsId ( void ) const ;

    /**
     * The nodes ID that this face belongs to (real geometric neighbours, excluding periodic boundary pairing)
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
     * return number of cells in this face (number of real geometric neighbours, excluding periodic boundary pairing)
     * @return _numberOfCells
     */
    int getNumberOfCells ( void ) const ;

    /**
     * return number of nodes in this face (number of real geometric neighbours, excluding periodic boundary pairing)
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
    void setGroupName(const std::string groupName, bool belongToInnerWall = false);

    /**
     * return true if the face is on an inner wall
     * else false
     */
    bool belongToInnerWall(void) const ;

    /**
     * return true if the face is associated to a periodic boundary condition
     * else false
     */
    bool isPeriodicFace(void) const;

    /**
     * return True if the face is on the geometric border of domain or physically on an inner wall  (faces associated to periodic boundary conditions are not considered proper borders)
     * else False
     */
    bool isBorder(void) ;

    /**
     * @param numCell : local index of cell to add in this face
     * @param cellId : global index of cell to add in this face
     */
    void addCellId (const int numCell, const int cellId ) ;

    /** Add a ghost cell ie a cell that is not geometrically connected to the face but is associated though a periodic boundary condition 
     * @param cellId : global index of the ghost cell to add in this face
     */
    void addPeriodicCellId ( const int cellId ) ;

    /**
     * @param numNode : local index of node to add in this face
     * @param nodeId : global index of node to add in this face
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

    /**
     * @param localId : local index of surrounding node
     * @return global index of neighbouring node with local index localId
     */
    int getNodeId(int localId) const ;
    
    /**
     * @param localId : local index of neighbouring cell
     * @return global index of neighbouring cell with local index localId
     */
    int getCellId(int localId) const ;

    /**
     * @return ghost cell associated to this face (return an exception of this face is not associated to a periodic boundary condition
     */
    int getPeriodicCellId() const ;

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
     * The group names of the face.
     */
    std::vector<std::string> _groupNames ;

    /*
     * This is to manage internal wall boundary conditions where a face is geometrically surrounded by two cells but is physically considered an inner wall.
     */
    bool _belongToInnerWall ;

    /*
     * This is to manage periodic boundary conditions where a face is geometrically surrounded by a single cell but has a ghost neighbour cell because of the periodic condition imposed at the boundary.
     * The ghost neighbour is another boundary cell stored at index numCell (not numCell-1) of the _cellsId vector
     * The vector _cellsId stores one or two cells that are the numCell geometric neighbours and possibly a ghost neighbour if a periodic boundary condition is applied on the face.
     */
    bool _isPeriodicFace ;

};

#endif /* FACE_HXX_ */
