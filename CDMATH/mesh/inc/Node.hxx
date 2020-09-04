/*
 * node.hxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#ifndef NODE_HXX_
#define NODE_HXX_

/**
 * Node class is defined by
 * - The coordinate of the node
 * - number of faces allocated for this node
 * - number of cells allocated for this node
 */

#include "Point.hxx"

#include <vector>
#include <string>

class Node
{
    public: //----------------------------------------------------------------
    /**
    * default constructor
    */
    Node ( void ) ;

    /**
     * constructor with data
     * @param numberOfCells : The number of cells allocated for this node
     * @param numberOfFaces : The number of faces allocated for this node
     * @param p : The barycenter of this node
     */
    Node ( const int numberOfCells, const int numberOfFaces, const int numberOfEdges, const Point p) ;

    /**
     * constructor by copy
     * @param node : The object to be copied
     */
    Node ( const Node & node ) ;

    /**
     * destructor
     */
    ~Node ( void ) ;

    /**
     * return number of cells in this node
     * @return _numberOfCells
     */
    int getNumberOfCells ( void ) const ;

    /**
     * return number of faces in this node
     * @return _numberOfFaces
     */
    int getNumberOfFaces ( void ) const ;

    /**
     * return number of edges in this node
     * @return _numberOfEdges
     */
    int getNumberOfEdges ( void ) const ;

    /**
     * return The coordinate of the Node
     * @return _point
     */
    Point getPoint ( void ) const ;

    /**
     * return cordinate x of this node
     */
    double x ( void ) const ;

    /**
     * return cordinate x of this node
     */
    double y ( void ) const ;

    /**
     * return cordinate x of this node
     */
    double z ( void ) const ;

    /**
     * return the list of group names of this node
     */
    std::vector<std::string> getGroupNames(void) const;

    /**
     * return a specific group name of this node
     */
    std::string getGroupName(int igroup=0) const;

    /**
     *  set a groupe name for this node
     */
    void setGroupName(const std::string groupName);


    /**
     * return 0 if the node is on the border of domain
     * else -1
     */
    int getRegion(void) const ;

    /**
     * return True if the node is on the border of domain
     * else False
     */
    bool isBorder(void) const ;

    /**
     * @param numFace : index local of face to add in this node
     * @param faceId : index global of face to add in this node
     */
    void addFaceId (const int numFace, const int faceId, const bool isBorder=false  ) ;

    /**
     * @param numCell : local index of cell to add in this node
     * @param cellId : global index of cell to add in this node
     */
    void addCellId (const int numCell, const int cellId ) ;

    /**
     * @param numNode : local index of node to add in this node neighbours
     * @param nodeId : global index of node to add in this node neighbours
     */
    void addNeighbourNodeId (const int numNode, const int nodeId ) ;

    /**
     * surcharge opertor =
     * @param node : The node object to be copied
     */
    const Node& operator= ( const Node& node ) ;

    /**
     * return Compute distance to given node
     * @param n : Node
     * @return The distance
     */
    double distance( const Node& n ) const ;

    /**
     * The cells ID that this Node belongs to
     * @return _cellsId
     */
    std::vector< int > getCellsId ( void ) const ;

    /**
     * The faces ID that this Node belongs to
     * @return _facesId
     */
    std::vector< int > getFacesId ( void ) const ;

    int getFaceId(int localId) const ;

    /**
     * The neighbour nodes ID : nodes connected to this node by an edge
     * @return _neighbourNodeId
     */
    std::vector< int > getNeighbourNodesId ( void ) const ;

    int getNeighbourNodeId(int localId) const ;

    private: //----------------------------------------------------------------

    /*
     * The coordinate of the Node.
     */
    Point _point ;

    /*
     * The cells ID that this Node belongs to.
     */
    std::vector< int > _cellsId ;

    /*
     * The faces ID that this Node belongs to.
     */
    std::vector< int > _facesId ;

    /*
     * The neighbour nodes ID that this Node belongs to.
     */
    std::vector< int > _neighbourNodesId ;

    /*
     * The number of cells allocated for this Node.
     */
    int _numberOfCells ;

    /*
     * The number of faces allocated for this Node.
     */
    int _numberOfFaces ;

    /*
     * The number of edges allocated for this Node.
     */
    int _numberOfEdges ;

    /*
     * The region of this face. -1 internal or number of edge that this face belongs to
     */
    int _region ;

    /*
     * The group names of the Node.
     */
    std::vector<std::string> _groupNames ;
    
    /* Is the node a border node ? */
    bool _isBorder;
};

#endif /* NODE_HXX_ */
