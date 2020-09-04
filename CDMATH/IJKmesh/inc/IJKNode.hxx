/*
 * IJKNode.hxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#ifndef IJKNODE_HXX_
#define IJKNODE_HXX_

/**
 * IJKNode class is defined by
 * - The coordinates of the node
 * - number of faces surrounding this node
 * - number of cells surrounding this node
 * - number of edges (or nodes) surrounding this node
 */

#include "Point.hxx"
#include "IntTab.hxx"

#include <vector>
#include <string>

class IJKNode
{
    public: //----------------------------------------------------------------
    /**
    * default constructor
    */
    IJKNode ( void ) ;

    /**
     * constructor with data
     * @param numberOfCells : The number of cells allocated for this node
     * @param numberOfFaces : The number of faces allocated for this node
     * @param p : The barycenter of this node
     */
    IJKNode ( const Point p, std::vector< int > IJKCoords, int nodeIndex, const std::vector< int > surroundingCells, const std::vector< int > surroundingFaces, const std::vector< int > surroundingNodes, std::vector<std::string> groupNames) ;

    /**
     * constructor by copy
     * @param node : The object to be copied
     */
    IJKNode ( const Node & node ) ;

    /**
     * destructor
     */
    ~IJKNode ( void ) ;

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
     * surcharge opertor =
     * @param node : The node object to be copied
     */
    const IJKNode& operator= ( const Node& node ) ;

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
    std::vector< int >  getCellsId ( void ) const ;

    /**
     * The faces ID that this Node belongs to
     * @return _facesId
     */
    std::vector< int >  getFacesId ( void ) const ;

    /**
     * The neighbour nodes ID : nodes connected to this node by an edge
     * @return _neighbourNodeId
     */
    std::vector< int >  getNeighbourNodesId ( void ) const ;

    private: //----------------------------------------------------------------

	/*
	 * The (i,j,k) coordinate of this node.
	 */
    std::vector< int > _IJKCoords;

	/*
	 * The index of this node.
	 */
    int _nodeIndex;
    
    /*
     * The coordinate of the Node.
     */
    Point _point ;

    /*
     * The number of cells surrounding this Node.
     */
    int _numberOfCells ;

    /*
     * The number of faces surrounding this Node.
     */
    int _numberOfFaces ;

    /*
     * The number of edges surrounding this Node.
     */
    int _numberOfNodes ;

    /*
     * The indices of cells surrounding this Node.
     */
    const std::vector< int >  _surroundingCells ;

    /*
     * The indices of faces surrounding this Node.
     */
    const std::vector< int >  _surroundingFaces ;

    /*
     * The indices of nodes surrounding this Node.
     */
    const std::vector< int >  _surroundingNodes ;

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

	/*
	 * The mesh MEDCouplingIMesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingIMesh> _mesh;

};

#endif /* IJKNODE_HXX_ */
