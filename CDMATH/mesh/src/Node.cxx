/*
 * Node.cxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#include "Node.hxx"
#include "CdmathException.hxx"
#include <algorithm> 

//----------------------------------------------------------------------
Node::Node( void )
//----------------------------------------------------------------------
{
	_numberOfCells = 0 ;
	_numberOfFaces = 0 ;
	_numberOfEdges = 0 ;
	_groupNames=std::vector<std::string>(0);
	_region=-1;
    _isBorder=false;
}

//----------------------------------------------------------------------
Node::~Node( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
Node::Node( const Node& node )
//----------------------------------------------------------------------
{
	_point = node.getPoint() ;
	_cellsId = node.getCellsId() ;
	_facesId = node.getFacesId() ;
	_numberOfCells = node.getNumberOfCells() ;
	_numberOfFaces = node.getNumberOfFaces() ;
	_numberOfEdges = node.getNumberOfEdges() ;
	_groupNames=node.getGroupNames();
	_region=node.getRegion();
    _isBorder=node.isBorder();
}

//----------------------------------------------------------------------
Node::Node( const int numberOfCells, const int numberOfFaces, const int numberOfEdges, const Point p)
//----------------------------------------------------------------------
{

	_point = p ;
	_numberOfCells = numberOfCells ;
	_numberOfFaces = numberOfFaces ;
	_numberOfEdges = numberOfEdges ;
	_cellsId = std::vector< int >(_numberOfCells,0);
	_facesId = std::vector< int >(_numberOfFaces,0);
	_neighbourNodesId = std::vector< int >(_numberOfEdges,0);
	_groupNames=std::vector<std::string>(0);
	_region=-1;
    _isBorder=false;
}

//----------------------------------------------------------------------
std::vector< int >
Node::getCellsId( void ) const 
//----------------------------------------------------------------------
{
	return _cellsId ;
}

//----------------------------------------------------------------------
std::vector< int >
Node::getFacesId( void ) const 
//----------------------------------------------------------------------
{
	return _facesId ;
}

//----------------------------------------------------------------------
std::vector< int >
Node::getNeighbourNodesId( void ) const 
//----------------------------------------------------------------------
{
	return _neighbourNodesId ;
}

//----------------------------------------------------------------------
int
Node::getNumberOfCells( void ) const 
//----------------------------------------------------------------------
{
	return _numberOfCells ;
}

//----------------------------------------------------------------------
int
Node::getNumberOfFaces( void ) const 
//----------------------------------------------------------------------
{
	return _numberOfFaces ;
}

//----------------------------------------------------------------------
int
Node::getNumberOfEdges( void ) const 
//----------------------------------------------------------------------
{
	return _numberOfEdges ;
}

//----------------------------------------------------------------------
Point
Node::getPoint( void ) const
//----------------------------------------------------------------------
{
	return _point ;
}

//----------------------------------------------------------------------
double
Node::x( void ) const 
//----------------------------------------------------------------------
{
	return _point.x() ;
}

//----------------------------------------------------------------------
double
Node::y( void ) const 
//----------------------------------------------------------------------
{
	return _point.y() ;
}

//----------------------------------------------------------------------
double
Node::z( void ) const 
//----------------------------------------------------------------------
{
	return _point.z() ;
}

std::vector<std::string>
Node::getGroupNames(void) const
{
	return _groupNames;
}

std::string
Node::getGroupName(int igroup) const
{
    if (igroup<_groupNames.size())
        return _groupNames[igroup];
    else
    {
        std::cout<<"Error Node::getGroupName(int igroup), group number "<<igroup+1<<" requested"<<std::endl;
        std::cout<<"Node belongs to "<< _groupNames.size() <<" groups"<<std::endl;
        throw CdmathException("Node has no group with number igroup");
    }
}

void
Node::setGroupName(const std::string groupName)
{
	if(std::find(_groupNames.begin(), _groupNames.end(), groupName) == _groupNames.end())//No group named groupName
	{
		_groupNames.insert(_groupNames.begin(),groupName);
		_region=0;
	}
	else
		std::cout<<"Warning Node::setGroupName, group name "<< groupName <<" is already present. No duplication"<<std::endl;
}

bool
Node::isBorder(void) const
{
	if (_region==0 | _isBorder)
		return true;
	else
		return false;
}

int
Node::getRegion(void) const
{
	return _region;
}

//----------------------------------------------------------------------
void
Node::addFaceId (const int numFace, const int faceId, bool isBorder  )
//----------------------------------------------------------------------
{
	_facesId[numFace] = faceId ;
    if(isBorder)
    {
        if(std::find(_groupNames.begin(), _groupNames.end(), "Boundary") == _groupNames.end())//No group named Boundary
            _groupNames.insert(_groupNames.begin(),"Boundary");
        _isBorder=true;
    }
}

//----------------------------------------------------------------------
void
Node::addCellId (const int numCell, const int cellId )
//----------------------------------------------------------------------
{
	_cellsId[numCell] = cellId ;
}

//----------------------------------------------------------------------
void
Node::addNeighbourNodeId (const int numNode, const int nodeId ) 
//----------------------------------------------------------------------
{
	_neighbourNodesId[numNode] = nodeId ;
}

//----------------------------------------------------------------------
double
Node::distance( const Node& n ) const
//----------------------------------------------------------------------
{
	double distance=_point.distance(n.getPoint());
	return distance ;
}

//----------------------------------------------------------------------
int
//----------------------------------------------------------------------
Node::getFaceId(int localId) const
//----------------------------------------------------------------------
{
    if(localId<_numberOfFaces)
        return _facesId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<", total number of faces= "<<_numberOfFaces<<std::endl;
        throw CdmathException("Node::getFaceId : incorrect face local id");
    }
}

//----------------------------------------------------------------------
int
//----------------------------------------------------------------------
Node::getNeighbourNodeId(int localId) const
//----------------------------------------------------------------------
{
    if(localId<_numberOfEdges)
        return _neighbourNodesId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<", total number of edges= "<<_numberOfEdges<<std::endl;
        throw CdmathException("Node::getNeighbourNodeId : incorrect node local id");
    }
}

const Node& 
Node::operator= ( const Node& node )
//----------------------------------------------------------------------
{
   _point = node.getPoint()  ;
   _cellsId = node.getCellsId() ;
   _facesId = node.getFacesId() ;
   _neighbourNodesId = node.getNeighbourNodesId() ;
   _numberOfCells = node.getNumberOfCells() ;
   _numberOfFaces = node.getNumberOfFaces() ;
   _numberOfEdges = node.getNumberOfEdges() ;
   _groupNames = node.getGroupNames();	
   _isBorder=node.isBorder();
	return *this;
}
