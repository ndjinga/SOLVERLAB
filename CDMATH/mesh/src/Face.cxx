/*
 * face.cxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#include "Face.hxx"
#include "CdmathException.hxx"

#include <algorithm> 

using namespace std;

//----------------------------------------------------------------------
Face::Face( void )
//----------------------------------------------------------------------
{
	_measure = 0.0 ;
	_belongToInnerWall = false;
	_isPeriodicFace = false;
	_groupNames=std::vector<std::string>(0);
	_numberOfCells = 0 ;
	_numberOfNodes = 0 ;
	_xN=0.;
	_yN=0.;
	_zN=0.;
}

//----------------------------------------------------------------------
Face::Face( const Face& face )
//----------------------------------------------------------------------
{
	_measure = face.getMeasure() ;
	_belongToInnerWall = face.belongToInnerWall();
	_isPeriodicFace = face.isPeriodicFace();
	_groupNames=face.getGroupNames();
	_point = face.getBarryCenter();
	_numberOfCells = face.getNumberOfCells() ;
	_numberOfNodes = face.getNumberOfNodes() ;
	_nodesId=face.getNodesId();
	_cellsId=face.getCellsId();
	_xN=face.getXN();
	_yN=face.getYN();
	_zN=face.getZN();
}

//----------------------------------------------------------------------
Face::Face( const int numberOfNodes, const int numberOfCells, const double measure, const Point p, double xN, double yN, double zN )
//----------------------------------------------------------------------
{
	_point        = p ;
	_numberOfNodes = numberOfNodes ;
	_numberOfCells   = numberOfCells ;
	_nodesId = std::vector< int >(_numberOfNodes,0);
	_cellsId = std::vector< int >(_numberOfCells,0);
	_measure = measure ;
	_belongToInnerWall = false;
	_isPeriodicFace = false;
	_xN=xN;
	_yN=yN;
	_zN=zN;

    if(numberOfCells==1)
        _groupNames=std::vector<std::string>(1,"Boundary");
    else if(numberOfCells>=2)//On a graph geometry (spaceDim>meshDim), a face is a node and can have more than 2 cell neighbour
        _groupNames=std::vector<std::string>(0);
    else
    {
        cout<<"Error, trying to create a face with "<< numberOfCells <<" cell neighbours"<<endl;
        throw CdmathException("A face must have at least one cell neighbour (MEDCoupling philosophy)");    
    }
}

//----------------------------------------------------------------------
Face::~Face( void )
//----------------------------------------------------------------------
{
}

bool
Face::belongToInnerWall(void) const
{
	return _belongToInnerWall;
}

bool
Face::isPeriodicFace(void) const
{
	return _isPeriodicFace;
}

double
Face::getXN(void) const
{
	return _xN;
}

double
Face::getYN(void) const
{
	return _yN;
}

double
Face::getZN(void) const
{
	return _zN;
}

std::vector<std::string>
Face::getGroupNames(void) const
{
	return _groupNames;
}

string
Face::getGroupName(int igroup) const
{
    if(igroup<_groupNames.size())
        return _groupNames[igroup];
    else
    {
        std::cout<<"Error Face::getGroupName(int igroup), group number "<<igroup+1<<" requested"<<std::endl;
        std::cout<<"Face belongs to "<< _groupNames.size() <<" groups"<<std::endl;
        throw CdmathException("Face has no group with number igroup");
    }
}

void
Face::setGroupName(const string groupName, bool belongToInnerWall)
{
	if(std::find(_groupNames.begin(), _groupNames.end(), groupName) == _groupNames.end())//No group named groupName
	{
		_groupNames.insert(_groupNames.begin(),groupName);
		_belongToInnerWall = belongToInnerWall;
	}
	else
		cout<<"Warning Face::setGroupName, group name "<< groupName <<" is already present. No duplication"<<endl;
}

bool
Face::isBorder(void)
{
	if ( _numberOfCells==1 || _belongToInnerWall)
		return true;
	else
		return false;
}
//----------------------------------------------------------------------
Point
Face::getBarryCenter( void ) const
//----------------------------------------------------------------------
{
	return _point ;
}

//----------------------------------------------------------------------
double
Face::x( void ) const
//----------------------------------------------------------------------
{
  return _point.x() ;
}

//----------------------------------------------------------------------
double
Face::y( void ) const
//----------------------------------------------------------------------
{
	return _point.y() ;
}

//----------------------------------------------------------------------
double
Face::z( void ) const
//----------------------------------------------------------------------
{
	return _point.z() ;
}

//----------------------------------------------------------------------
std::vector< int >
Face::getCellsId( void ) const
//----------------------------------------------------------------------
{
	return _cellsId ;
}

//----------------------------------------------------------------------
std::vector< int >
Face::getNodesId( void ) const
//----------------------------------------------------------------------
{
	return _nodesId ;
}

//----------------------------------------------------------------------
double
Face::getMeasure( void ) const
//----------------------------------------------------------------------
{
	return _measure ;
}

//----------------------------------------------------------------------
int
Face::getNumberOfCells( void ) const
//----------------------------------------------------------------------
{
	return _numberOfCells ;
}

//----------------------------------------------------------------------
int
Face::getNumberOfNodes( void ) const
//----------------------------------------------------------------------
{
	return _numberOfNodes ;
}

//----------------------------------------------------------------------
void
Face::addCellId(const int numCell, const int cellId )
//----------------------------------------------------------------------
{
    if(numCell<_numberOfCells)
       _cellsId[numCell] = cellId ;
    else
    {
        std::cout<< "Cell index : "<< numCell<<" should be strictly smaller than total number of cells "<<_numberOfCells<<std::endl;
        throw CdmathException("Face::addCellId : incorrect cell index");
    }
}

void
Face::addPeriodicCellId ( const int cellId )
{
	int size = _cellsId.size();
    if( _numberOfCells == size ) 
        _cellsId.push_back( cellId );// Increase the size of _cellsId to store the ghost cell index
    else if( _numberOfCells == size - 1 ) 
        _cellsId[size] = cellId;
    else
    {
        std::cout<< "Cell index : "<< cellId<<" can not be added as a ghost cell because face in in wrong state : _numberOfCells = "<<_numberOfCells<<", _cellsId.size() = "<<_cellsId.size()<<std::endl;
        throw CdmathException("Face::addPeriodicCellId : Face is in wrong state");
    }
}

int 
Face::getPeriodicCellId() const 
{
    if( _numberOfCells == _cellsId.size() - 1 ) 
        return _cellsId[_numberOfCells];
    else
    {
        std::cout<< "Face is is not associated to any ghost cell : _numberOfCells = "<<_numberOfCells<<", _cellsId.size() = "<<_cellsId.size()<<std::endl;
        throw CdmathException("Face::getPeriodicCellId : Face is is not associated to any ghost cell");
    }
}

//----------------------------------------------------------------------
void
Face::addNodeId(const int numNode, const int nodeId )
//----------------------------------------------------------------------
{
    if(numNode<_numberOfNodes)
        _nodesId[numNode] = nodeId ;
    else
    {
        std::cout<< "Node index : "<< numNode<<" should be strictly smaller than total number of nodes "<<_numberOfNodes<<std::endl;
        throw CdmathException("Face::addNodeId : incorrect node index");
    }
}

//----------------------------------------------------------------------
const Face&
Face::operator= ( const Face& face )
//----------------------------------------------------------------------
{
	_measure = face.getMeasure() ;
	_point = face.getBarryCenter() ;
	_numberOfCells = face.getNumberOfCells() ;
	_numberOfNodes = face.getNumberOfNodes() ;
	_nodesId=face.getNodesId();
	_cellsId=face.getCellsId();
	_groupNames=face.getGroupNames();
	_belongToInnerWall=face.belongToInnerWall();
	_xN=face.getXN();
	_yN=face.getYN();
	_zN=face.getZN();
	return *this;
}
//------------------------------------------------------------------------
int
Face::getNodeId(int localId) const
//------------------------------------------------------------------------
{
    if(localId<_numberOfNodes)
        return _nodesId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<" total number of nodes= "<<_numberOfNodes<<std::endl;
        throw CdmathException("Face::getNodeId : incorrect node local id");
    }
}
//------------------------------------------------------------------------
int
Face::getCellId(int localId) const
//------------------------------------------------------------------------
{
    if(localId<_numberOfCells)
        return _cellsId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<" total number of cells= "<<_numberOfCells<<std::endl;
        throw CdmathException("Face::getCellId : incorrect cell local id");
    }
}
