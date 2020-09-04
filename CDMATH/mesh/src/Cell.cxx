/*
 * cell.cxx
 *
 *  Created on: 23 janv. 2012
 *      Authors: CDMAT
 */

#include "Cell.hxx"
#include "CdmathException.hxx"

#include <assert.h>

//----------------------------------------------------------------------
Cell::Cell( void )
//----------------------------------------------------------------------
{
	_measure = 0.0 ;
	_numberOfNodes  = 0 ;
	_numberOfFaces    = 0 ;
}

//----------------------------------------------------------------------
Cell::Cell( int numberOfNodes, int numberOfFaces, double measure, const Point p )
//----------------------------------------------------------------------
{
	_point        = p ;
	_numberOfNodes  = numberOfNodes ;
	_numberOfFaces    = numberOfFaces ;
	_nodesId = std::vector< int >(_numberOfNodes,0);
	_facesId = std::vector< int >(_numberOfFaces,0);
	_normalVectors = Vector(3*_numberOfFaces);
	_measure = measure ;
}

//----------------------------------------------------------------------
Cell::~Cell( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
Cell::Cell( const Cell& cell )
//----------------------------------------------------------------------
{
	_point       = cell.getBarryCenter() ;
	_numberOfNodes = cell.getNumberOfNodes();
	_numberOfFaces   = cell.getNumberOfFaces() ;
	_measure = cell.getMeasure() ;
	_nodesId=cell.getNodesId();
	_facesId=cell.getFacesId();
	_normalVectors=cell.getNormalVectors();
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------
std::vector< int >
Cell::getNodesId ( void ) const 
//----------------------------------------------------------------------
{
	return _nodesId ;
}

//----------------------------------------------------------------------
double
Cell::getMeasure ( void ) const
//----------------------------------------------------------------------
{
	return _measure ;
}

//----------------------------------------------------------------------
std::vector< int >
Cell::getFacesId ( void ) const 
//----------------------------------------------------------------------
{
	return _facesId ;
}

//----------------------------------------------------------------------
int
Cell::getNumberOfNodes ( void ) const 
//----------------------------------------------------------------------
{
	return _numberOfNodes ;
}

//----------------------------------------------------------------------
Point
Cell::getBarryCenter( void ) const
//----------------------------------------------------------------------
{
	return _point ;
}

//----------------------------------------------------------------------
double
Cell::x( void ) const 
//----------------------------------------------------------------------
{
  return _point.x() ;
}

//----------------------------------------------------------------------
double
Cell::y( void ) const 
//----------------------------------------------------------------------
{
	return _point.y() ;
}

//----------------------------------------------------------------------
double
Cell::z( void ) const 
//----------------------------------------------------------------------
{
	return _point.z() ;
}

//----------------------------------------------------------------------
int
Cell::getNumberOfFaces ( void ) const 
//----------------------------------------------------------------------
{
	return _numberOfFaces ;
}

//----------------------------------------------------------------------
void
Cell::addFaceId (int numFace, int faceId )
//----------------------------------------------------------------------
{
	_facesId[numFace] = faceId ;
}

//----------------------------------------------------------------------
void
Cell::addNormalVector (int numNormalVector, double x, double y, double z)
//----------------------------------------------------------------------
{
	_normalVectors(3*numNormalVector) = x ;
	_normalVectors(3*numNormalVector+1) = y ;
	_normalVectors(3*numNormalVector+2) = z ;
}

//----------------------------------------------------------------------
double
Cell::getNormalVector( int numNormalVector, int numComposant ) const
//----------------------------------------------------------------------
{
	if (numComposant==0)
		return _normalVectors(3*numNormalVector);
	else if (numComposant==1)
		return _normalVectors(3*numNormalVector+1);
	else if (numComposant==2)
		return _normalVectors(3*numNormalVector+2);
	else
		throw CdmathException("Cell::getNormalVector, numComposant should be 0, 1 or 2");
}

//----------------------------------------------------------------------
Vector
Cell::getNormalVectors (void) const
//----------------------------------------------------------------------
{
		return _normalVectors;
}

//----------------------------------------------------------------------
void
Cell::addNodeId (int numNode, int nodeId )
//----------------------------------------------------------------------
{
	_nodesId[numNode] = nodeId ;
}

//----------------------------------------------------------------------
const Cell&
Cell::operator= ( const Cell& cell )
//----------------------------------------------------------------------
{
	_point       = cell.getBarryCenter() ;
	_numberOfNodes = cell.getNumberOfNodes() ;
	_numberOfFaces   = cell.getNumberOfFaces() ;
	_measure = cell.getMeasure() ;

	_nodesId = cell.getNodesId();
	_facesId = cell.getFacesId();
	_normalVectors=cell.getNormalVectors();
	return *this;
}

//----------------------------------------------------------------------
int 
Cell::getFaceId(int localId) const
//----------------------------------------------------------------------
{
    if(localId<_numberOfFaces)
        return _facesId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<", total number of faces= "<<_numberOfFaces<<std::endl;
        throw CdmathException("Cell::getfaceId : incorrect face local id");
    }
}

//----------------------------------------------------------------------
int 
Cell::getNodeId(int localId) const
//----------------------------------------------------------------------
{
    if(localId<_numberOfNodes)
        return _nodesId[localId];
    else
    {
        std::cout<< "Local id requested : "<< localId<<", total number of nodes= "<<_numberOfNodes<<std::endl;
        throw CdmathException("Cell::getNodeId : incorrect node local id");
    }
}


