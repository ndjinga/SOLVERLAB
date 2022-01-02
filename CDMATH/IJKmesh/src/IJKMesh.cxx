/*
 * IJKmesh.cxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm> 

#include "IJKMesh.hxx"
#include "IJKNode.hxx"
#include "IJKCell.hxx"
#include "IJKFace.hxx"

#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDFileCMesh.hxx"
#include "MEDLoader.hxx"

#include "CdmathException.hxx"

using namespace MEDCoupling;
using namespace std;

//----------------------------------------------------------------------
IJKMesh::IJKMesh( void )
//----------------------------------------------------------------------
{
	_mesh=NULL;
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(0)
	_measureField=NULL;
	_faceGroupNames.resize(0);
	_faceGroups.resize(0);
	_nodeGroupNames.resize(0);
	_nodeGroups.resize(0);
    _indexFacePeriodicSet=false;
    _name="";
    _epsilon=1e-6;
}

//----------------------------------------------------------------------
IJKMesh::~IJKMesh( void )
//----------------------------------------------------------------------
{
	_measureField->decrRef();
}

std::string 
IJKMesh::getName( void ) const
{
    return _name;
}

IJKMesh::IJKMesh( const MEDCoupling::MEDCouplingStructuredMesh* mesh )
{
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
	_measureField = mesh->getMeasureField(true);    
	
	_mesh=mesh->clone(false);//No deep copy : it is assumed node coordinates and cell connectivity will not change
}

//----------------------------------------------------------------------
IJKMesh::IJKMesh( const IJKMesh& m )
//----------------------------------------------------------------------
{
	_measureField = m->getMeasureField(true);    
	_faceMeshes=m.getFaceMeshes();
	_faceGroupNames = m.getNameOfFaceGroups() ;
	_faceGroups = m.getFaceGroups() ;
	_nodeGroupNames = m.getNameOfNodeGroups() ;
	_nodeGroups = m.getNodeGroups() ;
    _indexFacePeriodicSet= m.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=m.getIndexFacePeriodic();
    
	MCAuto<MEDCouplingIMesh> m1=m.getMEDCouplingIMesh()->deepCopy();
	_mesh=m1;
}

//----------------------------------------------------------------------
IJKMesh::IJKMesh( const std::string filename, const std::string & meshName, int meshLevel )
//----------------------------------------------------------------------
{
	readMeshMed(filename, meshName, meshLevel);
}

//----------------------------------------------------------------------
void
IJKMesh::readMeshMed( const std::string filename, const std::string & meshName , int meshLevel )
//----------------------------------------------------------------------
{
	MEDFileCMesh *m;//Here we would like a MEDFileStructuredMesh but that class does not exist
	
	if( meshName == "" )
		m=MEDFileCMesh::New(filename.c_str());//reads the first mesh encountered in the file, otherwise call New (const char *fileName, const char *mName, int dt=-1, int it=-1)
	else
		m=MEDFileCMesh::New(filename.c_str(), meshName.c_str());//seeks the mesh named meshName in the file
		
	_mesh=m->getMeshAtLevel(meshLevel);
    _mesh->checkConsistencyLight();
	_mesh->setName(_mesh->getName());
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
	_measureField = mesh->getMeasureField(true);    
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(_mesh->getMeshDimension())
    
	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Structured Mesh name= "<<m->getName()<<", mesh dim="<< _meshDim<< ", space dim="<< _spaceDim<< ", nb cells= "<<getNumberOfCells()<< ", nb nodes= "<<getNumberOfNodes()<<endl;

	m->decrRef();
}

//----------------------------------------------------------------------
IJKMesh::IJKMesh( double xmin, double xmax, int nx, std::string meshName )
//----------------------------------------------------------------------
{
	if(nx<=0)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx) : nx <= 0");
	if(xmin>=xmax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx) : xmin >= xmax");

	double dx = (xmax - xmin)/nx ;

	double spaceDim = 1 ;
    _epsilon=1e-6;
    _isStructured = true;
    _indexFacePeriodicSet=false;

	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;

	double originPtr[_spaceDim];
	double dxyzPtr[_spaceDim];
	mcIdType nodeStrctPtr[_spaceDim];

	originPtr[0]=xmin;
	nodeStrctPtr[0]=nx+1;
	dxyzPtr[0]=dx;

	_mesh=MEDCouplingIIJKMesh::New(meshName,
			spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+spaceDim,
			originPtr,
			originPtr+spaceDim,
			dxyzPtr,
			dxyzPtr+spaceDim);
	_measureField = _mesh->getMeasureField(true);    
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(spaceDim)
}

//----------------------------------------------------------------------
IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : nx <= 0 or ny <= 0");
	if(xmin>=xmax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : ymin >= ymax");

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;

	double spaceDim = 2 ;
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;

	double originPtr[_spaceDim];
	double dxyzPtr[_spaceDim];
	mcIdType nodeStrctPtr[_spaceDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;

	_mesh=MEDCouplingIMesh::New(meshName,
			spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+spaceDim,
			originPtr,
			originPtr+spaceDim,
			dxyzPtr,
			dxyzPtr+spaceDim);
	_measureField = _mesh->getMeasureField(true);    
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(spaceDim)
}

//----------------------------------------------------------------------
IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0 || nz<=0)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : nx <= 0 or ny <= 0 or nz <= 0");
	if(xmin>=xmax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : ymin >= ymax");
	if(zmin>=zmax)
		throw CdmathException("IJKMesh::IJKMesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : zmin >= zmax");

	double spaceDim = 3;
    _epsilon=1e-6;

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;
	double dz = (zmax - zmin)/nz ;

	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;
	_nxyz[2]=nz;

	double originPtr[_spaceDim];
	double dxyzPtr[_spaceDim];
	mcIdType nodeStrctPtr[_spaceDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	originPtr[2]=zmin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	nodeStrctPtr[2]=nz+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;
	dxyzPtr[2]=dz;

	_mesh=MEDCouplingIMesh::New(meshName,
			spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+spaceDim,
			originPtr,
			originPtr+spaceDim,
			dxyzPtr,
			dxyzPtr+spaceDim);
	_measureField = _mesh->getMeasureField(true);    
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(spaceDim)
}

void
IJKMesh::setPeriodicFaces()
{
    _indexFacePeriodicSet=true;    
}

bool
IJKMesh::isBorderNode(int nodeid) const
{
	return getNode(nodeid).isBorder();
}

bool
IJKMesh::isTriangular() const
{
	return false;
}

bool
IJKMesh::isQuadrangular() const
{
	return _meshDim==2;
}
bool
IJKMesh::isHexahedral() const
{
	return _meshDim==3;
}
bool
IJKMesh::isStructured() const
{
	return true;
}

std::string 
IJKMesh::getElementTypes() const
{
	if( _meshDim==1 )
		return "Segments ";
	else if( _meshDim==2 )
		return "Quadrangles ";
	else if( _meshDim==3 )
		return "Hexahedra ";
	else
	{
		cout<< "Mesh " + _name + " does not have acceptable dimension. Dimension is " << _meshDim<<endl;
		throw CdmathException("IJKMesh::getElementTypes : wrong dimension");
	}
}

//----------------------------------------------------------------------
double
IJKMesh::getXMin( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[0] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getXMax( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[1] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getYMin( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim<2)
		throw CdmathException("IJKMesh::getYMin : dimension should be >=2");
		
	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[2] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getYMax( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim<2)
		throw CdmathException("IJKMesh::getYMax : dimension should be >=2");
		
	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[3] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getZMin( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim<3)
		throw CdmathException("IJKMesh::getZMin : dimension should be 3");
		
	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[4] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getZMax( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim<3)
		throw CdmathException("IJKMesh::getZMax : dimension should be 3");

	double Box0[2*_spaceDim];
    _mesh->getBoundingBox(Box0);

	return Box0[5] ;
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingStructuredMesh>
IJKMesh::getMEDCouplingStructuredMesh( void )  const
//----------------------------------------------------------------------
{
	return _mesh ;
}

std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> > 
IJKMesh::getFaceMeshes() const
{
	return _faceMeshes;
}	
//----------------------------------------------------------------------
int
IJKMesh::getNumberOfNodes ( void ) const
//----------------------------------------------------------------------
{
	return _mesh->getNumberOfNodes ();
}

//----------------------------------------------------------------------
int
IJKMesh::getNumberOfCells ( void ) const
//----------------------------------------------------------------------
{
	return _mesh->getNumberOfCells (); 
}

std::vector<int> 
IJKMesh::getCellGridStructure() const
{
	return _mesh->getCellGridStructure();
}

std::vector<int> 
IJKMesh::getNodeGridStructure() const
{
	return _mesh->getNodeGridStructure();
}

std::vector< vector<int> >
IJKMesh::getFaceGridStructures() const
{
	std::vector< vector<int> > result(_mesh->getSpaceDimension());
	for(int i = 0; i<_mesh->getSpaceDimension(); i++)
		result[i] = _faceMeshes[i]->getNodeGridStructure();
		
	return result;
}

//----------------------------------------------------------------------
int
IJKMesh::getNumberOfFaces ( void ) const
//----------------------------------------------------------------------
{
	std::vector<int> NxNyNz = _mesh->getCellGridStructure();
	
	switch( _mesh->getMeshDimension () )
	{
		case 1:
			return NxNyNz[0]+1;
		case 2:
			return NxNyNz[0]*(1+NxNyNz[1]) + NxNyNz[1]*(1 + NxNyNz[0]);
		case 2:
			return NxNyNz[0]*NxNyNz[1]*(1+NxNyNz[2]) + NxNyNz[0]*NxNyNz[2]*(1+NxNyNz[1]) + NxNyNz[1]*NxNyNz[2]*(1+NxNyNz[0]);
		default
			throw CdmathException("IJKMesh::getNumberOfFaces space dimension must be between 1 and 3");
	}
}

int 
IJKMesh::getNumberOfEdges ( void )  const 
{
	std::vector<int> NxNyNz = _mesh->getCellGridStructure();
	
	switch( _mesh->getMeshDimension () )
	{
		case 1:
			return NxNyNz[0];
		case 2:
			return NxNyNz[0]*(1+NxNyNz[1]) + NxNyNz[1]*(1 + NxNyNz[0]);
		case 2:
			return NxNyNz[0]*(1+NxNyNz[1])*(1+NxNyNz[2]) + NxNyNz[1]*(1 + NxNyNz[0])*(1+NxNyNz[2]) + NxNyNz[2]*(1 + NxNyNz[0])*(1+NxNyNz[1]);
		default
			throw CdmathException("IJKMesh::getNumberOfEdges space dimension must be between 1 and 3");
	}
};

vector<string>
IJKMesh::getNameOfFaceGroups( void )  const
{
	return _faceGroupNames;
}

vector<MEDCoupling::MEDCouplingIMesh *>
IJKMesh::getFaceGroups( void )  const
{
	return _faceGroups;
}

vector<string>
IJKMesh::getNameOfNodeGroups( void )  const
{
	return _nodeGroupNames;
}

vector<MEDCoupling::DataArrayIdType *>
IJKMesh::getNodeGroups( void )  const
{
	return _nodeGroups;
}

//----------------------------------------------------------------------
const IJKMesh&
IJKMesh::operator= ( const IJKMesh& mesh )
//----------------------------------------------------------------------
{
    _epsilon=mesh.getComparisonEpsilon();
    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();
    
    _nxyz = mesh.getCellGridStructure() ;

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	_mesh=mesh.getMEDCouplingStructuredMesh()->clone(false);
	
	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >(_mesh->getMeshDimension());
	return *this;
}
 
bool IJKMesh::isIndexFacePeriodicSet() const
{
 return    _indexFacePeriodicSet;
}
//----------------------------------------------------------------------
double 
IJKMesh::minRatioVolSurf()
{
    return dx_min;
}
int 
IJKMesh::getMaxNbNeighbours(EntityType type) const
{
    return 2*_meshDim;
}
//----------------------------------------------------------------------
void
IJKMesh::writeVTK ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
IJKMesh::writeMED ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	MEDCoupling::WriteCMesh(fname.c_str(),_mesh,true);

	mu->decrRef();
}

std::vector< double >   
IJKMesh::getCellCenterCoordinates (mcIdType cellId) const 
{ 
	std::vector< double > result(_spaceDim,0), coo_node; 
	std::vector< mcIdType > conn=getNodeIdsOfCell(mcIdType cellId); 
	int nbNodes=conn.size();

	for (int i = 0; i< nbNodes; i++)
	{
		coo_node = getNodeCoordinates (conn[i]);
		for(int j=0; j< _spaceDim; j++)
			result[j]+=coo_node[j];
	}
	for(int j=0; j< _spaceDim; j++)
		result[j]/=nbNodes;
	return result; 
};

//----------------------------------------------------------------------
int
IJKMesh::getNx( void )  const
//----------------------------------------------------------------------
{
	return _mesh->getCellGridStructure()[0];
}

//----------------------------------------------------------------------
int
IJKMesh::getNy( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension () < 2)
		throw CdmathException("int IJKMesh::getNy( void ) : Ny is not defined in dimension < 2!");
	else
		return _mesh->getCellGridStructure()[1];
}

//----------------------------------------------------------------------
int
IJKMesh::getNz( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension () < 3)
		throw CdmathException("int IJKMesh::getNz( void ) : Nz is not defined in dimension < 3!");
	else
		return _mesh->getCellGridStructure()[2];
}
