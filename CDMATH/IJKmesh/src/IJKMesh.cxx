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
	
	//_mesh=mesh if _mesh is declared const
	_mesh=mesh->clone(false);//No deep copy : it is assumed node coordinates and cell connectivity will not change
	setFaceMeshes();
}

void 
IJKMesh::setFaceMeshes()
{
	int meshDim = _mesh->getMeshDimension();

	_faceMeshes=std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> >( meshDim )

	std:vector< int > nodeStr = _mesh->getNodeGridStructure();
	std:vector< int > dxyz = _mesh->getDXYZ();
	std:vector< double > origin = _mesh->getOrigin();

	double originPtr[meshDim];
	double dxyzPtr[meshDim];
	mcIdType nodeStrctPtr[meshDim];

	/* Prepare the creation of face meshes, and the filling of face normals, face measures and cell measure */
	for(int i=0; i<meshDim; i++)
	{
		originPtr[i]=origin[i];
		nodeStrctPtr[i]=nodeStr[i];
		dxyzPtr[i]=dxyz[i];
	}
	_cellMeasure=1;
	_faceNormals=std:vector< std:vector< double > > (meshDim, std:vector< double >(meshDim,0));
	/* Creation of face meshes, and filling of face normals, face measures and cell measure */
	for(int i=0; i<meshDim; i++)
	{
		_cellMeasure*=dxyz[i];
		_faceMeasures[i]=1;
		_faceNormals[i][i]=1;
		for(int j=0; j<meshDim; j++)
			if(j != i)
			{
				nodeStrctPtr[j]-=1;
				originPtr[j]+=dxyz[j]/2;
				_faceMeasures[i]*=dxyz[j];
			}
		_faceMesh[i]=MEDCouplingIMesh::New(_mesh->getName()+"_faces_"+std::to_string(i),
				_mesh->getMeshDimension(),
				nodeStrctPtr,
				nodeStrctPtr+meshDim,
				originPtr,
				originPtr+meshDim,
				dxyzPtr,
				dxyzPtr+meshDim);
		for(int j=0; j<meshDim; j++)
			if(j != i)
			{
				nodeStrctPtr[j]+=1;
				originPtr[j]-=dxyz[j]/2;
			}
	}
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
    
	//_mesh=m if _mesh is declared const
	_mesh=m.getMEDCouplingIMesh()->clone(false);
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
	MEDFileCMesh * m;//Here we would like a MEDFileStructuredMesh but that class does not exist
	
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

	setFaceMeshes();
    
	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Structured Mesh name= "<<_mesh->getName()<<", mesh dim="<< _mesh->getMeshDimension()<< ", space dim="<< _mesh->getSpaceDimension()<< ", nb cells= "<<_mesh->getNumberOfCells()<< ", nb nodes= "<<_mesh->getNumberOfNodes()<<endl;

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

	double meshDim = 1 ;
    _epsilon=1e-6;
    _isStructured = true;
    _indexFacePeriodicSet=false;

	double originPtr[meshDim];
	double dxyzPtr[meshDim];
	mcIdType nodeStrctPtr[meshDim];

	originPtr[0]=xmin;
	nodeStrctPtr[0]=nx+1;
	dxyzPtr[0]=dx;

	_mesh=MEDCouplingIIJKMesh::New(meshName,
			meshDim,
			nodeStrctPtr,
			nodeStrctPtr+meshDim,
			originPtr,
			originPtr+meshDim,
			dxyzPtr,
			dxyzPtr+meshDim);
	_measureField = _mesh->getMeasureField(true);    
	setFaceMeshes();
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

	double meshDim = 2 ;
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;

	double originPtr[meshDim];
	double dxyzPtr[meshDim];
	mcIdType nodeStrctPtr[meshDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;

	_mesh=MEDCouplingIMesh::New(meshName,
			meshDim,
			nodeStrctPtr,
			nodeStrctPtr+meshDim,
			originPtr,
			originPtr+meshDim,
			dxyzPtr,
			dxyzPtr+meshDim);
	_measureField = _mesh->getMeasureField(true);    
	setFaceMeshes();
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

	double meshDim = 3;
    _epsilon=1e-6;

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;
	double dz = (zmax - zmin)/nz ;

	double originPtr[meshDim];
	double dxyzPtr[meshDim];
	mcIdType nodeStrctPtr[meshDim];

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
			meshDim,
			nodeStrctPtr,
			nodeStrctPtr+meshDim,
			originPtr,
			originPtr+meshDim,
			dxyzPtr,
			dxyzPtr+meshDim);
	_measureField = _mesh->getMeasureField(true);    
	setFaceMeshes();
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
	return _mesh->getMeshDimension()==2;
}
bool
IJKMesh::isHexahedral() const
{
	return _mesh->getMeshDimension()==3;
}
bool
IJKMesh::isStructured() const
{
	return true;
}

std::string 
IJKMesh::getElementTypes() const
{
	if( _mesh->getMeshDimension()==1 )
		return "Segments ";
	else if( _mesh->getMeshDimension()==2 )
		return "Quadrangles ";
	else if( _mesh->getMeshDimension()==3 )
		return "Hexahedra ";
	else
	{
		cout<< "Mesh " + _name + " does not have acceptable dimension. Mesh dimension is " << meshDim<<endl;
		throw CdmathException("IJKMesh::getElementTypes : wrong dimension");
	}
}

//----------------------------------------------------------------------
double
IJKMesh::getXMin( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_mesh->getMeshDimension()];
    _mesh->getBoundingBox(Box0);

	return Box0[0] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getXMax( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_mesh->getMeshDimension()];
    _mesh->getBoundingBox(Box0);

	return Box0[1] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getYMin( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension()<2)
		throw CdmathException("IJKMesh::getYMin : dimension should be >=2");
		
	double Box0[2*_mesh->getMeshDimension()];
    _mesh->getBoundingBox(Box0);

	return Box0[2] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getYMax( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension()<2)
		throw CdmathException("IJKMesh::getYMax : dimension should be >=2");
		
	double Box0[2*_mesh->getMeshDimension()];
    _mesh->getBoundingBox(Box0);

	return Box0[3] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getZMin( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension()<3)
		throw CdmathException("IJKMesh::getZMin : dimension should be 3");
		
	double Box0[2*_mesh->getMeshDimension()];
    _mesh->getBoundingBox(Box0);

	return Box0[4] ;
}

//----------------------------------------------------------------------
double
IJKMesh::getZMax( void )  const
//----------------------------------------------------------------------
{
	if(_mesh->getMeshDimension()<3)
		throw CdmathException("IJKMesh::getZMax : dimension should be 3");

	double Box0[2*_mesh->getMeshDimension()];
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
	std::vector< vector<int> > result(_mesh->getmeshDimension());
	for(int i = 0; i<_mesh->getmeshDimension(); i++)
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
    
	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	//_mesh=mesh if _mesh is declared const
	_mesh=mesh.getMEDCouplingStructuredMesh()->clone(false);
	
	_faceMeshes=mesh.getFaceMeshes());
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
    return 2*_mesh->getMeshDimension();
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
IJKMesh::writeVTKAllMeshes ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;

	for(int i=0; i< _mesh->getMeshDimension(); i++)
		_faceMeshes[i]->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
IJKMesh::writeMED ( const std::string fileName, bool fromScratch ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	//Save cell mesh after checking mesh is imesh
	const MEDCoupling::MEDCouplingIMesh* iMesh = dynamic_cast< const MEDCoupling::MEDCouplingIMesh* > ((const MEDCoupling::MEDCouplingStructuredMesh*) _mesh);
	if(iMesh)//medcouplingimesh : Use convertToCartesian in order to write mesh
		MEDCoupling::WriteMesh(fname.c_str(),iMesh->convertToCartesian(), fromScratch);
	else//medcouplingcmesh : save directly
		MEDCoupling::WriteMesh(fname.c_str(),_mesh, fromScratch);
}
//----------------------------------------------------------------------
void
IJKMesh::writeMEDAllMeshes ( const std::string fileName, bool fromScratch ) const
//----------------------------------------------------------------------
{
	//Save cell mesh after checking mesh is imesh
	string fname=fileName+".med";
	const MEDCoupling::MEDCouplingIMesh* iMesh = dynamic_cast< const MEDCoupling::MEDCouplingIMesh* > ((const MEDCoupling::MEDCouplingStructuredMesh*) _mesh);
	if(iMesh)//medcouplingimesh : Use convertToCartesian in order to write mesh
		MEDCoupling::WriteMesh(fname.c_str(),iMesh->convertToCartesian(), fromScratch);
	else//medcouplingcmesh : save directly
		MEDCoupling::WriteMesh(fname.c_str(),_mesh, fromScratch);

	//Save face meshes after checking mesh is imesh
	std::vector< const MEDCoupling::MEDCouplingIMesh* > iMeshes(_mesh->getMeshDimension());
	for(int i=0; i< _mesh->getMeshDimension(); i++)
	{
		const MEDCoupling::MEDCouplingIMesh* iMeshes[i] = dynamic_cast< const MEDCoupling::MEDCouplingIMesh* > ((const MEDCoupling::MEDCouplingStructuredMesh*) _faceMeshes[i]);
		if(iMesh)//medcouplingimesh : Use convertToCartesian in order to write mesh
			MEDCoupling::WriteMesh(fname.c_str(),iMeshes[i]->convertToCartesian(), fromScratch);
		else//medcouplingcmesh : save directly
			MEDCoupling::WriteMesh(fname.c_str(), _faceMeshes[i], fromScratch);
	}
}

std::vector< double >   
IJKMesh::getCellCenterCoordinates (mcIdType cellId) const 
{ 
	std::vector< double > result(_mesh->getMeshDimension(),0), coo_node; 
	std::vector< mcIdType > conn=getNodeIdsOfCell(mcIdType cellId); 
	int nbNodes=conn.size();

	for (int i = 0; i< nbNodes; i++)
	{
		coo_node = getNodeCoordinates (conn[i]);
		for(int j=0; j< _mesh->getMeshDimension(); j++)
			result[j]+=coo_node[j];
	}
	for(int j=0; j< _mesh->getMeshDimension(); j++)
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
