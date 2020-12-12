/*
 * IJKmesh.cxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#include "IJKMesh.hxx"
#include "IJKNode.hxx"
#include "IJKCell.hxx"
#include "IJKFace.hxx"

#include "MEDFileCMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingIMesh.hxx"

#include "CdmathException.hxx"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm> 

using namespace MEDCoupling;
using namespace std;

//----------------------------------------------------------------------
Mesh::Mesh( void )
//----------------------------------------------------------------------
{
	_mesh=NULL;
	_spaceDim = 0 ;
	_meshDim  = 0 ;
	_numberOfNodes = 0;
	_numberOfFaces = 0;
	_numberOfCells = 0;
	_xMin=0.;
	_xMax=0.;
	_yMin=0.;
	_yMax=0.;
	_zMin=0.;
	_zMax=0.;
    _nxyz.resize(0);
    _dxyz.resize(0.);
	_faceGroupNames.resize(0);
	_faceGroups.resize(0);
	_nodeGroupNames.resize(0);
	_nodeGroups.resize(0);
    _indexFacePeriodicSet=false;
    _name="";
}

//----------------------------------------------------------------------
Mesh::~Mesh( void )
//----------------------------------------------------------------------
{
}

std::string 
Mesh::getName( void ) const
{
    return _name;
}

Mesh::Mesh( const MEDCoupling::MEDCouplingIMesh* mesh )
{
	_spaceDim=mesh->getSpaceDimension();
	_meshDim=mesh->getMeshDimension();
	_dxyz=mesh->getDXYZ();
	_nxyz=mesh->getCellGridStructure();
	double* Box0=new double[2*_spaceDim];
	mesh->getBoundingBox(Box0);
    _name=mesh->getName();
    _indexFacePeriodicSet=false;
    
	_xMin=Box0[0];
	_xMax=Box0[1];
	if (_spaceDim>=2)
	{
		_yMin=Box0[2];
		_yMax=Box0[3];
	}
	if (_spaceDim>=3)
	{
		_zMin=Box0[4];
		_zMax=Box0[5];
	}

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	for(int i=0;i<_spaceDim;i++)
	{
		originPtr[i]=Box0[2*i];
		nodeStrctPtr[i]=_nxyz[i]+1;
		dxyzPtr[i]=_dxyz[i];
	}
	_mesh=MEDCouplingIMesh::New(_name,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
	delete [] Box0 ;
}

//----------------------------------------------------------------------
Mesh::Mesh( const IJKMesh& m )
//----------------------------------------------------------------------
{
	_spaceDim = m.getSpaceDimension() ;
	_meshDim = m.getMeshDimension() ;
    _name=m.getName();
    _xMax=m.getXMax();
    _yMin=m.getYMin();
    _yMax=m.getYMax();
    _zMin=m.getZMin();
    _zMax=m.getZMax();
    _nxyz = m.getCellGridStructure() ;
    _dxyz=m.getDXYZ();

	_numberOfNodes = m.getNumberOfNodes();
	_numberOfFaces = m.getNumberOfFaces();
	_numberOfCells = m.getNumberOfCells();
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
Mesh::Mesh( const std::string filename, int meshLevel )
//----------------------------------------------------------------------
{
	readMeshMed(filename, meshLevel);
}

//----------------------------------------------------------------------
void
Mesh::readMeshMed( const std::string filename, const int meshLevel)
//----------------------------------------------------------------------
{
	MEDFileCMesh *m=MEDFileCMesh::New(filename.c_str());//reads the first mesh encountered in the file, otherwise call New (const char *fileName, const char *mName, int dt=-1, int it=-1)
	_mesh=m->getMeshAtLevel(meshLevel);
    _mesh->checkConsistencyLight();
	_mesh->setName(_mesh->getName());
	_meshDim=_mesh->getMeshDimension();
	_spaceDim=_mesh->getSpaceDimension();
    _name=_mesh->getName();
    _indexFacePeriodicSet=false;
    MEDCoupling::MEDCouplingIMesh* structuredMesh = dynamic_cast<MEDCoupling::MEDCouplingIMesh*> (_mesh.retn());
    if(structuredMesh)
    {
        _dxyz=structuredMesh->getDXYZ();
        _nxyz=structuredMesh->getCellGridStructure();
        double* Box0=new double[2*_spaceDim];
        structuredMesh->getBoundingBox(Box0);
    
        _xMin=Box0[0];
        _xMax=Box0[1];
        std::cout<<"nx= "<<_nxyz[0];
        if (_spaceDim>=2)
        {
            _yMin=Box0[2];
            _yMax=Box0[3];
            std::cout<<", "<<"ny= "<<_nxyz[1];
        }
        if (_spaceDim>=3)
        {
            _zMin=Box0[4];
            _zMax=Box0[5];
            std::cout<<", "<<"nz= "<<_nxyz[2];
        }
    }
    else
        throw CdmathException("Mesh::readMeshMed med file does not contain a structured MedcouplingIMesh mesh");
    
	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Structured Mesh name= "<<m->getName()<<", mesh dim="<< _meshDim<< ", space dim="<< _spaceDim<< ", nb cells= "<<getNumberOfCells()<< ", nb nodes= "<<getNumberOfNodes()<<endl;

	m->decrRef();
}

void
Mesh::setPeriodicFaces()
{
    _indexFacePeriodicSet=true;    
}

bool
Mesh::isBorderNode(int nodeid) const
{
	return getNode(nodeid).isBorder();
}

bool
Mesh::isTriangular() const
{
	return false;
}

bool
Mesh::isQuadrangular() const
{
	return _meshDim==2;
}
bool
Mesh::isHexahedral() const
{
	return _meshDim==3;
}
bool
Mesh::isStructured() const
{
	return true;
}

std::string 
Mesh::getElementTypes() const
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
		throw CdmathException("Mesh::getElementTypes : wrong dimension");
	}
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, std::string meshName )
//----------------------------------------------------------------------
{
	if(nx<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx) : nx <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx) : xmin >= xmax");

	double dx = (xmax - xmin)/nx ;

	_spaceDim = 1 ;
	_meshDim  = 1 ;
    _name=meshName;
    _indexFacePeriodicSet=false;

	_xMin=xmin;
	_xMax=xmax;
	_yMin=0.;
	_yMax=0.;
	_zMin=0.;
	_zMax=0.;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	originPtr[0]=xmin;
	nodeStrctPtr[0]=nx+1;
	dxyzPtr[0]=dx;

	_mesh=MEDCouplingIMesh::New(meshName,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;

	_numberOfCells = _mesh->getNumberOfCells() ;

	_numberOfNodes = _mesh->getNumberOfNodes() ;

	_numberOfFaces = _numberOfNodes;
    
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : nx <= 0 or ny <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : ymin >= ymax");

	_xMin=xmin;
	_xMax=xmax;
	_yMin=ymin;
	_yMax=ymax;
	_zMin=0.;
	_zMax=0.;


	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;

	_spaceDim = 2 ;
	_meshDim  = 2 ;
    _name=meshName;
    _indexFacePeriodicSet=false;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_dxyz[1]=dy;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;

	_mesh=MEDCouplingIMesh::New(meshName,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);

	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;

	_numberOfCells = _mesh->getNumberOfCells() ;

	_numberOfNodes = _mesh->getNumberOfNodes() ;

	_numberOfFaces = nx*(ny+1)+ny*(nx+1);
    
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0 || nz<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : nx <= 0 or ny <= 0 or nz <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : ymin >= ymax");
	if(zmin>=zmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : zmin >= zmax");

	_spaceDim = 3;
	_meshDim  = 3;
    _name=meshName;
    _indexFacePeriodicSet=false;
	_xMin=xmin;
	_xMax=xmax;
	_yMin=ymin;
	_yMax=ymax;
	_zMin=zmin;
	_zMax=zmax;

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;
	double dz = (zmax - zmin)/nz ;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_dxyz[1]=dy;
	_dxyz[2]=dz;

	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;
	_nxyz[2]=nz;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

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
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);

	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;

	_numberOfCells = _mesh->getNumberOfCells() ;

	_numberOfNodes = _mesh->getNumberOfNodes() ;

	_numberOfFaces = nx*ny*(nz+1)+nx*nz*(ny+1)+ny*nz*(nx+1);
	
}

//----------------------------------------------------------------------
int
Mesh::getSpaceDimension( void )  const
//----------------------------------------------------------------------
{
	return _spaceDim ;
}

//----------------------------------------------------------------------
int
Mesh::getMeshDimension( void )  const
//----------------------------------------------------------------------
{
	return _meshDim ;
}

std::vector<double>
Mesh::getDXYZ() const
{
	return _dxyz;
}

std::vector<int>
Mesh::getCellGridStructure() const
{
	return _nxyz;
}

//----------------------------------------------------------------------
int
Mesh::getNx( void )  const
//----------------------------------------------------------------------
{
	return _nxyz[0];
}

//----------------------------------------------------------------------
int
Mesh::getNy( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim < 2)
		throw CdmathException("int double& Field::operator[ielem] : Ny is not defined in dimension < 2!");

	return _nxyz[1];
}

//----------------------------------------------------------------------
int
Mesh::getNz( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim < 3)
		throw CdmathException("int Mesh::getNz( void ) : Nz is not defined in dimension < 3!");

	return _nxyz[2];
}

//----------------------------------------------------------------------
double
Mesh::getXMin( void )  const
//----------------------------------------------------------------------
{        
	return _xMin ;
}

//----------------------------------------------------------------------
double
Mesh::getXMax( void )  const
//----------------------------------------------------------------------
{
	return _xMax ;
}

//----------------------------------------------------------------------
double
Mesh::getYMin( void )  const
//----------------------------------------------------------------------
{
	return _yMin ;
}

//----------------------------------------------------------------------
double
Mesh::getYMax( void )  const
//----------------------------------------------------------------------
{
	return _yMax ;
}

//----------------------------------------------------------------------
double
Mesh::getZMin( void )  const
//----------------------------------------------------------------------
{
	return _zMin ;
}

//----------------------------------------------------------------------
double
Mesh::getZMax( void )  const
//----------------------------------------------------------------------
{
	return _zMax ;
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingIMesh>
Mesh::getMEDCouplingIMesh( void )  const
//----------------------------------------------------------------------
{
	return _mesh ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfNodes ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfNodes ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfCells ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfCells ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfFaces ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfFaces ;
}

//----------------------------------------------------------------------
Cell&
Mesh::getCell ( int i ) const
//----------------------------------------------------------------------
{
	return _cells[i] ;
}

//----------------------------------------------------------------------
Face&
Mesh::getFace ( int i ) const
//----------------------------------------------------------------------
{
	return _faces[i] ;
}

//----------------------------------------------------------------------
Node&
Mesh::getNode ( int i ) const
//----------------------------------------------------------------------
{
	return _nodes[i] ;
}

vector<string>
Mesh::getNameOfFaceGroups( void )  const
{
	return _faceGroupNames;
}

vector<MEDCoupling::MEDCouplingIMesh *>
Mesh::getFaceGroups( void )  const
{
	return _faceGroups;
}

vector<string>
Mesh::getNameOfNodeGroups( void )  const
{
	return _nodeGroupNames;
}

vector<MEDCoupling::DataArrayIdType *>
Mesh::getNodeGroups( void )  const
{
	return _nodeGroups;
}

//----------------------------------------------------------------------
const IJKMesh&
Mesh::operator= ( const IJKMesh& mesh )
//----------------------------------------------------------------------
{
	_spaceDim = mesh.getSpaceDimension() ;
	_meshDim  = mesh.getMeshDimension() ;
    _name = mesh.getName();
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();
    
        _nxyz = mesh.getCellGridStructure() ;
        _dxyz=mesh.getDXYZ();
        _xMin=mesh.getXMin();
        _xMax=mesh.getXMax();
        _yMin=mesh.getYMin();
        _yMax=mesh.getYMax();
        _zMin=mesh.getZMin();
        _zMax=mesh.getZMax();

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	MCAuto<MEDCouplingIMesh> m1=mesh.getMEDCouplingIMesh()->deepCopy();
	_mesh=m1;
	return *this;
}

bool Mesh::isIndexFacePeriodicSet() const
{
 return    _indexFacePeriodicSet;
}
//----------------------------------------------------------------------
double 
Mesh::minRatioVolSurf()
{
    return dx_min;
}
int 
Mesh::getMaxNbNeighbours(EntityType type) const
{
    return 2*_meshDim;
}
//----------------------------------------------------------------------
void
Mesh::writeVTK ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
Mesh::writeMED ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	MEDCoupling::WriteCMesh(fname.c_str(),_mesh,true);

	mu->decrRef();
}
