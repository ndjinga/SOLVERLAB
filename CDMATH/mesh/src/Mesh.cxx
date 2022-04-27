/*
 * mesh.cxx
 *
 *  Created on: 22 janv. 2012
 *      Authors: CDMATH
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm> 
#include <cassert> 

#include "Mesh.hxx"
#include "Node.hxx"
#include "Cell.hxx"
#include "Face.hxx"

#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"

#include "CdmathException.hxx"

using namespace MEDCoupling;
using namespace std;

//----------------------------------------------------------------------
Mesh::Mesh( void )
//----------------------------------------------------------------------
{
	_mesh=NULL;
    _meshNotDeleted=false;
	_cells=NULL;
	_nodes=NULL;
	_faces=NULL;
	_spaceDim = 0 ;
	_meshDim  = 0 ;
	_numberOfNodes = 0;
	_numberOfFaces = 0;
	_numberOfCells = 0;
	_numberOfEdges = 0;
    _isStructured=false;
    _nxyz.resize(0);
	_boundaryMesh=NULL;
    _boundaryFaceIds.resize(0);
    _boundaryNodeIds.resize(0);
	_faceGroupNames.resize(0);
	_faceGroups.resize(0);
	_faceGroupsIds.resize(0);
	_nodeGroupNames.resize(0);
	_nodeGroups.resize(0);
	_nodeGroupsIds.resize(0);
    _indexFacePeriodicSet=false;
    _name="";
    _epsilon=1e-6;
}

//----------------------------------------------------------------------
Mesh::~Mesh( void )
//----------------------------------------------------------------------
{
	//for(int i=0; i< _faceGroups.size(); i++)
	//	_faceGroups[i]->decrRef();
	//	_nodeGroups[i]->decrRef();
	if( _meshNotDeleted)
		(_mesh.retn())->decrRef();
	//if( _boundaryMesh)
	//    _boundaryMesh->decrRef();
}

std::string 
Mesh::getName( void ) const
{
    return _name;
}

Mesh::Mesh( MEDCoupling::MCAuto<const MEDCoupling::MEDCouplingMesh> mesh )
{
	_spaceDim=mesh->getSpaceDimension();
	_meshDim=mesh->getMeshDimension();
    _name=mesh->getName();
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
	_meshNotDeleted=true;
    
	_mesh= mesh->clone(false);//Clone because you will need to buildUnstructured. No deep copy : it is assumed node coordinates and cell connectivity will not change

    MEDCoupling::MEDCouplingStructuredMesh* structuredMesh = dynamic_cast<MEDCoupling::MEDCouplingStructuredMesh*> (_mesh.retn());
    if(structuredMesh)
    {
        _isStructured=true;
        _nxyz=structuredMesh->getCellGridStructure();
    }
    else
        _isStructured=false;

	MEDCouplingUMesh*  mu = setMesh();
	mu->decrRef();
}

//----------------------------------------------------------------------
Mesh::Mesh( const std::string filename, const std::string & meshName, int meshLevel)
//----------------------------------------------------------------------
{
	readMeshMed(filename, meshName, meshLevel);
}

//----------------------------------------------------------------------
Mesh::Mesh( const Mesh& mesh )
//----------------------------------------------------------------------
{
	_spaceDim = mesh.getSpaceDimension() ;
	_meshDim = mesh.getMeshDimension() ;
    _name=mesh.getName();
    _epsilon=mesh.getComparisonEpsilon();
    _isStructured=mesh.isStructured();
    if(_isStructured)
        _nxyz = mesh.getCellGridStructure() ;
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
	_numberOfEdges = mesh.getNumberOfEdges();

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getMEDCouplingFaceGroups() ;
	_faceGroupsIds = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getMEDCouplingNodeGroups() ;
	_nodeGroupsIds = mesh.getNodeGroups() ;

	_nodes   = mesh.getNodes() ;
	_faces   = mesh.getFaces() ;
	_cells   = mesh.getCells() ;
    
    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();

    _boundaryFaceIds=mesh.getBoundaryFaceIds();
    _boundaryNodeIds=mesh.getBoundaryNodeIds();
	
    _boundaryMesh=mesh.getBoundaryMEDCouplingMesh();
		
    _eltsTypes=mesh.getElementTypes();
    _eltsTypesNames=mesh.getElementTypesNames();
    
	MCAuto<MEDCouplingMesh> m1=mesh.getMEDCouplingMesh()->clone(false);//Clone because you will need to buildUnstructured. No deep copy : it is assumed node coordinates and cell connectivity will not change

	_mesh=m1;
    _meshNotDeleted=mesh.meshNotDeleted();
}

//----------------------------------------------------------------------
void
Mesh::readMeshMed( const std::string filename, const std::string & meshName, int meshLevel)
//----------------------------------------------------------------------
{
	MEDFileMesh *m;
	if( meshName == "" )
		m=MEDFileMesh::New(filename.c_str());//reads the first mesh encountered in the file, otherwise call New (const char *fileName, const char *mName, int dt=-1, int it=-1)
	else
		m=MEDFileMesh::New(filename.c_str(), meshName.c_str());//seeks the mesh named meshName in the file

	_mesh=m->getMeshAtLevel(meshLevel);
    _mesh->checkConsistencyLight();
	_mesh->setName(_mesh->getName());
	_meshDim=_mesh->getMeshDimension();
	_spaceDim=_mesh->getSpaceDimension();
    _name=_mesh->getName();
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
	_meshNotDeleted=true;

    MEDCoupling::MEDCouplingStructuredMesh* structuredMesh = dynamic_cast<MEDCoupling::MEDCouplingStructuredMesh*> (_mesh.retn());
    if(structuredMesh)
    {
        _isStructured=true;
        _nxyz=structuredMesh->getCellGridStructure();
    }
    else
        _isStructured=false;
    
	MEDCouplingUMesh*  mu = setMesh();
	setNodeGroups(m, mu);//Works for both cartesan and unstructured meshes
	MEDFileUMesh *umedfile=dynamic_cast< MEDFileUMesh * > (m);
	if(umedfile)
		setFaceGroups(umedfile, mu);//Works only for unstructured meshes

	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Mesh name = "<<m->getName()<<", mesh dim = "<< _meshDim<< ", space dim = "<< _spaceDim<< ", nb cells= "<<getNumberOfCells()<< ", nb nodes= "<<getNumberOfNodes()<<endl;

	m->decrRef();
	mu->decrRef();
}

//----------------------------------------------------------------------
Mesh::Mesh( std::vector<double> points, std::string meshName )
//----------------------------------------------------------------------
{
    int nx=points.size();
    
	if(nx<2)
		throw CdmathException("Mesh::Mesh( vector<double> points, string meshName) : nx < 2, vector should contain at least two values");
    int i=0;
    while( i<nx-1 && points[i+1]>points[i] )
        i++;
	if( i!=nx-1 )
    {
        //cout<< points << endl;
		throw CdmathException("Mesh::Mesh( vector<double> points, string meshName) : vector values should be sorted");
    }
    
	_spaceDim = 1 ;
	_meshDim  = 1 ;
    _name=meshName;
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
    
    MEDCouplingUMesh * mesh1d = MEDCouplingUMesh::New(meshName, 1);
    mesh1d->allocateCells(nx - 1);
    double * coords = new double[nx];
    mcIdType nodal_con[2];
    coords[0]=points[0];
    for(int i=0; i<nx- 1 ; i++)
    {
        nodal_con[0]=i;
        nodal_con[1]=i+1;
        mesh1d->insertNextCell(INTERP_KERNEL::NORM_SEG2, 2, nodal_con);
        coords[i+1]=points[i + 1];
    }
    mesh1d->finishInsertingCells();

    DataArrayDouble * coords_arr = DataArrayDouble::New();
    coords_arr->alloc(nx,1);
    std::copy(coords,coords+nx,coords_arr->getPointer());
    mesh1d->setCoords(coords_arr);

    delete [] coords;
    coords_arr->decrRef();

    _mesh=mesh1d;//To enable writeMED. Because we declared the mesh as unstructured, we decide to build the unstructured data (not mandatory)
    _meshNotDeleted=true;
    _isStructured = false;

	MEDCouplingUMesh*  mu = setMesh();
	mu->decrRef();
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
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;

	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;

	double originPtr[_spaceDim];
	double dxyzPtr[_spaceDim];
	mcIdType nodeStrctPtr[_spaceDim];

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
    _meshNotDeleted=true;
    _isStructured = true;

	MEDCouplingUMesh*  mu = setMesh();
	mu->decrRef();
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, int split_to_triangles_policy, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : nx <= 0 or ny <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : ymin >= ymax");

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;

	_spaceDim = 2 ;
	_meshDim  = 2 ;
    _name=meshName;
    _epsilon=1e-6;
    _indexFacePeriodicSet=false;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;

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
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
    _meshNotDeleted=true;
    _isStructured = true;

    if(split_to_triangles_policy==0 || split_to_triangles_policy==1)
        {
            _mesh=_mesh->buildUnstructured();//simplexize is not available for structured meshes
            DataArrayIdType * o2n = _mesh->simplexize(split_to_triangles_policy);
            o2n->decrRef();
			_isStructured = false;
        }
    else if (split_to_triangles_policy != -1)
        {
            cout<< "split_to_triangles_policy = "<< split_to_triangles_policy << endl;
            throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : Unknown splitting policy");
        }
    
	MEDCouplingUMesh*  mu = setMesh();
	mu->decrRef();
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, int split_to_tetrahedra_policy, std::string meshName)
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
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
    _meshNotDeleted=true;
    _isStructured = true;

    if( split_to_tetrahedra_policy == 0 )
        {
            _mesh=_mesh->buildUnstructured();//simplexize is not available for structured meshes
            DataArrayIdType * o2n = _mesh->simplexize(INTERP_KERNEL::PLANAR_FACE_5);
            o2n->decrRef();
			_isStructured = false;
        }
    else if( split_to_tetrahedra_policy == 1 )
        {
            _mesh=_mesh->buildUnstructured();//simplexize is not available for structured meshes
            DataArrayIdType * o2n = _mesh->simplexize(INTERP_KERNEL::PLANAR_FACE_6);
            o2n->decrRef();
			_isStructured = false;
        }
    else if ( split_to_tetrahedra_policy != -1 )
        {
            cout<< "split_to_tetrahedra_policy = "<< split_to_tetrahedra_policy << endl;
            throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : splitting policy value should be 0 or 1");
        }
    
	MEDCouplingUMesh*  mu = setMesh();
	mu->decrRef();
}

//----------------------------------------------------------------------
MEDCouplingUMesh* 
Mesh::setMesh( void )
//----------------------------------------------------------------------
{
	/* This is the main function translating medcouplingumesh info into Mesh class to be used when designing numerical methods
	 * We need the level 0 mesh to extract the cell-node connectvity
	 * We need the level -1 mesh to extract the cell-face and face-node connectivities (use o build descending connectivity)
	 * Be careful : the nodes in the medcoupling mesh are not necessarily all conected to a cell/face. 
	 * Mesh class discard isolated nodes, hence the number of nodes in Mesh class can be lower than the number of nodes in medcouplingumesh.
	 */
	 
	DataArrayIdType *desc  = DataArrayIdType::New();
	DataArrayIdType *descI = DataArrayIdType::New();
	DataArrayIdType *revDesc  = DataArrayIdType::New();
	DataArrayIdType *revDescI = DataArrayIdType::New();
	MEDCouplingUMesh* mu = _mesh->buildUnstructured();
	MEDCouplingUMesh* mu2;//mesh of dimension N-1 containing the cell interfaces->cell/face connectivity
	
	mu->unPolyze();
    DataArrayIdType * o2n = mu->sortCellsInMEDFileFrmt( );
	o2n->decrRef();
	
	if(_meshDim<2)
		mu2=mu->buildDescendingConnectivity(desc,descI,revDesc,revDescI);
	else
		mu2=mu->buildDescendingConnectivity2(desc,descI,revDesc,revDescI);
	
    const mcIdType *tmp = desc->getConstPointer();//Lists the faces surrounding each cell
    const mcIdType *tmpI=descI->getConstPointer();

	const mcIdType *tmpA =revDesc->getConstPointer();//Lists the cells surrounding each face
	const mcIdType *tmpAI=revDescI->getConstPointer();

	//Test du type d'éléments contenus dans le maillage afin d'éviter les éléments contenant des points de gauss
	_eltsTypes=mu->getAllGeoTypesSorted();
	for(int i=0; i<_eltsTypes.size();i++)
	{
		if(
				   _eltsTypes[i]!= INTERP_KERNEL::NORM_POINT1 && _eltsTypes[i]!= INTERP_KERNEL::NORM_SEG2
				&& _eltsTypes[i]!= INTERP_KERNEL::NORM_TRI3   && _eltsTypes[i]!= INTERP_KERNEL::NORM_QUAD4
				&& _eltsTypes[i]!= INTERP_KERNEL::NORM_TETRA4 && _eltsTypes[i]!= INTERP_KERNEL::NORM_PYRA5
				&& _eltsTypes[i]!= INTERP_KERNEL::NORM_PENTA6 && _eltsTypes[i]!= INTERP_KERNEL::NORM_HEXA8
				&& _eltsTypes[i]!= INTERP_KERNEL::NORM_POLYGON&& _eltsTypes[i]!= INTERP_KERNEL::NORM_POLYHED
		)
		{
			cout<< "Mesh " + mu->getName() + " contains an element of type " <<endl;
			cout<< _eltsTypes[i]<<endl;
			throw CdmathException("Mesh::setMesh : in order to avoid gauss points, mesh should contain elements of type NORM_POINT1, NORM_SEG2, NORM_TRI3, NORM_QUAD4, NORM_TETRA4, NORM_PYRA5, NORM_PENTA6, NORM_HEXA8, NORM_POLYGON, NORM_POLYHED");
		}
	}

	DataArrayDouble *baryCell = mu->computeCellCenterOfMass() ;
	const double *coorBary=baryCell->getConstPointer();//Used for cell center coordinates

	MEDCouplingFieldDouble* fields=mu->getMeasureField(true);
	DataArrayDouble *surface = fields->getArray();
	const double *surf=surface->getConstPointer();//Used for cell lenght/surface/volume

	DataArrayDouble *coo = mu->getCoords() ;
	const double    *cood=coo->getConstPointer();//Used for nodes coordinates

	DataArrayIdType *revNode =DataArrayIdType::New();
	DataArrayIdType *revNodeI=DataArrayIdType::New();
	mu->getReverseNodalConnectivity(revNode,revNodeI) ;
	const mcIdType *tmpN =revNode->getConstPointer();//Used to know which cells surround a given node
	const mcIdType *tmpNI=revNodeI->getConstPointer();

	DataArrayIdType *revCell =DataArrayIdType::New();
	DataArrayIdType *revCellI=DataArrayIdType::New();
	mu2->getReverseNodalConnectivity(revCell,revCellI);
	const mcIdType *tmpC =revCell->getConstPointer();//Used to know which faces surround a given node
	const mcIdType *tmpCI=revCellI->getConstPointer();

	const DataArrayIdType *nodal  = mu2->getNodalConnectivity() ;
	const DataArrayIdType *nodalI = mu2->getNodalConnectivityIndex() ;
	const mcIdType *tmpNE =nodal->getConstPointer();//Used to know which nodes surround a given face
	const mcIdType *tmpNEI=nodalI->getConstPointer();

	_numberOfCells = mu->getNumberOfCells() ;
	_cells      = std::shared_ptr<Cell>(new Cell[_numberOfCells], std::default_delete<Cell[]>()) ;

	_numberOfNodes = mu->getNumberOfNodes() ;//This number may include isolated nodes that will not be loaded. The number will be updated during nodes constructions
	_nodes      = std::shared_ptr<Node>(new Node[_numberOfNodes], std::default_delete<Node[]>())  ;//This array may be resized if isolated nodes are found

	_numberOfFaces = mu2->getNumberOfCells();
	_faces       = std::shared_ptr<Face>(new Face[_numberOfFaces], std::default_delete<Face[]>()) ;

    _indexFacePeriodicSet=false;

    //Definition used if _meshDim =3 to determine the edges
    DataArrayIdType *desc2 =DataArrayIdType::New();
    DataArrayIdType *descI2=DataArrayIdType::New();
    DataArrayIdType *revDesc2 =DataArrayIdType::New();
    DataArrayIdType *revDescI2=DataArrayIdType::New();
    DataArrayIdType *revNode2 =DataArrayIdType::New();
    DataArrayIdType *revNodeI2=DataArrayIdType::New();
    const mcIdType *tmpN2 ;
    const mcIdType *tmpNI2;
    MEDCouplingUMesh* mu3;
    
	if (_meshDim == 1)
        _numberOfEdges = mu->getNumberOfCells();
    else if (_meshDim == 2)
        _numberOfEdges = mu2->getNumberOfCells();
    else
    {
        mu3=mu2->buildDescendingConnectivity(desc2,descI2,revDesc2,revDescI2);//1D mesh of segments
        _numberOfEdges = mu3->getNumberOfCells();
        mu3->getReverseNodalConnectivity(revNode2,revNodeI2) ;
        tmpN2 =revNode2->getConstPointer();
        tmpNI2=revNodeI2->getConstPointer();
    }    

	// _cells, _nodes and _faces initialization:
	if (_meshDim == 1)
	{
		double xn, yn=0., zn=0.;//Components of the normal vector at a cell interface
		double norm;
		for( int id=0;id<_numberOfCells;id++ )
		{
			Point p(0.0,0.0,0.0) ;
			for(int idim=0; idim<_spaceDim; idim++)
				p[idim]=coorBary[id*_spaceDim+idim];
	
			mcIdType nbVertices=mu->getNumberOfNodesInCell(id) ;//should be equal to 2
			assert( nbVertices==2);
			std::vector<mcIdType> nodeIdsOfCell ;
			mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
	
	        mcIdType nbFaces=tmpI[id+1]-tmpI[id];//should be equal to 2
			assert( nbFaces==2);
	        const mcIdType *work=tmp+tmpI[id];
	
			/* compute the normal to the face */
	            xn = cood[nodeIdsOfCell[0]*_spaceDim  ] - cood[nodeIdsOfCell[nbFaces-1]*_spaceDim  ];
	        if(_spaceDim>1)        
				yn = cood[nodeIdsOfCell[0]*_spaceDim+1] - cood[nodeIdsOfCell[nbFaces-1]*_spaceDim+1];
	        if(_spaceDim>2)        
				zn = cood[nodeIdsOfCell[0]*_spaceDim+2] - cood[nodeIdsOfCell[nbFaces-1]*_spaceDim+2];
			norm = sqrt(xn*xn+yn*yn+zn*zn);
			if(norm<_epsilon)
				throw CdmathException("!!! Mesh::setMesh Normal vector has norm 0 !!!");
			else
			{
				xn /= norm;
				yn /= norm;
				zn /= norm;
			}
	        
			Cell ci( nbVertices, nbFaces, surf[id], p ) ;//nbCells=nbFaces=2
	        for( int el=0;el<nbFaces;el++ )
			{
				ci.addNodeId(el,nodeIdsOfCell[el]) ;//global node number
				ci.addNormalVector(el,xn,yn,zn) ;
				ci.addFaceId(el,work[el]) ;
				xn = - xn; yn=-yn; zn=-zn;
			}
			_cells.get()[id] = ci ;
		}
	
		for( int id(0); id<_numberOfFaces; id++ )
		{
			const mcIdType *workv=tmpNE+tmpNEI[id]+1;
			mcIdType nbNodes= tmpNEI[id+1]-tmpNEI[id]-1;//Normally equal to 1.
			assert( nbNodes==1);
	
			std::vector<double> coo(0) ;
			mu2->getCoordinatesOfNode(workv[0],coo);
			Point p(0,0.0,0.0) ;
			for(int idim=0; idim<_spaceDim; idim++)
				p[idim]=coo[idim];
	
		    const mcIdType *workc=tmpA+tmpAI[id];
		    mcIdType nbCells=tmpAI[id+1]-tmpAI[id];
			assert( nbCells>0);//To make sure our face is not located on an isolated node
		    
			Face fi( nbNodes, nbCells, 1.0, p, 1., 0., 0. ) ;
			for(int node_id=0; node_id<nbNodes;node_id++)//This loop could be deleted since nbNodes=1. Trying to merge with setMesh
				fi.addNodeId(node_id,workv[node_id]) ;//global node number
	
			fi.addCellId(0,workc[0]) ;
			for(int cell_id=1; cell_id<nbCells;cell_id++)
			{
				int cell_idx=0;
				if (workc[cell_id]!=workc[cell_id-1])//For some meshes (bad ones) the same cell can appear several times
					{
					fi.addCellId(cell_idx+1,workc[cell_id]) ;
					cell_idx++;
					}                
			}
			if(nbCells==1)
				_boundaryFaceIds.push_back(id);
			_faces.get()[id] = fi ;
		}
	
		int correctNbNodes=0;
		for( int id=0;id<_numberOfNodes;id++ )
		{
			const mcIdType *workc=tmpN+tmpNI[id];
			mcIdType nbCells=tmpNI[id+1]-tmpNI[id];
			
			if( nbCells>0)//To make sure this is not an isolated node
			{
				correctNbNodes++;
				std::vector<double> coo(0) ;
				mu->getCoordinatesOfNode(id,coo);
				Point p(0,0.0,0.0) ;
				for(int idim=0; idim<_spaceDim; idim++)
					p[idim]=coo[idim];
		
				const mcIdType *workf=tmpC+tmpCI[id];
				mcIdType nbFaces=tmpCI[id+1]-tmpCI[id];
				assert( nbFaces==1);
		
			    const mcIdType *workn=tmpN+tmpNI[id];
			    mcIdType nbNeighbourNodes=tmpNI[id+1]-tmpNI[id];
		        
				Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;
		        for( int el=0;el<nbCells;el++ )
					vi.addCellId(el,workc[el]) ;
		        for( int el=0;el<nbNeighbourNodes;el++ )
					vi.addNeighbourNodeId(el,workn[el]) ;//global node number
				for( int el=0;el<nbFaces;el++ )
					vi.addFaceId(el,workf[el],_faces.get()[workf[el]].isBorder()) ;
		 		if(vi.isBorder())
					_boundaryNodeIds.push_back(id);
				_nodes.get()[id] = vi ;
			}
		}
		if( _numberOfNodes!=correctNbNodes)//To do : reduce the size of pointer _nodes
		{
			cout<<"Found isolated nodes : correctNbNodes= "<<correctNbNodes<<", _numberOfNodes= "<<_numberOfNodes<<endl;
			_numberOfNodes = correctNbNodes;
			//memcpy(_nodes,mesh.getNodes(),correctNbNodes*sizeof(*mesh.getNodes())) ;
		}
	}
	else if(_meshDim==2  || _meshDim==3)
	{
		DataArrayDouble *barySeg = mu2->computeIsoBarycenterOfNodesPerCell();//computeCellCenterOfMass() ;//Used as face center
		const double *coorBarySeg=barySeg->getConstPointer();

		MEDCouplingFieldDouble* fieldl=mu2->getMeasureField(true);
		DataArrayDouble *longueur = fieldl->getArray();
		const double *lon=longueur->getConstPointer();//The lenght/area of each face

		MEDCouplingFieldDouble* fieldn;//The normal to each face
		DataArrayDouble *normal;
		const double *tmpNormal;

		if(_spaceDim==_meshDim)
			fieldn = mu2->buildOrthogonalField();//Compute the normal to each cell interface
		else
			fieldn = mu->buildOrthogonalField();//compute the 3D normal vector to the 2D cell
		
		normal = fieldn->getArray();
		tmpNormal = normal->getConstPointer();

		/*Building mesh cells */
		for(int id(0), k(0); id<_numberOfCells; id++, k+=_spaceDim)
		{
            const mcIdType *work=tmp+tmpI[id];      
			mcIdType nbFaces=tmpI[id+1]-tmpI[id];
            
			mcIdType nbVertices=mu->getNumberOfNodesInCell(id) ;

			vector<double> coorBaryXyz(3,0);
			for (int d=0; d<_spaceDim; d++)
				coorBaryXyz[d] = coorBary[k+d];

			Point p(coorBaryXyz[0],coorBaryXyz[1],coorBaryXyz[2]) ;
			Cell ci( nbVertices, nbFaces, surf[id], p ) ;

			/* Filling cell nodes */
			std::vector<mcIdType> nodeIdsOfCell ;
			mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
			for( int el=0;el<nbVertices;el++ )
				ci.addNodeId(el,nodeIdsOfCell[el]) ;

			/* Filling cell faces */
			if(_spaceDim==_meshDim)//use the normal field generated by buildOrthogonalField()
				for( int el=0;el<nbFaces;el++ )
				{
                    mcIdType faceIndex=(abs(work[el])-1);//=work[el] since Fortran type numbering was used, and negative sign means anticlockwise numbering
					vector<double> xyzn(3,0);//Outer normal to the cell
					if (work[el]<0)
						for (int d=0; d<_spaceDim; d++)
							xyzn[d] = -tmpNormal[_spaceDim*faceIndex+d];
					else
						for (int d=0; d<_spaceDim; d++)
							xyzn[d] = +tmpNormal[_spaceDim*faceIndex+d];
					ci.addNormalVector(el,xyzn[0],xyzn[1],xyzn[2]) ;
					ci.addFaceId(el,faceIndex) ;
				}
			else//build normals associated to the couple (cell id, face el)
			{//Case _meshDim=1 should be moved upper since we are in the 2D/3D branch
				if(_meshDim==1)//we know in this case there are only two faces around the cell id, each face is composed of a single node
				{//work[0]= first face global number, work[1]= second face global number
                    mcIdType indexFace0=abs(work[0])-1;//=work[0] since Fortran type numbering was used, and negative sign means anticlockwise numbering
                    mcIdType indexFace1=abs(work[1])-1;//=work[1] since Fortran type numbering was used, and negative sign means anticlockwise numbering
					mcIdType idNodeA=(tmpNE+tmpNEI[indexFace0]+1)[0];//global number of the first  face node work[0]=(abs(work[0])-1)
					mcIdType idNodeB=(tmpNE+tmpNEI[indexFace1]+1)[0];//global number of the second face node work[1]=(abs(work[1])-1)
					Vector vecAB(3);
					for(int i=0;i<_spaceDim;i++)
						vecAB[i]=coo->getIJ(idNodeB,i) - coo->getIJ(idNodeA,i);
					vecAB/=vecAB.norm();
					ci.addNormalVector(0,-vecAB[0],-vecAB[1],-vecAB[2]) ;
					ci.addNormalVector(1,vecAB[0],vecAB[1],vecAB[2]) ;				
					ci.addFaceId(0,indexFace0) ;
					ci.addFaceId(1,indexFace1) ;	
				}
				else//_meshDim==2, number of faces around the cell id is variable, each face is composed of two nodes
				{
					Vector xyzn(3);
					for (int d=0; d<_spaceDim; d++)
						xyzn[d] = tmpNormal[_spaceDim*id+d];
					for( int el=0;el<nbFaces;el++ )
					{
                        int faceIndex=(abs(work[el])-1);//=work[el] since Fortran type numbering was used, and negative sign means anticlockwise numbering
						const mcIdType *workv=tmpNE+tmpNEI[faceIndex]+1;
						mcIdType nbNodes= tmpNEI[faceIndex+1]-tmpNEI[faceIndex]-1;
						if(nbNodes!=2)//We want to compute the normal to a straight line, not a curved interface composed of more thant 2 points
						{
							cout<<"Mesh name "<< mu->getName()<< " space dim= "<< _spaceDim <<" mesh dim= "<< _meshDim <<endl;
							cout<<"For cell id "<<id<<" and local face number "<<el<<", the number of nodes is "<< nbNodes<< ", total number of faces is "<< nbFaces <<endl;
							throw CdmathException("Mesh::setMesh number of nodes around a face should be 2");
						}

						mcIdType idNodeA=workv[0];
						mcIdType idNodeB=workv[1];
						vector<double> nodeA(_spaceDim), nodeB(_spaceDim), nodeP(_spaceDim);
						for(int i=0;i<_spaceDim;i++)
						{
							nodeA[i]=coo->getIJ(idNodeA,i);
							nodeB[i]=coo->getIJ(idNodeB,i);
							nodeP[i]=coorBary[_spaceDim*id+i];
						}
						//Let P be the barycenter of the cell id
						Vector vecAB(3), vecPA(3);
						for(int i=0;i<_spaceDim;i++)
						{
							vecAB[i]=coo->getIJ(idNodeB,i)       - coo->getIJ(idNodeA,i);
							vecPA[i]=coo->getIJ(idNodeA,i) - coorBary[_spaceDim*id+i];
						}

						Vector normale = xyzn % vecAB;//Normal to the edge
						normale/=normale.norm();
                        
						if(normale*vecPA<0)
							ci.addNormalVector(el,normale[0],normale[1],normale[2]) ;	
						else
							ci.addNormalVector(el,-normale[0],-normale[1],-normale[2]) ;	
						ci.addFaceId(el,faceIndex) ;
					}
				}
			}
			_cells.get()[id] = ci ;
		}

		if(_spaceDim!=_meshDim)
		{
			/* Since spaceDim!=meshDim, don't build normal to faces */
			fieldn->decrRef();
            normal=NULL;
            tmpNormal=NULL;
		}

		/*Building mesh faces */
		for(int id(0), k(0); id<_numberOfFaces; id++, k+=_spaceDim)
		{
			vector<double> coorBarySegXyz(3);
			for (int d=0; d<_spaceDim; d++)
				coorBarySegXyz[d] = coorBarySeg[k+d];
			Point p(coorBarySegXyz[0],coorBarySegXyz[1],coorBarySegXyz[2]) ;
			const mcIdType *workc=tmpA+tmpAI[id];
			mcIdType nbCells=tmpAI[id+1]-tmpAI[id];
            
            if (nbCells>2 && _spaceDim==_meshDim)
            {
                cout<<"Warning : nbCells>2, numberOfFaces="<<_numberOfFaces<<endl;
                cout<<"nbCells= "<<nbCells<<", _spaceDim="<<_spaceDim<<", _meshDim="<<_meshDim<<endl;
                for(int icell=0; icell<nbCells; icell++)
                    cout<<workc[icell]<<", ";
                cout<<endl;
                throw CdmathException("Wrong mesh : nbCells>2 and spaceDim==meshDim");
            }
            if (nbCells==1)
                _boundaryFaceIds.push_back(id);
                
			const mcIdType *workv=tmpNE+tmpNEI[id]+1;
			mcIdType nbNodes= tmpNEI[id+1]-tmpNEI[id]-1;

			Face fi;
			if(_spaceDim==_meshDim)//Euclidean flat mesh geometry
                if(_spaceDim==2)
                    fi=Face( nbNodes, nbCells, lon[id], p, tmpNormal[k], tmpNormal[k+1], 0.0) ;
                else
                    fi=Face( nbNodes, nbCells, lon[id], p, tmpNormal[k], tmpNormal[k+1], tmpNormal[k+2]) ;
			else//Curved mesh geometry
				fi=Face( nbNodes, nbCells, lon[id], p, 0.0, 0.0, 0.0) ;//Since spaceDim!=meshDim, normal to face is not defined

			for(int node_id=0; node_id<nbNodes;node_id++)
				fi.addNodeId(node_id,workv[node_id]) ;

			fi.addCellId(0,workc[0]) ;
			for(int cell_id=1; cell_id<nbCells;cell_id++)
            {
                int cell_idx=0;
                if (workc[cell_id]!=workc[cell_id-1])//For some meshes (bad ones) the same cell can appear several times
                    {
                    fi.addCellId(cell_idx+1,workc[cell_id]) ;
                    cell_idx++;
                    }                
            }
            
			_faces.get()[id] = fi ;
		}

		/*Building mesh nodes, should be done after building mesh faces in order to detect boundary nodes*/
		int correctNbNodes=0;
		for(int id(0), k(0); id<_numberOfNodes; id++, k+=_spaceDim)
		{
			const mcIdType *workc=tmpN+tmpNI[id];
			mcIdType nbCells=tmpNI[id+1]-tmpNI[id];

			if( nbCells>0)//To make sure this is not an isolated node
			{
				correctNbNodes++;
				vector<double> coorP(3);
				for (int d=0; d<_spaceDim; d++)
					coorP[d] = cood[k+d];
				Point p(coorP[0],coorP[1],coorP[2]) ;
		
				const mcIdType *workf=tmpC+tmpCI[id];
				mcIdType nbFaces=tmpCI[id+1]-tmpCI[id];
				const mcIdType *workn;
				mcIdType nbNeighbourNodes;
				if (_meshDim == 1)
				{
					workn=tmpA+tmpAI[id];
					nbNeighbourNodes=tmpAI[id+1]-tmpAI[id];
				}
				else if (_meshDim == 2)
				{
					workn=tmpC+tmpCI[id];
					nbNeighbourNodes=tmpCI[id+1]-tmpCI[id];
				}
				else//_meshDim == 3
				{
					workn=tmpN2+tmpNI2[id];
					nbNeighbourNodes=tmpNI2[id+1]-tmpNI2[id];
				}    
				Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;
		
				for( int el=0;el<nbCells;el++ )
					vi.addCellId(el,workc[el]) ;
				for( int el=0;el<nbNeighbourNodes;el++ )
					vi.addNeighbourNodeId(el,workn[el]) ;
				//Detection of border nodes    
				for( int el=0;el<nbFaces;el++ )
					vi.addFaceId(el,workf[el],_faces.get()[workf[el]].isBorder()) ;
				if(vi.isBorder())
					_boundaryNodeIds.push_back(id);
				_nodes.get()[id] = vi ;
			}
		}
		if( _numberOfNodes!=correctNbNodes)//To do : reduce the size of pointer _nodes
		{
			cout<<"Found isolated nodes : correctNbNodes= "<<correctNbNodes<<", _numberOfNodes= "<<_numberOfNodes<<endl;
			_numberOfNodes = correctNbNodes;
		}
		
		if(_spaceDim==_meshDim)
			fieldn->decrRef();
		fieldl->decrRef();
		barySeg->decrRef();
	}
	else
		throw CdmathException("Mesh::setMesh space dimension should be 1, 2 or 3");

    //Set boundary groups
    _faceGroupNames.push_back("Boundary");
    _nodeGroupNames.push_back("Boundary");
    _faceGroupsIds.push_back(_boundaryFaceIds);
    _nodeGroupsIds.push_back(_boundaryNodeIds);
    if( _meshDim>1 )//Set face boundary group
		_boundaryMesh = mu->computeSkin();
    else//in dimension 1 computeSkin leads to a degenerate mesh
		_boundaryMesh = NULL;

    _faceGroups.push_back(_boundaryMesh);
    _nodeGroups.push_back(NULL);

    desc->decrRef();
	descI->decrRef();
	revDesc->decrRef();
	revDescI->decrRef();
	mu2->decrRef();
	baryCell->decrRef();
	fields->decrRef();
	revNode->decrRef();
	revNodeI->decrRef();
	revCell->decrRef();
	revCellI->decrRef();

    // Arrays used if _meshDim =3 to determine the edges
	revNode2->decrRef();
	revNodeI2->decrRef();
	desc2->decrRef();
	descI2->decrRef();
	revDesc2->decrRef();
	revDescI2->decrRef();

    if (_meshDim == 3)
        mu3->decrRef();
    	
    return mu;
}

void
Mesh::setGroupAtFaceByCoords(double x, double y, double z, double eps, std::string groupName, bool isBoundaryGroup)
{
	std::vector< int > faceIds(0);
	double FX, FY, FZ;
	Face Fi;
	
	/* Construction of the face group */
	if(isBoundaryGroup)        
        for(int i=0; i<_boundaryFaceIds.size(); i++)
        {
			Fi=_faces.get()[_boundaryFaceIds[i]];
			FX=Fi.x();
			FY=Fi.y();
			FZ=Fi.z();
			if (abs(FX-x)<eps && abs(FY-y)<eps && abs(FZ-z)<eps)
			{
				faceIds.insert(faceIds.end(),_boundaryFaceIds[i]);
				_faces.get()[_boundaryFaceIds[i]].setGroupName(groupName);
			}
        }
	else
		for (int iface=0;iface<_numberOfFaces;iface++)
		{
			Fi=_faces.get()[iface];
			FX=Fi.x();
			FY=Fi.y();
			FZ=Fi.z();
			if (abs(FX-x)<eps && abs(FY-y)<eps && abs(FZ-z)<eps)
			{
				faceIds.insert(faceIds.end(),iface);
				_faces.get()[iface].setGroupName(groupName);
			}
		}

	if (faceIds.size()>0)
    {
		std::vector< std::string >::iterator it = std::find(_faceGroupNames.begin(), _faceGroupNames.end(), groupName);
		if(it == _faceGroupNames.end())//No group named groupName
		{
			_faceGroupNames.insert(_faceGroupNames.end(),groupName);
			_faceGroupsIds.insert(  _faceGroupsIds.end(),faceIds);
			_faceGroups.insert(    _faceGroups.end(), NULL);//No mesh created. Create one ?
		}
		else
		{
			std::vector< int > faceGroupIds = _faceGroupsIds[it-_faceGroupNames.begin()];
			faceGroupIds.insert( faceGroupIds.end(), faceIds.begin(), faceIds.end());
			/* Detect and erase duplicates face ids */
			sort( faceGroupIds.begin(), faceGroupIds.end() );
			faceGroupIds.erase( unique( faceGroupIds.begin(), faceGroupIds.end() ), faceGroupIds.end() );
			_faceGroupsIds[it-_faceGroupNames.begin()] = faceGroupIds;
		}
	}
}

void
Mesh::setGroupAtNodeByCoords(double x, double y, double z, double eps, std::string groupName, bool isBoundaryGroup)
{
	std::vector< int > nodeIds(0);
	double NX, NY, NZ;
	Node Ni;
	
	/* Construction of the node group */
	if(isBoundaryGroup)        
        for(int i=0; i<_boundaryNodeIds.size(); i++)
        {
			Ni=_nodes.get()[_boundaryNodeIds[i]];
			NX=Ni.x();
			NY=Ni.y();
			NZ=Ni.z();
			if (abs(NX-x)<eps && abs(NY-y)<eps && abs(NZ-z)<eps)
			{
				nodeIds.insert(nodeIds.end(),_boundaryNodeIds[i]);
				_nodes.get()[_boundaryNodeIds[i]].setGroupName(groupName);
			}
        }
	else
		for (int inode=0;inode<_numberOfNodes;inode++)
		{
			NX=_nodes.get()[inode].x();
			NY=_nodes.get()[inode].y();
			NZ=_nodes.get()[inode].z();
			if (abs(NX-x)<eps && abs(NY-y)<eps && abs(NZ-z)<eps)
			{
				nodeIds.insert(nodeIds.end(),inode);
				_nodes.get()[inode].setGroupName(groupName);
			}
		}

	if (nodeIds.size()>0)
    {
		std::vector< std::string >::iterator it = std::find(_nodeGroupNames.begin(), _nodeGroupNames.end(), groupName);
		if(it == _nodeGroupNames.end())//No group named groupName
		{
			_nodeGroupNames.insert(_nodeGroupNames.end(),groupName);
			_nodeGroupsIds.insert(  _nodeGroupsIds.end(),nodeIds);
			_nodeGroups.insert(    _nodeGroups.end(), NULL);//No mesh created. Create one ?
		}
		else
		{
			std::vector< int > nodeGroupIds = _nodeGroupsIds[it-_nodeGroupNames.begin()];
			nodeGroupIds.insert( nodeGroupIds.end(), nodeIds.begin(), nodeIds.end());
			/* Detect and erase duplicates node ids */
			sort( nodeGroupIds.begin(), nodeGroupIds.end() );
			nodeGroupIds.erase( unique( nodeGroupIds.begin(), nodeGroupIds.end() ), nodeGroupIds.end() );
			_nodeGroupsIds[it-_nodeGroupNames.begin()] = nodeGroupIds;
		}
    }
}

void
Mesh::setGroupAtPlan(double value, int direction, double eps, std::string groupName, bool isBoundaryGroup)
{
	std::vector< int > faceIds(0), nodeIds(0);
	double cord;
	
	/* Construction of the face group */	
	if(isBoundaryGroup)        
        for(int i=0; i<_boundaryFaceIds.size(); i++)
        {
			cord=_faces.get()[_boundaryFaceIds[i]].getBarryCenter()[direction];
			if (abs(cord-value)<eps)
			{
				faceIds.insert(faceIds.end(),_boundaryFaceIds[i]);
				_faces.get()[_boundaryFaceIds[i]].setGroupName(groupName);
			}
        }
	else
		for (int iface=0;iface<_numberOfFaces;iface++)
		{
			cord=_faces.get()[iface].getBarryCenter()[direction];
			if (abs(cord-value)<eps)
			{
				faceIds.insert(faceIds.end(),iface);
				_faces.get()[iface].setGroupName(groupName);
			}
		}

	/* Construction of the node group */
	if(isBoundaryGroup)        
        for(int i=0; i<_boundaryNodeIds.size(); i++)
        {
			cord=_nodes.get()[_boundaryNodeIds[i]].getPoint()[direction];
			if (abs(cord-value)<eps)
			{
				nodeIds.insert(nodeIds.end(),_boundaryNodeIds[i]);
				_nodes.get()[_boundaryNodeIds[i]].setGroupName(groupName);
			}
        }
	else
		for (int inode=0;inode<_numberOfNodes;inode++)
		{
			cord=_nodes.get()[inode].getPoint()[direction];
			if (abs(cord-value)<eps)
			{
				nodeIds.insert(nodeIds.end(),inode);
				_nodes.get()[inode].setGroupName(groupName);
			}
		}

	if (faceIds.size()>0)
    {
		std::vector< std::string >::iterator it = std::find(_faceGroupNames.begin(), _faceGroupNames.end(), groupName);
		if(it == _faceGroupNames.end())//No group named groupName
		{
			_faceGroupNames.insert(_faceGroupNames.end(),groupName);
			_faceGroupsIds.insert(  _faceGroupsIds.end(),faceIds);
			_faceGroups.insert(    _faceGroups.end(), NULL);//No mesh created. Create one ?
		}
		else
		{
			std::vector< int > faceGroupIds = _faceGroupsIds[it-_faceGroupNames.begin()];
			faceGroupIds.insert( faceGroupIds.end(), faceIds.begin(), faceIds.end());
			/* Detect and erase duplicates face ids */
			sort( faceGroupIds.begin(), faceGroupIds.end() );
			faceGroupIds.erase( unique( faceGroupIds.begin(), faceGroupIds.end() ), faceGroupIds.end() );
			_faceGroupsIds[it-_faceGroupNames.begin()] = faceGroupIds;
		}
	}
	if (nodeIds.size()>0)
    {
		std::vector< std::string >::iterator it = std::find(_nodeGroupNames.begin(), _nodeGroupNames.end(), groupName);
		if(it == _nodeGroupNames.end())//No group named groupName
		{
			_nodeGroupNames.insert(_nodeGroupNames.end(),groupName);
			_nodeGroupsIds.insert(  _nodeGroupsIds.end(),nodeIds);
			_nodeGroups.insert(    _nodeGroups.end(), NULL);//No mesh created. Create one ?
		}
		else
		{
			std::vector< int > nodeGroupIds = _nodeGroupsIds[it-_nodeGroupNames.begin()];
			nodeGroupIds.insert( nodeGroupIds.end(), nodeIds.begin(), nodeIds.end());
			/* Detect and erase duplicates node ids */
			sort( nodeGroupIds.begin(), nodeGroupIds.end() );
			nodeGroupIds.erase( unique( nodeGroupIds.begin(), nodeGroupIds.end() ), nodeGroupIds.end() );
			_nodeGroupsIds[it-_nodeGroupNames.begin()] = nodeGroupIds;
		}
    }
}

void
Mesh::setBoundaryNodesFromFaces()
{
    for (int iface=0;iface<_boundaryFaceIds.size();iface++)
    {
        std::vector< int > nodesID= _faces.get()[_boundaryFaceIds[iface]].getNodesId();
        int nbNodes = _faces.get()[_boundaryFaceIds[iface]].getNumberOfNodes();
        for(int inode=0 ; inode<nbNodes ; inode++)
        {
            std::vector<int>::const_iterator  it = std::find(_boundaryNodeIds.begin(),_boundaryNodeIds.end(),nodesID[inode]);
            if( it != _boundaryNodeIds.end() )
                _boundaryNodeIds.push_back(nodesID[inode]);
        }
    }
}

std::map<int,int>
Mesh::getIndexFacePeriodic( void ) const
{
    return _indexFacePeriodicMap;
}

void
Mesh::setPeriodicFaces(bool check_groups, bool use_central_inversion)
{
    if(_indexFacePeriodicSet)
        return;
        
    for (int indexFace=0;indexFace<_boundaryFaceIds.size() ; indexFace++)
    {
        Face my_face=_faces.get()[_boundaryFaceIds[indexFace]];
        int iface_perio=-1;
        if(_meshDim==1)
        {
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
                if(iface!=indexFace)
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
        }
        else if(_meshDim==2)
        {
            double x=my_face.x();
            double y=my_face.y();
            
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
            {
                Face face_i=_faces.get()[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                if (   (abs(y-yi)<_epsilon || abs(x-xi)<_epsilon )// Case of a square geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi)<_epsilon ) // Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<_epsilon
                    && fabs(my_face.getXN()      + face_i.getXN())<_epsilon
                    && fabs(my_face.getYN()      + face_i.getYN())<_epsilon
                    && fabs(my_face.getZN()      + face_i.getZN())<_epsilon )
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
            }
        }
        else if(_meshDim==3)
        {
            double x=my_face.x();
            double y=my_face.y();
            double z=my_face.z();
        
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
            {
                Face face_i=_faces.get()[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                double zi=face_i.z();
                if ( ((abs(y-yi)<_epsilon && abs(x-xi)<_epsilon) || (abs(x-xi)<_epsilon && abs(z-zi)<_epsilon) || (abs(y-yi)<_epsilon && abs(z-zi)<_epsilon))// Case of a cube geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi) + abs(z+zi)<_epsilon )// Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<_epsilon
                    && fabs(my_face.getXN()      + face_i.getXN())<_epsilon
                    && fabs(my_face.getYN()      + face_i.getYN())<_epsilon
                    && fabs(my_face.getZN()      + face_i.getZN())<_epsilon )
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
            }  
        }
        else
            throw CdmathException("Mesh::setPeriodicFaces: Mesh dimension should be 1, 2 or 3");
        
        if (iface_perio==-1)
            throw CdmathException("Mesh::setPeriodicFaces: periodic face not found, iface_perio==-1 " );
        else
            _indexFacePeriodicMap[_boundaryFaceIds[indexFace]]=iface_perio;
    }
    _indexFacePeriodicSet=true;    
}

int
Mesh::getIndexFacePeriodic(int indexFace, bool check_groups, bool use_central_inversion)
{
	if (!_faces.get()[indexFace].isBorder())
        {
            cout<<"Pb with indexFace= "<<indexFace<<endl;
            throw CdmathException("Mesh::getIndexFacePeriodic: not a border face" );
        }
        
    if(!_indexFacePeriodicSet)
        setPeriodicFaces(check_groups, use_central_inversion);

    std::map<int,int>::const_iterator  it = _indexFacePeriodicMap.find(indexFace);
    if( it != _indexFacePeriodicMap.end() )
        return it->second;
    else
    {
        cout<<"Pb with indexFace= "<<indexFace<<endl;
        throw CdmathException("Mesh::getIndexFacePeriodic: not a periodic face" );
    }
}

bool
Mesh::isBorderNode(int nodeid) const
{
	return _nodes.get()[nodeid].isBorder();
}

bool
Mesh::isBorderFace(int faceid) const
{
	return _faces.get()[faceid].isBorder();
}

std::vector< int > 
Mesh::getBoundaryFaceIds() const
{
    return _boundaryFaceIds;
}

std::vector< int > 
Mesh::getBoundaryNodeIds() const
{
    return _boundaryNodeIds;
}

bool
Mesh::isTriangular() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_TRI3;
}
bool
Mesh::isTetrahedral() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_TETRA4;
}
bool
Mesh::isQuadrangular() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_QUAD4;
}
bool
Mesh::isHexahedral() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_HEXA8;
}
bool
Mesh::isStructured() const
{
	return _isStructured;
}

std::vector< INTERP_KERNEL::NormalizedCellType > 
Mesh::getElementTypes() const
{
  return _eltsTypes;  
}

std::vector< string >
Mesh::getElementTypesNames() const
{
    std::vector< string > result(0);    
    for(int i=0; i< _eltsTypes.size(); i++)
    {
        if( _eltsTypes[i]==INTERP_KERNEL::NORM_POINT1)
            result.push_back("Points");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_SEG2)
            result.push_back("Segments");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_TRI3)
            result.push_back("Triangles");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_QUAD4)
            result.push_back("Quadrangles");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_POLYGON)
            result.push_back("Polygons");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_TETRA4)
            result.push_back("Tetrahedra");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_PYRA5)
            result.push_back("Pyramids");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_PENTA6)
            result.push_back("Pentahedra");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_HEXA8)
            result.push_back("Hexahedra");
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_POLYHED)
            result.push_back("Polyhedrons");
        else
		{
			cout<< "Mesh " + _name + " contains an element of type " <<endl;
			cout<< _eltsTypes[i]<<endl;
			throw CdmathException("Mesh::getElementTypes : recognised cell med types are NORM_POINT1, NORM_SEG2, NORM_TRI3, NORM_QUAD4, NORM_TETRA4, NORM_PYRA5, NORM_PENTA6, NORM_HEXA8, NORM_POLYGON, NORM_POLYHED");
        }
    }
    return result;
}

void
Mesh::setFaceGroups( const MEDFileUMesh* medmesh, MEDCouplingUMesh*  mu)
{
	int nbCellsSubMesh;
	bool foundFace;
	
	/* Searching for face groups */
	vector<string> faceGroups=medmesh->getGroupsNames() ;

	for (unsigned int i=0;i<faceGroups.size();i++ )
	{
		string groupName=faceGroups[i];
		vector<mcIdType> nonEmptyGrp(medmesh->getGrpNonEmptyLevels(groupName));
		//We check if the group has a relative dimension equal to -1 
		//before call to the function getGroup(-1,groupName.c_str())
		vector<mcIdType>::iterator it = find(nonEmptyGrp.begin(), nonEmptyGrp.end(), -1);
		if (it != nonEmptyGrp.end())
		{
			MEDCouplingUMesh *m=medmesh->getGroup(-1,groupName.c_str());
            m->unPolyze();
			nbCellsSubMesh=m->getNumberOfCells();
			
			_faceGroups.insert(_faceGroups.end(),m);//Vector of group meshes
			_faceGroupNames.insert(_faceGroupNames.end(),groupName);//Vector of group names
			_faceGroupsIds.insert(_faceGroupsIds.end(),std::vector<int>(nbCellsSubMesh));//Vector of group face Ids. The filling of the face ids takes place below.
			
			DataArrayDouble *baryCell = m->computeIsoBarycenterOfNodesPerCell();//computeCellCenterOfMass() ;
			const double *coorBary=baryCell->getConstPointer();
			// Face identification 
			for (int ic(0), k(0); ic<nbCellsSubMesh; ic++, k+=_spaceDim)
			{
				vector<double> coorBaryXyz(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorBaryXyz[d] = coorBary[k+d];
				Point p1(coorBaryXyz[0],coorBaryXyz[1],coorBaryXyz[2]) ;

				foundFace=false;
				for (int iface=0;iface<_numberOfFaces;iface++ )
				{
					Point p2=_faces.get()[iface].getBarryCenter();
					if(p1.distance(p2)<_epsilon)
					{
						_faces.get()[iface].setGroupName(groupName);
						_faceGroupsIds[_faceGroupsIds.size()-1][ic]=iface;
						foundFace=true;
						break;
					}
				}
				if (not foundFace)
					throw CdmathException("No face found for group " + groupName );
			}
			baryCell->decrRef();
			//m->decrRef();
		}
	}
}
void
Mesh::setNodeGroups( const MEDFileMesh* medmesh, MEDCouplingUMesh*  mu)
{
	int  nbNodesSubMesh;
	bool foundNode;
	
	/* Searching for node groups */
	vector<string> nodeGroups=medmesh->getGroupsOnSpecifiedLev(1) ;

	for (unsigned int i=0;i<nodeGroups.size();i++ )
	{
		string groupName=nodeGroups[i];
		DataArrayIdType * nodeGroup=medmesh->getNodeGroupArr( groupName );
		const mcIdType *nodeids=nodeGroup->getConstPointer();

		if(nodeids!=NULL)
		{
			_nodeGroups.insert(_nodeGroups.end(),nodeGroup);
			_nodeGroupNames.insert(_nodeGroupNames.end(),groupName);

			nbNodesSubMesh=nodeGroup->getNumberOfTuples();//nodeGroup->getNbOfElems();

			DataArrayDouble *coo = mu->getCoords() ;
			const double *cood=coo->getConstPointer();

			_nodeGroupsIds.insert(_nodeGroupsIds.end(),std::vector<int>(nbNodesSubMesh));//Vector of boundary faces
			/* Node identification */
			for (int ic(0); ic<nbNodesSubMesh; ic++)
			{
				vector<double> coorP(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorP[d] = cood[nodeids[ic]*_spaceDim+d];
				Point p1(coorP[0],coorP[1],coorP[2]) ;

				foundNode=false;
				for (int inode=0;inode<_numberOfNodes;inode++ )
				{
					Point p2=_nodes.get()[inode].getPoint();
					if(p1.distance(p2)<_epsilon)
					{
						_nodes.get()[inode].setGroupName(groupName);
						_nodeGroupsIds[_nodeGroupsIds.size()-1][ic]=inode;
						foundNode=true;
						break;
					}
				}
				if (not foundNode)
					throw CdmathException("No node found for group " + groupName );
			}
		}
	nodeGroup->decrRef();
	}
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

std::vector<mcIdType>
Mesh::getCellGridStructure() const
{
    if(!_isStructured)
		throw CdmathException("std::vector<int> Mesh::getCellGridStructure() : nx, ny and nz are defined only for structured meshes !");

	return _nxyz;
}

//----------------------------------------------------------------------
int
Mesh::getNx( void )  const
//----------------------------------------------------------------------
{
    if(!_isStructured)
		throw CdmathException("int Mesh::getNx( void ) : Nx is defined only for structured meshes !");

	return _nxyz[0];
}

//----------------------------------------------------------------------
int
Mesh::getNy( void )  const
//----------------------------------------------------------------------
{
    if(!_isStructured)
		throw CdmathException("int Mesh::getNy( void ) : Ny is defined only for structured meshes !");
	if(_meshDim < 2)
		throw CdmathException("int Mesh::getNy( void ) : Ny is not defined for mesh dimension < 2!");

	return _nxyz[1];
}

//----------------------------------------------------------------------
int
Mesh::getNz( void )  const
//----------------------------------------------------------------------
{
    if(!_isStructured)
		throw CdmathException("int Mesh::getNz( void ) : Nz is defined only for structured meshes !");
	if(_meshDim < 3)
		throw CdmathException("int Mesh::getNz( void ) : Nz is not defined for mesh dimension < 3!");

	return _nxyz[2];
}

//----------------------------------------------------------------------
double
Mesh::getXMin( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[0] ;
}

//----------------------------------------------------------------------
double
Mesh::getXMax( void )  const
//----------------------------------------------------------------------
{
	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[1] ;
}

//----------------------------------------------------------------------
double
Mesh::getYMin( void )  const
//----------------------------------------------------------------------
{
	if(_meshDim<2)
		throw CdmathException("Mesh::getYMin : dimension should be >=2");
		
	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[2] ;
}

//----------------------------------------------------------------------
double
Mesh::getYMax( void )  const
//----------------------------------------------------------------------
{
	if(_meshDim<2)
		throw CdmathException("Mesh::getYMax : dimension should be >=2");
		
	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[3] ;
}

//----------------------------------------------------------------------
double
Mesh::getZMin( void )  const
//----------------------------------------------------------------------
{
	if(_meshDim<3)
		throw CdmathException("Mesh::getZMin : dimension should be 3");
		
	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[4] ;
}

//----------------------------------------------------------------------
double
Mesh::getZMax( void )  const
//----------------------------------------------------------------------
{
	if(_meshDim<3)
		throw CdmathException("Mesh::getZMax : dimension should be 3");

	double Box0[2*_meshDim];
    _mesh->getBoundingBox(Box0);

	return Box0[5] ;
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingMesh>
Mesh::getMEDCouplingMesh( void )  const
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
int
Mesh::getNumberOfEdges ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfEdges ;
}

//----------------------------------------------------------------------
std::shared_ptr<Face>
Mesh::getFaces ( void )  const
//----------------------------------------------------------------------
{
	return _faces ;
}

//----------------------------------------------------------------------
std::shared_ptr<Cell>
Mesh::getCells ( void ) const
//----------------------------------------------------------------------
{
	return _cells ;
}

//----------------------------------------------------------------------
Cell&
Mesh::getCell ( int i ) const
//----------------------------------------------------------------------
{
	return _cells.get()[i] ;
}

//----------------------------------------------------------------------
Face&
Mesh::getFace ( int i ) const
//----------------------------------------------------------------------
{
	return _faces.get()[i] ;
}

//----------------------------------------------------------------------
Node&
Mesh::getNode ( int i ) const
//----------------------------------------------------------------------
{
	return _nodes.get()[i] ;
}

//----------------------------------------------------------------------
std::shared_ptr<Node>
Mesh::getNodes ( void )  const
//----------------------------------------------------------------------
{
	return _nodes ;
}

vector<string>
Mesh::getNameOfFaceGroups( void )  const
{
	return _faceGroupNames;
}

vector< std::vector<int> >
Mesh::getFaceGroups( void )  const
{
	return _faceGroupsIds;
}

vector<MEDCoupling::MEDCouplingUMesh *>
Mesh::getMEDCouplingFaceGroups( void )  const
{
	return _faceGroups;
}

vector<string>
Mesh::getNameOfNodeGroups( void )  const
{
	return _nodeGroupNames;
}

vector< std::vector<int> >
Mesh::getNodeGroups( void )  const
{
	return _nodeGroupsIds;
}


vector<MEDCoupling::DataArrayIdType *>
Mesh::getMEDCouplingNodeGroups( void )  const
{
	return _nodeGroups;
}

//----------------------------------------------------------------------
const Mesh&
Mesh::operator= ( const Mesh& mesh )
//----------------------------------------------------------------------
{
	_spaceDim = mesh.getSpaceDimension() ;
	_meshDim  = mesh.getMeshDimension() ;
    _name = mesh.getName();
    _epsilon=mesh.getComparisonEpsilon();
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
	_numberOfEdges = mesh.getNumberOfEdges();
    
    _isStructured = mesh.isStructured();
    _meshNotDeleted = mesh.meshNotDeleted();
    
    if(_isStructured)
        _nxyz = mesh.getCellGridStructure() ;

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroupsIds = mesh.getFaceGroups() ;
	_faceGroups = mesh.getMEDCouplingFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroupsIds = mesh.getNodeGroups() ;
	_nodeGroups = mesh.getMEDCouplingNodeGroups() ;

	_nodes   = mesh.getNodes() ;
	_faces   = mesh.getFaces() ;
	_cells   = mesh.getCells() ;

    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();
    
    _boundaryFaceIds=mesh.getBoundaryFaceIds();
    _boundaryNodeIds=mesh.getBoundaryNodeIds();

    _eltsTypes=mesh.getElementTypes();
    _eltsTypesNames=mesh.getElementTypesNames();

	_mesh=mesh.getMEDCouplingMesh()->clone(false);//Clone because you will need to buildUnstructured. No deep copy : it is assumed node coordinates and cell connectivity will not change

	return *this;
}

bool Mesh::isIndexFacePeriodicSet() const
{
 return    _indexFacePeriodicSet;
}
//----------------------------------------------------------------------
double 
Mesh::minRatioVolSurf() const
{
    double dx_min  = 1e30;
    for(int i=0; i<_numberOfCells; i++)
    {
        Cell Ci = getCell(i);
        if (_meshDim > 1)
        {
            double perimeter=0;
            for(int k=0; k< Ci.getNumberOfFaces(); k++)
            {
                int indexFace=Ci.getFacesId()[k];
                Face Fk = getFace(indexFace);
                perimeter+=Fk.getMeasure();
            }
            dx_min = min(dx_min,Ci.getMeasure()/perimeter);
        }
        else
            dx_min = min(dx_min,Ci.getMeasure());
    }
    
    return dx_min;
}

Mesh
Mesh::getBoundaryMesh ( void )  const 
{
	if( _boundaryMesh )
		return Mesh(_boundaryMesh);
	else
		throw CdmathException("Mesh::getBoundaryMesh Boundary mesh is empty");
}

MEDCoupling::MEDCouplingUMesh *
Mesh::getBoundaryMEDCouplingMesh ( void )  const 
{
	return _boundaryMesh;
}

Mesh 
Mesh::getBoundaryGroupMesh ( std::string groupName, int nth_occurence )  const 
{
	//count occurences of groupName in known group name list
	int count_occurences = std::count (_faceGroupNames.begin(),_faceGroupNames.end(),groupName);
	
    if( count_occurences ==0 )//No group found
    {
        cout<<"Mesh::getBoundaryGroupMesh Error : face group " << groupName << " does not exist"<<endl;
        cout<<"Known face group names are " ;
        for(int i=0; i<_faceGroupNames.size(); i++)
			cout<< _faceGroupNames[i]<<" ";
		cout<< endl;
        throw CdmathException("Required face group does not exist");
    }
    else if ( count_occurences <= nth_occurence)//Too many groups found
    {
        cout<<"Mesh::getBoundaryGroupMesh Error : "<<count_occurences<<" groups have name " << groupName<<", but you asked fo occurencer number "<<nth_occurence<<"which is too large"<<endl;
        cout<<"Call function getBoundaryGroupMesh ( string groupName, int nth_group_match ) with nth_group_match between 0 and "<<count_occurences-1<<" to discriminate them "<<endl ;
        throw CdmathException("Several face groups have the same name but you asked for an occurence that does not exsit");
    }
    else if( count_occurences >1 )//Wrning several groups found
		cout<<"Warning : "<<count_occurences<<" groups have name " << groupName<<". Searching occurence number "<<nth_occurence<<endl;
		
	//search occurence of group name in known group name list
	std::vector<std::string>::const_iterator it = _faceGroupNames.begin();
	for (int i = 0; i<nth_occurence+1; i++)
		it = std::find(it,_faceGroupNames.end(),groupName);
		
	return Mesh(_faceGroups[it - _faceGroupNames.begin()]);	
}

int 
Mesh::getMaxNbNeighbours(EntityType type) const
{
    int result=0;
    
    if (type==CELLS)
	{
		int nbNeib;//local number of neighbours
        for(int i=0; i<_numberOfCells; i++)
        {
            Cell Ci = _cells.get()[i];
            //Careful with mesh with junctions : do not just take Ci.getNumberOfFaces()
            nbNeib=0;
            for(int j=0; j<Ci.getNumberOfFaces(); j++)
                nbNeib+=_faces.get()[Ci.getFacesId()[j]].getNumberOfCells()-1;//Without junction this would be +=1
            
            if(result < nbNeib)
                result=nbNeib;
		}
	}
    else if(type==NODES)
	{
        for(int i=0; i<_numberOfNodes; i++)
            if(result < _nodes.get()[i].getNumberOfEdges())
                result=_nodes.get()[i].getNumberOfEdges();
	}
    else
		throw CdmathException("Mesh::getMaxNbNeighbours : entity type is not accepted. Should be CELLS or NODES");

    return result;
}
//----------------------------------------------------------------------
void
Mesh::writeVTK ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	if( !_meshNotDeleted )
		throw CdmathException("Mesh::writeVTK : Cannot save mesh : no MEDCouplingMesh loaded (may be deleted)");
		
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
Mesh::writeMED ( const std::string fileName, bool fromScratch ) const
//----------------------------------------------------------------------
{
	if( !_meshNotDeleted )
		throw CdmathException("Mesh::writeMED : Cannot save mesh : no MEDCouplingMesh loaded (may be deleted)");

	string fname=fileName+".med";
	if(_isStructured)//Check if we have a medcouplingimesh that can't be written directly
	{
		const MEDCoupling::MEDCouplingIMesh* iMesh = dynamic_cast< const MEDCoupling::MEDCouplingIMesh* > ((const MEDCoupling::MEDCouplingMesh*) _mesh);
		if(iMesh)//medcouplingimesh : Use convertToCartesian in order to write mesh
		{
			MEDCouplingCMesh * cMesh = iMesh->convertToCartesian();
			MEDCoupling::WriteMesh(fname.c_str(), cMesh, fromScratch);
			cMesh->decrRef();
		}
		else//medcouplingcmesh : save directly
			MEDCoupling::WriteMesh(fname.c_str(),_mesh, fromScratch);
	}
	else
		MEDCoupling::WriteMesh(fname.c_str(),_mesh, fromScratch);

	/* Try to save mesh groups */
	//MEDFileUMesh meshMEDFile;
	//meshMEDFile.setMeshAtLevel(0,mu);
	//for(int i=0; i< _faceGroups.size(); i++)
	//meshMEDFile.setMeshAtLevel(-1,_faceGroups[i]);
	//if (fromScratch)
	//MEDCoupling::meshMEDFile.write(fname.c_str(),2)	;
	//else
	//MEDCoupling::meshMEDFile.write(fname.c_str(),1)	;
}

std::vector< int > 
Mesh::getFaceGroupIds(std::string groupName, bool isBoundaryGroup) const
{
	std::vector<std::string>::const_iterator  it = std::find(_faceGroupNames.begin(),_faceGroupNames.end(),groupName);
    if( it == _faceGroupNames.end() )
    {
        cout<<"Mesh::getGroupFaceIds Error : face group " << groupName << " does not exist"<<endl;
        throw CdmathException("Required face group does not exist");
    }
    else
    {   
		return _faceGroupsIds[it-_faceGroupNames.begin()];
	}
}

std::vector< int > 
Mesh::getNodeGroupIds(std::string groupName, bool isBoundaryGroup) const
{
	std::vector<std::string>::const_iterator  it = std::find(_nodeGroupNames.begin(),_nodeGroupNames.end(),groupName);
    if( it == _nodeGroupNames.end() )
    {
        cout<<"Mesh::getGroupNodeIds Error : node group " << groupName << " does not exist"<<endl;
        throw CdmathException("Required node group does not exist");
    }
    else
        return _nodeGroupsIds[it-_nodeGroupNames.begin()];
}

void 
Mesh::setFaceGroupByIds(std::vector< int > faceIds, std::string groupName) 
{
    for(int i=0; i< faceIds.size(); i++)
        getFace(faceIds[i]).setGroupName(groupName);
        
    _faceGroupsIds.insert(  _faceGroupsIds.end(),faceIds);
    _faceGroupNames.insert(_faceGroupNames.end(), groupName);
    _faceGroups.insert(    _faceGroups.end(), NULL);
}

void 
Mesh::setNodeGroupByIds(std::vector< int > nodeIds, std::string groupName) 
{
    for(int i=0; i< nodeIds.size(); i++)
        getNode(nodeIds[i]).setGroupName(groupName);
}

void Mesh::deleteMEDCouplingMesh()
{ 
	if(_meshNotDeleted) 
	{
		(_mesh.retn())->decrRef(); 
		_meshNotDeleted=false;
	} 
	else 
		throw CdmathException("Mesh::deleteMEDCouplingMesh() : mesh is not loaded");
};
