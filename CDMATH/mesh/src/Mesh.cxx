/*
 * mesh.cxx
 *
 *  Created on: 22 janv. 2012
 *      Authors: CDMATH
 */

#include "Mesh.hxx"
#include "Node.hxx"
#include "Cell.hxx"
#include "Face.hxx"

#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "CdmathException.hxx"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm> 
#include <string.h>

using namespace MEDCoupling;
using namespace std;

//----------------------------------------------------------------------
Mesh::Mesh( void )
//----------------------------------------------------------------------
{
	_mesh=NULL;
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
	delete [] _cells;
	delete [] _nodes;
	delete [] _faces;
	//for(int i=0; i< _faceGroups.size(); i++)
	//	_faceGroups[i]->decrRef();
	//	_nodeGroups[i]->decrRef();
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
    _isStructured=true;
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
	setMesh();
}

Mesh::Mesh( const MEDCoupling::MEDCouplingUMesh* mesh )
{
	_spaceDim=mesh->getSpaceDimension();
	_meshDim=mesh->getMeshDimension();
    _isStructured=false;
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

	_mesh=mesh->deepCopy();
	delete [] Box0 ;
	setMesh();
}

//----------------------------------------------------------------------
Mesh::Mesh( const Mesh& mesh )
//----------------------------------------------------------------------
{
	_spaceDim = mesh.getSpaceDimension() ;
	_meshDim = mesh.getMeshDimension() ;
    _name=mesh.getName();
    _xMax=mesh.getXMax();
    _yMin=mesh.getYMin();
    _yMax=mesh.getYMax();
    _zMin=mesh.getZMin();
    _zMax=mesh.getZMax();
    _isStructured=mesh.isStructured();
    if(_isStructured)
    {
        _nxyz = mesh.getCellGridStructure() ;
        _dxyz=mesh.getDXYZ();
    }
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
	_numberOfEdges = mesh.getNumberOfEdges();

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	_nodes   = new Node[_numberOfNodes] ;
	_faces   = new Face[_numberOfFaces] ;
	_cells   = new Cell[_numberOfCells] ;
    
    //memcpy(_nodes,mesh.getNodes(),sizeof(*mesh.getNodes())) ;
    //memcpy(_cells,mesh.getCells(),sizeof(*mesh.getCells())) ;
    //memcpy(_faces,mesh.getFaces(),sizeof(*mesh.getFaces())) ;

	for (int i=0;i<_numberOfNodes;i++)
		_nodes[i]=mesh.getNode(i);

	for (int i=0;i<_numberOfFaces;i++)
		_faces[i]=mesh.getFace(i);

	for (int i=0;i<_numberOfCells;i++)
		_cells[i]=mesh.getCell(i);

    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();

    _boundaryFaceIds=mesh.getBoundaryFaceIds();
    _boundaryNodeIds=mesh.getBoundaryNodeIds();

    _eltsTypes=mesh.getElementTypes();
    _eltsTypesNames=mesh.getElementTypesNames();
    
	MCAuto<MEDCouplingMesh> m1=mesh.getMEDCouplingMesh()->deepCopy();
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
	MEDFileUMesh *m=MEDFileUMesh::New(filename.c_str());//reads the first mesh encountered in the file, otherwise call New (const char *fileName, const char *mName, int dt=-1, int it=-1)
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
        _isStructured=true;
        _dxyz=structuredMesh->getDXYZ();
        _nxyz=structuredMesh->getCellGridStructure();
        double* Box0=new double[2*_spaceDim];
        structuredMesh->getBoundingBox(Box0);
    
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
    }
    else
        _isStructured=false;
    
	MEDCouplingUMesh*  mu = setMesh();
	setGroups(m, mu);

	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Mesh name= "<<m->getName()<<", mesh dim="<< _meshDim<< ", space dim="<< _spaceDim<< ", nb cells= "<<getNumberOfCells()<< ", nb nodes= "<<getNumberOfNodes()<<endl;

	m->decrRef();
	mu->decrRef();
}

void
Mesh::setGroupAtFaceByCoords(double x, double y, double z, double eps, std::string groupName)
{
	int nbFace=getNumberOfFaces();
	bool flag=false;
	for (int iface=0;iface<nbFace;iface++)
	{
		double FX=_faces[iface].x();
		double FY=_faces[iface].y();
		double FZ=_faces[iface].z();
		if (abs(FX-x)<eps && abs(FY-y)<eps && abs(FZ-z)<eps)
		{
			_faces[iface].setGroupName(groupName);
			std::vector< int > nodesID= _faces[iface].getNodesId();
			int nbNodes = _faces[iface].getNumberOfNodes();
			for(int inode=0 ; inode<nbNodes ; inode++)
				_nodes[nodesID[inode]].setGroupName(groupName);

			flag=true;
		}
	}
	if (flag)
    {
		_faceGroupNames.insert(_faceGroupNames.begin(),groupName);
		_nodeGroupNames.insert(_nodeGroupNames.begin(),groupName);
        //To do : update _faceGroups and _nodeGroups
    }
}

void
Mesh::setGroupAtPlan(double value, int direction, double eps, std::string groupName)
{
	int nbFace=getNumberOfFaces();
	bool flag=false;
	for (int iface=0;iface<nbFace;iface++)
	{
		double cord=_faces[iface].getBarryCenter()[direction];
		if (abs(cord-value)<eps)
		{
			_faces[iface].setGroupName(groupName);
			std::vector< int > nodesID= _faces[iface].getNodesId();
			int nbNodes = _faces[iface].getNumberOfNodes();
			for(int inode=0 ; inode<nbNodes ; inode++)
                {
				_nodes[nodesID[inode]].setGroupName(groupName);
                }

			flag=true;
		}
	}
	if (flag)
    {
		_faceGroupNames.insert(_faceGroupNames.begin(),groupName);
		_nodeGroupNames.insert(_nodeGroupNames.begin(),groupName);
        //To do : update _faceGroups, _nodeGroups
    }
}

void
Mesh::setBoundaryNodesFromFaces()
{
    for (int iface=0;iface<_boundaryFaceIds.size();iface++)
    {
        std::vector< int > nodesID= _faces[_boundaryFaceIds[iface]].getNodesId();
        int nbNodes = _faces[_boundaryFaceIds[iface]].getNumberOfNodes();
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
        
    double eps=1.E-10;

    for (int indexFace=0;indexFace<_boundaryFaceIds.size() ; indexFace++)
    {
        Face my_face=_faces[_boundaryFaceIds[indexFace]];
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
                Face face_i=_faces[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                if (   (abs(y-yi)<eps || abs(x-xi)<eps )// Case of a square geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi)<eps ) // Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<eps
                    && fabs(my_face.getXN()      + face_i.getXN())<eps
                    && fabs(my_face.getYN()      + face_i.getYN())<eps
                    && fabs(my_face.getZN()      + face_i.getZN())<eps )
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
                Face face_i=_faces[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                double zi=face_i.z();
                if ( ((abs(y-yi)<eps && abs(x-xi)<eps) || (abs(x-xi)<eps && abs(z-zi)<eps) || (abs(y-yi)<eps && abs(z-zi)<eps))// Case of a cube geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi) + abs(z+zi)<eps )// Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<eps
                    && fabs(my_face.getXN()      + face_i.getXN())<eps
                    && fabs(my_face.getYN()      + face_i.getYN())<eps
                    && fabs(my_face.getZN()      + face_i.getZN())<eps )
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
	if (!_faces[indexFace].isBorder())
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
	return _nodes[nodeid].isBorder();
}

bool
Mesh::isBorderFace(int faceid) const
{
	return _faces[faceid].isBorder();
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
Mesh::setGroups( const MEDFileUMesh* medmesh, MEDCouplingUMesh*  mu)
{
	//Searching for face groups
	vector<string> faceGroups=medmesh->getGroupsNames() ;

	for (unsigned int i=0;i<faceGroups.size();i++ )
	{
		string groupName=faceGroups[i];
		vector<int> nonEmptyGrp(medmesh->getGrpNonEmptyLevels(groupName));
		//We check if the group has a relative dimension equal to -1 
		//before call to the function getGroup(-1,groupName.c_str())
		vector<int>::iterator it = find(nonEmptyGrp.begin(), nonEmptyGrp.end(), -1);
		if (it != nonEmptyGrp.end())
		{
			cout<<"Boundary face group named "<< groupName << " found"<<endl;
			MEDCouplingUMesh *m=medmesh->getGroup(-1,groupName.c_str());
            mu->unPolyze();
			_faceGroups.insert(_faceGroups.begin(),m);
			_faceGroupNames.insert(_faceGroupNames.begin(),groupName);
			DataArrayDouble *baryCell = m->computeIsoBarycenterOfNodesPerCell();//computeCellCenterOfMass() ;
			const double *coorBary=baryCell->getConstPointer();

			int nbCellsSubMesh=m->getNumberOfCells();
			for (int ic(0), k(0); ic<nbCellsSubMesh; ic++, k+=_spaceDim)
			{
				vector<double> coorBaryXyz(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorBaryXyz[d] = coorBary[k+d];
				Point p1(coorBaryXyz[0],coorBaryXyz[1],coorBaryXyz[2]) ;

				int flag=0;
                /* double min=INFINITY;
                Point p3;
                int closestFaceNb;*/
				for (int iface=0;iface<_numberOfFaces;iface++ )
				{
					Point p2=_faces[iface].getBarryCenter();
                    /*if(groupName=="Wall" and p1.distance(p2)<min)
                    {
                        min=p1.distance(p2);
                        p3=p2;
                        closestFaceNb=iface;
                    }*/
					if(p1.distance(p2)<1.E-10)
					{
						_faces[iface].setGroupName(groupName);
						flag=1;
						break;
					}
				}
                /* if(groupName=="Wall" and min>1.E-10)
                    {
                        cout.precision(15);
                        cout<<"Subcell number "<< ic <<" Min distance to Wall = "<<min <<" p= "<<p1[0]<<" , "<<p1[1]<<" , "<<p1[2]<<endl;
                        cout<<" attained at face "<< closestFaceNb << " p_min= "<<p3[0]<<" , "<<p3[1]<<" , "<<p3[2]<<endl;
                    }*/
				if (flag==0)
					throw CdmathException("No face belonging to group " + groupName + " found");
			}
			baryCell->decrRef();
			//m->decrRef();
		}
	}

	//Searching for node groups
	vector<string> nodeGroups=medmesh->getGroupsOnSpecifiedLev(1) ;

	for (unsigned int i=0;i<nodeGroups.size();i++ )
	{
		string groupName=nodeGroups[i];
		DataArrayInt * nodeGroup=medmesh->getNodeGroupArr( groupName );
		const int *nodeids=nodeGroup->getConstPointer();

		if(nodeids!=NULL)
		{
			cout<<"Boundary node group named "<< groupName << " found"<<endl;

			_nodeGroups.insert(_nodeGroups.begin(),nodeGroup);
			_nodeGroupNames.insert(_nodeGroupNames.begin(),groupName);

			int nbNodesSubMesh=nodeGroup->getNumberOfTuples();//nodeGroup->getNbOfElems();

			DataArrayDouble *coo = mu->getCoords() ;
			const double *cood=coo->getConstPointer();

			for (int ic(0); ic<nbNodesSubMesh; ic++)
			{
				vector<double> coorP(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorP[d] = cood[nodeids[ic]*_spaceDim+d];
				Point p1(coorP[0],coorP[1],coorP[2]) ;

				int flag=0;
				for (int inode=0;inode<_numberOfNodes;inode++ )
				{
					Point p2=_nodes[inode].getPoint();
					if(p1.distance(p2)<1.E-10)
					{
						_nodes[inode].setGroupName(groupName);
						flag=1;
						break;
					}
				}
				if (flag==0)
					throw CdmathException("No node belonging to group " + groupName + " found");
			}
		}
	}
}

//----------------------------------------------------------------------
MEDCouplingUMesh* 
Mesh::setMesh( void )
//----------------------------------------------------------------------
{
	DataArrayInt *desc  = DataArrayInt::New();
	DataArrayInt *descI = DataArrayInt::New();
	DataArrayInt *revDesc  = DataArrayInt::New();
	DataArrayInt *revDescI = DataArrayInt::New();
	MEDCouplingUMesh* mu = _mesh->buildUnstructured();

	mu->unPolyze();
	/* Save the boundary mesh for later use*/
	_boundaryMesh = mu->computeSkin();
	
	MEDCouplingUMesh* mu2=mu->buildDescendingConnectivity2(desc,descI,revDesc,revDescI);//mesh of dimension N-1 containing the cell interfaces
    
    const int *tmp = desc->getConstPointer();
    const int *tmpI=descI->getConstPointer();

	const int *tmpA =revDesc->getConstPointer();
	const int *tmpAI=revDescI->getConstPointer();

    //const int *work=tmp+tmpI[id];//corresponds to buildDescendingConnectivity

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
	const double *coorBary=baryCell->getConstPointer();

	MEDCouplingFieldDouble* fields=mu->getMeasureField(true);
	DataArrayDouble *surface = fields->getArray();
	const double *surf=surface->getConstPointer();

	DataArrayDouble *coo = mu->getCoords() ;
	const double    *cood=coo->getConstPointer();

	DataArrayInt *revNode =DataArrayInt::New();
	DataArrayInt *revNodeI=DataArrayInt::New();
	mu->getReverseNodalConnectivity(revNode,revNodeI) ;
	const int *tmpN =revNode->getConstPointer();
	const int *tmpNI=revNodeI->getConstPointer();

	DataArrayInt *revCell =DataArrayInt::New();
	DataArrayInt *revCellI=DataArrayInt::New();
	mu2->getReverseNodalConnectivity(revCell,revCellI) ;
	const int *tmpC =revCell->getConstPointer();
	const int *tmpCI=revCellI->getConstPointer();

	const DataArrayInt *nodal  = mu2->getNodalConnectivity() ;
	const DataArrayInt *nodalI = mu2->getNodalConnectivityIndex() ;
	const int *tmpNE =nodal->getConstPointer();
	const int *tmpNEI=nodalI->getConstPointer();

	_numberOfCells = mu->getNumberOfCells() ;
	_cells      = new Cell[_numberOfCells] ;

	_numberOfNodes = mu->getNumberOfNodes() ;
	_nodes      = new Node[_numberOfNodes] ;

	_numberOfFaces = mu2->getNumberOfCells();
	_faces       = new Face[_numberOfFaces] ;

    _indexFacePeriodicSet=false;

    //Definition used if _meshDim =3 to determine the edges
    DataArrayInt *desc2 =DataArrayInt::New();
    DataArrayInt *descI2=DataArrayInt::New();
    DataArrayInt *revDesc2 =DataArrayInt::New();
    DataArrayInt *revDescI2=DataArrayInt::New();
    DataArrayInt *revNode2 =DataArrayInt::New();
    DataArrayInt *revNodeI2=DataArrayInt::New();
    const int *tmpN2 ;
    const int *tmpNI2;
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
	if (_spaceDim == 1)
	{
		for( int id=0;id<_numberOfCells;id++ )
		{
			const int *work=tmp+tmpI[id];
            int nbFaces=tmpI[id+1]-tmpI[id];
			
			int nbVertices=mu->getNumberOfNodesInCell(id) ;
            
			Cell ci( nbVertices, nbFaces, surf[id], Point(coorBary[id], 0.0, 0.0) ) ;

			std::vector<int> nodeIdsOfCell ;
			mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
			for( int el=0;el<nbVertices;el++ )
            {
				ci.addNodeId(el,nodeIdsOfCell[el]) ;
                ci.addFaceId(el,nodeIdsOfCell[el]) ;
            }
            
            //This assumes that in a 1D mesh the cell numbers form a consecutive sequence 
            //going either from left to right (xn=-1) or right to left (xn=1)
			double xn = (cood[nodeIdsOfCell[nbVertices-1]] - cood[nodeIdsOfCell[0]] > 0.0) ? -1.0 : 1.0;

			for( int el=0;el<nbFaces;el++ )
			{
				ci.addNormalVector(el,xn,0.0,0.0) ;
				int indexFace=abs(work[el])-1;
                ci.addFaceId(el,indexFace) ;
				xn = - xn;
			}
			_cells[id] = ci ;
		}

		for( int id(0), k(0); id<_numberOfNodes; id++, k+=_spaceDim)
		{
			Point p(cood[k], 0.0, 0.0) ;
			const int *workc=tmpN+tmpNI[id];
			int nbCells=tmpNI[id+1]-tmpNI[id];

			const int *workf=tmpC+tmpCI[id];
			int nbFaces=tmpCI[id+1]-tmpCI[id];
            const int *workn=tmpA+tmpAI[id];
            int nbNeighbourNodes=tmpAI[id+1]-tmpAI[id];
			Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;

			for( int el=0;el<nbCells;el++ )
				vi.addCellId(el,workc[el]) ;
			for( int el=0;el<nbFaces;el++ )
				vi.addFaceId(el,workf[el]) ;
			for( int el=0;el<nbNeighbourNodes;el++ )
				vi.addNeighbourNodeId(el,workn[el]) ;
			_nodes[id] = vi ;
		}
        _boundaryNodeIds.push_back(0);
        _boundaryNodeIds.push_back(_numberOfNodes-1);

		for(int id(0), k(0); id<_numberOfFaces; id++, k+=_spaceDim)
		{
			Point p(cood[k], 0.0, 0.0) ;
			const int *workc=tmpA+tmpAI[id];
			int nbCells=tmpAI[id+1]-tmpAI[id];

			const int *workv=tmpNE+tmpNEI[id]+1;
			Face fi( 1, nbCells, 1.0, p, 1.0, 0.0, 0.0) ;
			fi.addNodeId(0,workv[0]) ;

			for(int idCell=0; idCell<nbCells; idCell++)
				fi.addCellId(idCell,workc[idCell]) ;

			_faces[id] = fi ;
		}
        _boundaryFaceIds.push_back(0);
        _boundaryFaceIds.push_back(_numberOfFaces-1);
	}
	else if(_spaceDim==2  || _spaceDim==3)
	{
		DataArrayDouble *barySeg = mu2->computeIsoBarycenterOfNodesPerCell();//computeCellCenterOfMass() ;
		const double *coorBarySeg=barySeg->getConstPointer();

		MEDCouplingFieldDouble* fieldn;
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
            const int *work=tmp+tmpI[id];      
			int nbFaces=tmpI[id+1]-tmpI[id];
            
			int nbVertices=mu->getNumberOfNodesInCell(id) ;

			vector<double> coorBaryXyz(3,0);
			for (int d=0; d<_spaceDim; d++)
				coorBaryXyz[d] = coorBary[k+d];

			Point p(coorBaryXyz[0],coorBaryXyz[1],coorBaryXyz[2]) ;
			Cell ci( nbVertices, nbFaces, surf[id], p ) ;

			/* Filling cell nodes */
			std::vector<int> nodeIdsOfCell ;
			mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
			for( int el=0;el<nbVertices;el++ )
				ci.addNodeId(el,nodeIdsOfCell[el]) ;

			/* Filling cell faces */
			if(_spaceDim==_meshDim)//use the normal field generated by buildOrthogonalField()
				for( int el=0;el<nbFaces;el++ )
				{
                    int faceIndex=(abs(work[el])-1);//=work[el] since Fortran type numbering was used, and negative sign means anticlockwise numbering
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
			{
				if(_meshDim==1)//we know in this case there are only two faces around the cell id, each face is composed of a single node
				{//work[0]= first face global number, work[1]= second face global number
                    int indexFace0=abs(work[0])-1;//=work[0] since Fortran type numbering was used, and negative sign means anticlockwise numbering
                    int indexFace1=abs(work[1])-1;//=work[1] since Fortran type numbering was used, and negative sign means anticlockwise numbering
					int idNodeA=(tmpNE+tmpNEI[indexFace0]+1)[0];//global number of the first  face node work[0]=(abs(work[0])-1)
					int idNodeB=(tmpNE+tmpNEI[indexFace1]+1)[0];//global number of the second face node work[1]=(abs(work[1])-1)
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
						const int *workv=tmpNE+tmpNEI[faceIndex]+1;
						int nbNodes= tmpNEI[faceIndex+1]-tmpNEI[faceIndex]-1;
						if(nbNodes!=2)//We want to compute the normal to a straight line, not a curved interface composed of more thant 2 points
						{
							cout<<"Mesh name "<< mu->getName()<< " space dim= "<< _spaceDim <<" mesh dim= "<< _meshDim <<endl;
							cout<<"For cell id "<<id<<" and local face number "<<el<<", the number of nodes is "<< nbNodes<< ", total number of faces is "<< nbFaces <<endl;
							throw CdmathException("Mesh::setMesh number of nodes around a face should be 2");
						}

						int idNodeA=workv[0];
						int idNodeB=workv[1];
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
			_cells[id] = ci ;
		}

		MEDCouplingFieldDouble* fieldl=mu2->getMeasureField(true);
		DataArrayDouble *longueur = fieldl->getArray();
		const double *lon=longueur->getConstPointer();

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
			vector<double> coorBarySegXyz(3,0);
			for (int d=0; d<_spaceDim; d++)
				coorBarySegXyz[d] = coorBarySeg[k+d];
			Point p(coorBarySegXyz[0],coorBarySegXyz[1],coorBarySegXyz[2]) ;
			const int *workc=tmpA+tmpAI[id];
			int nbCells=tmpAI[id+1]-tmpAI[id];
            
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
                
			const int *workv=tmpNE+tmpNEI[id]+1;
			int nbNodes= tmpNEI[id+1]-tmpNEI[id]-1;

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
            
			_faces[id] = fi ;
		}

		/*Building mesh nodes, should be done after building mesh faces in order to detect boundary nodes*/
		for(int id(0), k(0); id<_numberOfNodes; id++, k+=_spaceDim)
		{
			vector<double> coorP(3,0);
			for (int d=0; d<_spaceDim; d++)
				coorP[d] = cood[k+d];
			Point p(coorP[0],coorP[1],coorP[2]) ;

			const int *workc=tmpN+tmpNI[id];
			int nbCells=tmpNI[id+1]-tmpNI[id];
			const int *workf=tmpC+tmpCI[id];
			int nbFaces=tmpCI[id+1]-tmpCI[id];
			const int *workn;
			int nbNeighbourNodes;
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
            bool isBorder=false;
			for( int el=0;el<nbFaces;el++ )
            {
				vi.addFaceId(el,workf[el],_faces[workf[el]].isBorder()) ;
                isBorder= isBorder || _faces[workf[el]].isBorder();
            }
            if(isBorder)
                _boundaryNodeIds.push_back(id);
			_nodes[id] = vi ;
		}

		if(_spaceDim==_meshDim)
			fieldn->decrRef();
		fieldl->decrRef();
		barySeg->decrRef();
	}
	else
		throw CdmathException("Mesh::setMesh space dimension should be 1, 2 or 3");

    //definition of the bounding box for unstructured meshes
    if(!_isStructured)//Structured case is treated in function readMeshMed
    {
        double Box0[2*_spaceDim];
        coo->getMinMaxPerComponent(Box0);

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
    }

    //Set boundary groups
	_faceGroupNames.push_back("Boundary");
    _nodeGroupNames.push_back("Boundary");
    
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

    if (_meshDim == 3)
    {
        revNode2->decrRef();
        revNodeI2->decrRef();
        desc2->decrRef();
        descI2->decrRef();
        revDesc2->decrRef();
        revDescI2->decrRef();
        mu3->decrRef();
    }
	
    return mu;
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
    _isStructured = true;
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

	DataArrayInt *desc=DataArrayInt::New();
	DataArrayInt *descI=DataArrayInt::New();
	DataArrayInt *revDesc=DataArrayInt::New();
	DataArrayInt *revDescI=DataArrayInt::New();
	MEDCouplingUMesh* mu=_mesh->buildUnstructured();
	MEDCouplingUMesh *mu2=mu->buildDescendingConnectivity(desc,descI,revDesc,revDescI);

	const int *tmp=desc->getConstPointer();
	const int *tmpI=descI->getConstPointer();

	const int *tmpA =revDesc->getConstPointer();
	const int *tmpAI=revDescI->getConstPointer();

	_eltsTypes=mu->getAllGeoTypesSorted();

	DataArrayDouble *baryCell = mu->computeCellCenterOfMass() ;
	const double *coorBary=baryCell->getConstPointer();

	_numberOfCells = _mesh->getNumberOfCells() ;
	_cells    = new Cell[_numberOfCells] ;

	_numberOfNodes = mu->getNumberOfNodes() ;
	_nodes    = new Node[_numberOfNodes] ;

	_numberOfFaces = _numberOfNodes;
	_faces    = new Face[_numberOfFaces] ;

    _numberOfEdges = _numberOfCells;
    
	MEDCouplingFieldDouble* fieldl=mu->getMeasureField(true);
	DataArrayDouble *longueur = fieldl->getArray();
	const double *lon=longueur->getConstPointer();

	DataArrayDouble *coo = mu->getCoords() ;
	const double *cood=coo->getConstPointer();

	int comp=0;
	for( int id=0;id<_numberOfCells;id++ )
	{
		int nbVertices=mu->getNumberOfNodesInCell(id) ;
		Point p(coorBary[id],0.0,0.0) ;
		Cell ci( nbVertices, nbVertices, lon[id], p ) ;

		std::vector<int> nodeIdsOfCell ;
		mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
		for( int el=0;el<nbVertices;el++ )
		{
			ci.addNodeId(el,nodeIdsOfCell[el]) ;
			ci.addFaceId(el,nodeIdsOfCell[el]) ;
		}

        double xn = (cood[nodeIdsOfCell[nbVertices-1]] - cood[nodeIdsOfCell[0]] > 0.0) ? -1.0 : 1.0;

        int nbFaces=tmpI[id+1]-tmpI[id];
        const int *work=tmp+tmpI[id];
		
        for( int el=0;el<nbFaces;el++ )
		{
			ci.addNormalVector(el,xn,0.0,0.0) ;
			ci.addFaceId(el,work[el]) ;
			xn = - xn;
		}

		_cells[id] = ci ;

		comp=comp+2;
	}

    //Suppress the following since tmpN=tmpA
	DataArrayInt *revNode=DataArrayInt::New();
	DataArrayInt *revNodeI=DataArrayInt::New();
	mu->getReverseNodalConnectivity(revNode,revNodeI) ;
	const int *tmpN=revNode->getConstPointer();
	const int *tmpNI=revNodeI->getConstPointer();

	for( int id=0;id<_numberOfNodes;id++ )
	{
		std::vector<double> coo ;
		mu->getCoordinatesOfNode(id,coo);
		Point p(coo[0],0.0,0.0) ;
		const int *workc=tmpN+tmpNI[id];
		int nbCells=tmpNI[id+1]-tmpNI[id];
		int nbFaces=1;
        const int *workn=tmpA+tmpAI[id];
        int nbNeighbourNodes=tmpAI[id+1]-tmpAI[id];
        
		Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;
        for( int el=0;el<nbCells;el++ )
			vi.addCellId(el,workc[el]) ;
		for( int el=0;el<nbFaces;el++ )
			vi.addFaceId(el,id) ;
        for( int el=0;el<nbNeighbourNodes;el++ )
			vi.addNeighbourNodeId(el,workn[el]) ;
		_nodes[id] = vi ;


		int nbVertices=1;
		Face fi( nbVertices, nbCells, 1.0, p, 1., 0., 0. ) ;
        for( int el=0;el<nbVertices;el++ )
			fi.addNodeId(el,id) ;

		for( int el=0;el<nbCells;el++ )
			fi.addCellId(el,workc[el]) ;
		_faces[id] = fi ;
	}

    //Set boundary groups
	_faceGroupNames.push_back("Boundary");
    _nodeGroupNames.push_back("Boundary");
    _boundaryFaceIds.push_back(0);
    _boundaryFaceIds.push_back(_numberOfFaces-1);
    _boundaryNodeIds.push_back(0);
    _boundaryNodeIds.push_back(_numberOfNodes-1);

	fieldl->decrRef();
	baryCell->decrRef();
	desc->decrRef();
	descI->decrRef();
	revDesc->decrRef();
	revDescI->decrRef();
	revNode->decrRef();
	revNodeI->decrRef();
	mu2->decrRef();
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
	_xMin=points[0];
	_xMax=points[nx-1];
	_yMin=0.;
	_yMax=0.;
	_zMin=0.;
	_zMax=0.;

    _isStructured = false;
    _indexFacePeriodicSet=false;
    
    MEDCouplingUMesh * mesh1d = MEDCouplingUMesh::New(meshName, 1);
    mesh1d->allocateCells(nx - 1);
    double * coords = new double[nx];
    int * nodal_con = new int[2];
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

    delete [] coords, nodal_con;
    coords_arr->decrRef();

    _mesh=mesh1d;
    
	DataArrayInt *desc=DataArrayInt::New();
	DataArrayInt *descI=DataArrayInt::New();
	DataArrayInt *revDesc=DataArrayInt::New();
	DataArrayInt *revDescI=DataArrayInt::New();
    MEDCouplingUMesh* mu=_mesh->buildUnstructured();
	MEDCouplingUMesh *mu2=mu->buildDescendingConnectivity(desc,descI,revDesc,revDescI);

	DataArrayDouble *baryCell = mu->computeCellCenterOfMass() ;
	const double *coorBary=baryCell->getConstPointer();

	_numberOfCells = _mesh->getNumberOfCells() ;
	_cells    = new Cell[_numberOfCells] ;

    _numberOfEdges = _numberOfCells;

	_eltsTypes=mu->getAllGeoTypesSorted();

	MEDCouplingFieldDouble* fieldl=mu->getMeasureField(true);
	DataArrayDouble *longueur = fieldl->getArray();
	const double *lon=longueur->getConstPointer();

	int comp=0;
	for( int id=0;id<_numberOfCells;id++ )
	{
		int nbVertices=mu->getNumberOfNodesInCell(id) ;
		Point p(coorBary[id],0.0,0.0) ;
		Cell ci( nbVertices, nbVertices, lon[id], p ) ;

		std::vector<int> nodeIdsOfCell ;
		mu->getNodeIdsOfCell(id,nodeIdsOfCell) ;
		for( int el=0;el<nbVertices;el++ )
		{
			ci.addNodeId(el,nodeIdsOfCell[el]) ;
			ci.addFaceId(el,nodeIdsOfCell[el]) ;
		}
		_cells[id] = ci ;
		comp=comp+2;
	}


    //Suppress the following since tmpN=tmpA
	DataArrayInt *revNode=DataArrayInt::New();
	DataArrayInt *revNodeI=DataArrayInt::New();
	mu->getReverseNodalConnectivity(revNode,revNodeI) ;
	const int *tmpN=revNode->getConstPointer();
	const int *tmpNI=revNodeI->getConstPointer();

	_numberOfNodes = mu->getNumberOfNodes() ;
	_nodes    = new Node[_numberOfNodes] ;
	_numberOfFaces = _numberOfNodes;
	_faces    = new Face[_numberOfFaces] ;

	for( int id=0;id<_numberOfNodes;id++ )
	{
		std::vector<double> coo ;
		mu->getCoordinatesOfNode(id,coo);
		Point p(coo[0],0.0,0.0) ;
		const int *workc=tmpN+tmpNI[id];
		int nbCells=tmpNI[id+1]-tmpNI[id];
		int nbFaces=1;
        const int *workn=tmpN+tmpNI[id];
        int nbNeighbourNodes=tmpNI[id+1]-tmpNI[id];
		Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;
		int nbVertices=1;
		/* provisoire !!!!!!!!!!!!*/
		//        Point pf(0.0,0.0,0.0) ;
		Face fi( nbVertices, nbCells, 0.0, p, 0., 0., 0. ) ;

		for( int el=0;el<nbCells;el++ )
			vi.addCellId(el,workc[el]) ;
		for( int el=0;el<nbFaces;el++ )
			vi.addFaceId(el,id) ;
        for( int el=0;el<nbNeighbourNodes;el++ )
			vi.addNeighbourNodeId(el,workn[el]) ;
		_nodes[id] = vi ;

		for( int el=0;el<nbVertices;el++ )
			fi.addNodeId(el,id) ;

		for( int el=0;el<nbCells;el++ )
			fi.addCellId(el,workc[el]) ;
		_faces[id] = fi ;
	}
    //Set boundary groups
	_faceGroupNames.push_back("Boundary");
    _nodeGroupNames.push_back("Boundary");
    _boundaryFaceIds.push_back(0);
    _boundaryFaceIds.push_back(_numberOfFaces-1);
    _boundaryNodeIds.push_back(0);
    _boundaryNodeIds.push_back(_numberOfNodes-1);

	fieldl->decrRef();
	baryCell->decrRef();
	desc->decrRef();
	descI->decrRef();
	revDesc->decrRef();
	revDescI->decrRef();
	revNode->decrRef();
	revNodeI->decrRef();
	mu2->decrRef();
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
    _isStructured = true;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_dxyz[1]=dy;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr   = new double[_spaceDim];
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

    if(split_to_triangles_policy==0 || split_to_triangles_policy==1)
        {
            _mesh=_mesh->buildUnstructured();
            _mesh->simplexize(split_to_triangles_policy);
        }
    else if (split_to_triangles_policy != -1)
        {
            cout<< "split_to_triangles_policy = "<< split_to_triangles_policy << endl;
            throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : Unknown splitting policy");
        }
        
	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
    
	setMesh();
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

    _isStructured = true;
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

    if( split_to_tetrahedra_policy == 0 )
        {
            _mesh=_mesh->buildUnstructured();
            _mesh->simplexize(INTERP_KERNEL::PLANAR_FACE_5);
        }
    else if( split_to_tetrahedra_policy == 1 )
        {
            _mesh=_mesh->buildUnstructured();
            _mesh->simplexize(INTERP_KERNEL::PLANAR_FACE_6);
        }
    else if ( split_to_tetrahedra_policy != -1 )
        {
            cout<< "split_to_tetrahedra_policy = "<< split_to_tetrahedra_policy << endl;
            throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : splitting policy value should be 0 or 1");
        }

	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
    
	setMesh();
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
    if(!_isStructured)
		throw CdmathException("std::vector<double> Mesh::getDXYZ() : dx,dy and dz are defined only for structured meshes !");

	return _dxyz;
}

std::vector<int>
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
	if(_spaceDim < 2)
		throw CdmathException("int double& Field::operator[ielem] : Ny is not defined in dimension < 2!");

	return _nxyz[1];
}

//----------------------------------------------------------------------
int
Mesh::getNz( void )  const
//----------------------------------------------------------------------
{
    if(!_isStructured)
		throw CdmathException("int Mesh::getNz( void ) : Nz is defined only for structured meshes !");
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
Face*
Mesh::getFaces ( void )  const
//----------------------------------------------------------------------
{
	return _faces ;
}

//----------------------------------------------------------------------
Cell*
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

//----------------------------------------------------------------------
Node*
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

vector<MEDCoupling::MEDCouplingUMesh *>
Mesh::getFaceGroups( void )  const
{
	return _faceGroups;
}

vector<string>
Mesh::getNameOfNodeGroups( void )  const
{
	return _nodeGroupNames;
}

vector<MEDCoupling::DataArrayInt *>
Mesh::getNodeGroups( void )  const
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
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
	_numberOfEdges = mesh.getNumberOfEdges();
    
    _isStructured = mesh.isStructured();
    if(_isStructured)
    {
        _nxyz = mesh.getCellGridStructure() ;
        _dxyz=mesh.getDXYZ();
        _xMin=mesh.getXMin();
        _xMax=mesh.getXMax();
        _yMin=mesh.getYMin();
        _yMax=mesh.getYMax();
        _zMin=mesh.getZMin();
        _zMax=mesh.getZMax();
    }
	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	if (_nodes)
	{
		delete [] _nodes ;
		_nodes=NULL;
	}
	if (_faces)
	{
		delete [] _faces ;
		_faces=NULL;
	}
	if (_cells)
	{
		delete [] _cells ;
		_cells=NULL;
	}

	_nodes   = new Node[_numberOfNodes] ;
	_faces   = new Face[_numberOfFaces] ;
	_cells   = new Cell[_numberOfCells] ;

	for (int i=0;i<_numberOfNodes;i++)
		_nodes[i]=mesh.getNode(i);

	for (int i=0;i<_numberOfFaces;i++)
		_faces[i]=mesh.getFace(i);

	for (int i=0;i<_numberOfCells;i++)
		_cells[i]=mesh.getCell(i);

    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();
    
    _boundaryFaceIds=mesh.getBoundaryFaceIds();
    _boundaryNodeIds=mesh.getBoundaryNodeIds();

    _eltsTypes=mesh.getElementTypes();
    _eltsTypesNames=mesh.getElementTypesNames();

	MCAuto<MEDCouplingMesh> m1=mesh.getMEDCouplingMesh()->deepCopy();
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
	return Mesh(_boundaryMesh);
}

int 
Mesh::getMaxNbNeighbours(EntityType type) const
{
    double result=0;
    
    if (type==CELLS)
	{
        for(int i=0; i<_numberOfCells; i++)
            if(result < _cells[i].getNumberOfFaces())
                result=_cells[i].getNumberOfFaces();
	}
    else if(type==NODES)
	{
        for(int i=0; i<_numberOfNodes; i++)
            if(result < _nodes[i].getNumberOfEdges())
                result=_nodes[i].getNumberOfEdges();
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
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
Mesh::writeMED ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	MEDCouplingUMesh* mu=_mesh->buildUnstructured();
	MEDCoupling::WriteUMesh(fname.c_str(),mu,true);

	//MEDFileUMesh meshMEDFile;
	//meshMEDFile.setMeshAtLevel(0,mu);
	//for(int i=0; i< _faceGroups.size(); i++)
	//meshMEDFile.setMeshAtLevel(-1,_faceGroups[i]);
	//if (fromScratch)
	//MEDCoupling::meshMEDFile.write(fname.c_str(),2)	;
	//else
	//MEDCoupling::meshMEDFile.write(fname.c_str(),1)	;


	mu->decrRef();
}

std::vector< int > 
Mesh::getGroupFaceIds(std::string groupName) const
{
    if( std::find(_faceGroupNames.begin(),_faceGroupNames.end(),groupName) == _faceGroupNames.end() )
    {
        cout<<"Mesh::getGroupFaceIds Error : face group " << groupName << " does not exist"<<endl;
        throw CdmathException("Required face group does not exist");
    }
    else
    {
        std::vector< int > result(0);
        for(int i=0; i<_boundaryFaceIds.size(); i++)
        {
            vector< std::string > groupNames = getFace(_boundaryFaceIds[i]).getGroupNames();
            if( std::find(groupNames.begin(), groupNames.end(),groupName) != groupNames.end() )
                result.push_back(_boundaryFaceIds[i]);
        }
        return result;
    }
}

std::vector< int > 
Mesh::getGroupNodeIds(std::string groupName) const
{
    if( std::find(_nodeGroupNames.begin(),_nodeGroupNames.end(),groupName) == _nodeGroupNames.end() )
    {
        cout<<"Mesh::getGroupNodeIds Error : node group " << groupName << " does not exist"<<endl;
        throw CdmathException("Required node group does not exist");
    }
    else
    {
        std::vector< int > result(0);
        for(int i=0; i<_boundaryFaceIds.size(); i++)
        {
            vector< std::string > groupNames = getNode(_boundaryFaceIds[i]).getGroupNames();
            if( std::find(groupNames.begin(), groupNames.end(),groupName) != groupNames.end() )
                result.push_back(_boundaryFaceIds[i]);
        }
        return result;
    }
}

void 
Mesh::setFaceGroupByIds(std::vector< int > faceIds, std::string groupName) 
{
    for(int i=0; i< faceIds.size(); i++)
        getFace(faceIds[i]).setGroupName(groupName);
}

void 
Mesh::setNodeGroupByIds(std::vector< int > nodeIds, std::string groupName) 
{
    for(int i=0; i< nodeIds.size(); i++)
        getNode(nodeIds[i]).setGroupName(groupName);
}

