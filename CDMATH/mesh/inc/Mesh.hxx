/*
 * mesh.hxx
 *
 *  Created on: 22 janv. 2012
 *      Authors: CDMAT
 */

#ifndef MESH_HXX_
#define MESH_HXX_

#include <memory>

#include <MCAuto.hxx>
#include "NormalizedGeometricTypes"

#include "MEDCouplingUMesh.hxx"

/**
 * Mesh class is defined by
 * - case 1: file name of mesh med file (general unstructured)
 * - case 2: 1D cartesian, xmin and xmax and number of cells
 * - case 3: 2D cartesian, xmin, xmax, ymin and ymax and numbers of cells in x direction and y direction
 * - case 4: 3D cartesian, xmin, xmax, ymin, ymax, zmin and zmax and numbers of cells in x direction, y direction and z direction
 * - case 5: 2D regular triangular mesh
 * - case 6: 3D regular hexahedral mesh
 * - case 7: 1D unstructured
 */

namespace MEDCoupling
{
class MEDFileMesh;
class MEDFileUMesh;
class MEDCouplingMesh;
class DataArrayIdType;
}

class Node;
class Cell;
class Face;

typedef enum
  {
    CELLS = 11,
    NODES = 7,
    FACES = 5,
    GAUSS_PT =3
  } EntityType;

#include <vector>
#include <string>
#include <map>

class Mesh
{

public: //----------------------------------------------------------------
	/**
	 * \brief default constructor
	 */
	Mesh ( void ) ;

	/**
	 * \brief constructor with data from a medcoupling mesh
	 * @param medcoupling mesh 
	 */
	Mesh( MEDCoupling::MCAuto<const MEDCoupling::MEDCouplingMesh> mesh ) ;

	/**
	 * \brief constructor with data to load a general unstructured mesh
	 * @param filename name of mesh file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	Mesh ( const std::string filename, const std::string & meshName="" , int meshLevel=0) ;

	/**
	 * \brief constructor with data for a regular 1D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param nx : Number of cells in x direction
	 */
	Mesh( double xmin, double xmax, int nx, std::string meshName="MESH1D_Regular_Grid" ) ;

	/**
	 * \brief constructor with data for an unstructured 1D mesh
	 * @param points : abscissas of the mesh nodes
	 */
	Mesh( std::vector<double> points, std::string meshName="MESH1D_unstructured" ) ;

	/**
	 * \brief constructor with data for a regular 2D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
     * @param split_to_triangles_policy : each rectangle will be split into 2 triangles with orientation of the cut depending if value is 0 or 1
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, int split_to_triangles_policy=-1, std::string meshName="MESH2D_Regular_Rectangle_Grid") ;

	/**
	 * \brief constructor with data for a regular 3D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param zmin : minimum z
	 * @param zmax : maximum z
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 * @param nz : Number of cells in z direction
     * @param split_to_tetrahedra_policy : each cuboid will be split into 5 tetrahedra if value is 0 or 6 tetrahedra if the value is 1
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, int split_to_tetrahedra_policy=-1, std::string meshName="MESH3D_Regular_Cuboid_Grid") ;

	/**
	 * \brief constructor with data
	 * @param filename : file name of mesh med file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	void readMeshMed( const std::string filename, const std::string & meshName="" , int meshLevel=0) ;

	/**
	 * \brief constructor by copy
	 * @param mesh : The Mesh object to be copied
	 */
	Mesh ( const Mesh & mesh ) ;

	/**
	 * \brief destructor
	 */
	~Mesh( void ) ;

	/**
	 * \brief return mesh name
	 * @return _name
	 */
	std::string getName( void ) const ;

	/**
	 * \brief return Space dimension
	 * @return _spaceDim
	 */
	int getSpaceDimension( void ) const ;

	/**
	 * \brief return Mesh dimension
	 * @return _meshDim
	 */
	int getMeshDimension( void ) const ;

	/**
	 * \brief return the number of nodes in this mesh
	 * @return _numberOfNodes
	 */
	int getNumberOfNodes ( void )  const ;

	/**
	 * \brief return the number of faces in this mesh
	 * @return _numberOfFaces
	 */
	int getNumberOfFaces ( void )  const ;

	/**
	 * \brief return the number of cells in this mesh
	 * @return _numberOfCells
	 */
	int getNumberOfCells ( void )  const ;

	/**
	 * \brief return the number of edges in this mesh
	 * @return _numberOfEdges
	 */
	int getNumberOfEdges ( void )  const ;

	/**
	 * \brief return the cell i in this mesh
	 * @return _cells[i]
	 */
	Cell& getCell ( int i )  const ;

	/**
	 * return The face i in this mesh
	 * @return _faces[i]
	 */
	Face& getFace ( int i )  const ;

	/**
	 * \brief return The node i in this mesh
	 * @return _nodes[i]
	 */
	Node& getNode ( int i )  const ;

	/**
	 * \brief return number of cell in x direction (structured mesh)
	 * @return _nX
	 */
	int getNx( void )  const ;

	/**
	 * \brief return number of cell in y direction (structured mesh)
	 * @return _nY
	 */
	int getNy( void )  const ;

	/**
	 * \brief return number of cell in z direction (structured mesh)
	 * @return _nZ
	 */
	int getNz( void )  const ;

	double getXMin( void )  const ;

	double getXMax( void )  const ;

	double getYMin( void )  const ;

	double getYMax( void )  const ;

	double getZMin( void )  const ;

	double getZMax( void )  const ;

	std::vector<mcIdType> getCellGridStructure() const;// for structured meshes

	/**
	 * \brief overload operator =
	 * @param mesh : The Mesh object to be copied
	 */
	const Mesh& operator= ( const Mesh& mesh ) ;

	/**
	 * \brief return the mesh MEDCoupling
	 * @return _mesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> getMEDCouplingMesh ( void )  const ;

	/**
	 * \brief return the skin surrounding the mesh
	 */
	Mesh getBoundaryMesh ( void )  const ;

	/**
	 * \brief return the skin surrounding the mesh (MEDCouplingMesh)
	 */
	MEDCoupling::MEDCouplingUMesh * getBoundaryMEDCouplingMesh ( void )  const ;

	/**
	 * \brief return a group surrounding the mesh
	 */
	Mesh getBoundaryGroupMesh ( std::string groupName, int nth_group_match = 0 )  const ;

	/**
	 * \brief return the list of face group names
	 * @return _faceGroupNames
	 */
	std::vector<std::string> getNameOfFaceGroups( void )  const ;

	/**
	 * \brief return the list of node group names
	 * @return _nodeGroupNames
	 */
	std::vector<std::string> getNameOfNodeGroups( void )  const ;

	/**
	 * \brief return the list of face groups ids
	 * @return _faceGroupsIds
	 */
	std::vector< std::vector<int> > getFaceGroups( void )  const ;

	/**
	 * \brief return the list of node groups Ids
	 * @return _nodeGroupsIds
	 */
	std::vector< std::vector<int> > getNodeGroups( void )  const ;

	/**
	 * \brief return the list of face groups
	 * @return _faceGroups
	 */
	std::vector<MEDCoupling::MEDCouplingUMesh *> getMEDCouplingFaceGroups( void )  const ;

	/**
	 * \brief return the list of node groups
	 * @return _nodeGroups
	 */
	std::vector<MEDCoupling::DataArrayIdType *> getMEDCouplingNodeGroups( void )  const ;

    /**
     * \brief Functions to extract boundary nodes and faces Ids
     */
     /**
      *  \brief return the list of boundary faces Ids
      *  @return _boundaryFaceIds
      */
    std::vector< int > getBoundaryFaceIds() const;
    /**
     * \brief list of boundary nodes Ids
     * @return _boundaryNodeIds
     */
    std::vector< int > getBoundaryNodeIds() const;
    /**
     * \brief Functions to extract group nodes and faces ids
     */
     /** 
      * @return list of face group Ids
      */
    std::vector< int > getFaceGroupIds(std::string groupName, bool isBoundaryGroup=true) const;
    /**
     * @return list of node group Ids
     * */
    std::vector< int > getNodeGroupIds(std::string groupName, bool isBoundaryGroup=true) const;
 
	/**
	 * \brief write mesh in the VTK format
	 */
	void writeVTK ( const std::string fileName ) const ;

	/**
	 * \brief write mesh in the MED format
	 */
	void writeMED ( const std::string fileName, bool fromScratch = true ) const;

	void setGroupAtPlan(double value, int direction, double eps, std::string groupName, bool isBoundaryGroup=true) ;

	void setGroupAtFaceByCoords(double x, double y, double z, double eps, std::string groupName, bool isBoundaryGroup=true) ;

	void setGroupAtNodeByCoords(double x, double y, double z, double eps, std::string groupName, bool isBoundaryGroup=true) ;

	void setFaceGroupByIds(std::vector< int > faceIds, std::string groupName) ;

	void setNodeGroupByIds(std::vector< int > nodeIds, std::string groupName) ;

	/*
     * Functions to manage periodic boundary condition in square/cubic geometries 
     */
    void setPeriodicFaces(bool check_groups= false, bool use_central_inversion=false) ;
    int getIndexFacePeriodic(int indexFace, bool check_groups= false, bool use_central_inversion=false);
    void setBoundaryNodesFromFaces();
    std::map<int,int> getIndexFacePeriodic( void ) const;
    bool isIndexFacePeriodicSet() const ;
    
	bool isBorderNode(int nodeid) const ;
	bool isBorderFace(int faceid) const ;
	
	bool isTriangular() const ;
	bool isTetrahedral() const ;
	bool isQuadrangular() const ;
	bool isHexahedral() const ;
    bool isStructured() const ;

	// epsilon used in mesh comparisons
	double getComparisonEpsilon() const {return _epsilon;};
	void setComparisonEpsilon(double epsilon){ _epsilon=epsilon;};
    // Quick comparison of two meshes to see if they are identical with high probability (three cells are compared)
    void checkFastEquivalWith( Mesh m) const { return getMEDCouplingMesh()->checkFastEquivalWith(m.getMEDCouplingMesh(),_epsilon);};
    // Deep comparison of two meshes to see if they are identical Except for their names and units
    bool isEqualWithoutConsideringStr( Mesh m) const { return getMEDCouplingMesh()->isEqualWithoutConsideringStr(m.getMEDCouplingMesh(),_epsilon);};

    std::vector< std::string > getElementTypesNames() const ;
	/**
	 * \brief Compute the minimum value over all cells of the ratio cell perimeter/cell volume
	 */
    double minRatioVolSurf() const;
    
	/**
	 * \brief Compute the maximum number of neighbours around an element (cells around a cell or nodes around a node)
	 */
    int getMaxNbNeighbours(EntityType type) const;
    
    /** 
     * \brief Delete the medcoupling mesh to save memory space
     */
    void deleteMEDCouplingMesh();
    
    /** 
     * \brief Returns true iff an unstructured mesh has been loaded
     */
     bool meshNotDeleted() const {return _meshNotDeleted;}
    
private: //----------------------------------------------------------------

	MEDCoupling::MEDCouplingUMesh*  setMesh( void ) ;
	void setFaceGroups( const MEDCoupling::MEDFileUMesh* medmesh, MEDCoupling::MEDCouplingUMesh*  mu) ;//Read all face groups
	void setNodeGroups( const MEDCoupling::MEDFileMesh*  medmesh, MEDCoupling::MEDCouplingUMesh*  mu) ;//Read all node groups
	void addNewFaceGroup( const MEDCoupling::MEDCouplingUMesh *m);//adds one face group in the vectors _faceGroups, _faceGroupNames and _faceGroupIds
	
	/**
	 * \brief The MEDCoupling mesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> _mesh;// This is either a MEDCouplingUMesh or a MEDCouplingStructuredMesh

	bool _meshNotDeleted;
	
    std::string _name;
    
	/**
	 * \brief Space dimension
	 */
	int _spaceDim ;

	/**
	 * \brief Mesh dimension
	 */
	int _meshDim ;
    
    /**
     * \brief Signal a structured mesh
     */
	bool _isStructured;
    /**
     * \brief Number of cells in each direction (Structured meshes)
     */
	std::vector<mcIdType> _nxyz;

	/**
	 * \brief The nodes in this mesh.
	 */
	std::shared_ptr<Node> _nodes;

	/**
	 * \brief The number of nodes in this mesh.
	 */
	int _numberOfNodes;

	/**
	 * \brief The faces in this mesh.
	 */
	std::shared_ptr<Face> _faces;

	/**
	 * \brief The numbers of faces in this mesh.
	 */
	int _numberOfFaces;

	/**
	 * \brief The cells in this mesh.
	 */
	std::shared_ptr<Cell> _cells;

	/**
	 * \brief The number of cells in this mesh.
	 */
	int _numberOfCells;

	/**
	 * \brief return The nodes in this mesh
	 * @return _nodes
	 */
	std::shared_ptr<Node> getNodes ( void ) const ;

	/**
	 * \brief return The cells in this mesh
	 * @return _vertices
	 */
	std::shared_ptr<Cell> getCells ( void ) const ;

	/**
	 * \brief return the faces in this mesh
	 * @return _vertices
	 */
	std::shared_ptr<Face> getFaces ( void ) const ;

	/**
	 * \brief The number of edges in this mesh.
	 */
	int _numberOfEdges;//Useful to deduce the number of non zero coefficients in a finite element matrix 

	/**
	 * \brief The names of face groups.
	 */
	std::vector<std::string> _faceGroupNames;

	/**
	 * \brief The names of node groups.
	 */
	std::vector<std::string> _nodeGroupNames;

	/**
	 * \brief The list of face groups.
	 */
	std::vector<MEDCoupling::MEDCouplingUMesh *> _faceGroups;
	/**
	 * \brief The list of node groups.
	 */
	std::vector<MEDCoupling::DataArrayIdType *> _nodeGroups;
	
	/**
	 * \brief The list of face id in each face groups.
	 */
	std::vector< std::vector<int> > _faceGroupsIds;
	
	/**
	 * \brief The list of node id in each node groups.
	 */
	std::vector< std::vector<int> > _nodeGroupsIds;
	
	/**
	 * \brief Elements types (SEG2, TRI3, QUAD4, HEXA6 ...)
	 */
	std::vector< INTERP_KERNEL::NormalizedCellType > _eltsTypes;//List of cell types contained in the mesh
	std::vector< std::string > _eltsTypesNames;//List of cell types contained in the mesh
    std::vector< INTERP_KERNEL::NormalizedCellType > getElementTypes() const;    
    
    /**
     * \brief Tools to manage periodic boundary conditions in square/cube geometries
     */
     bool _indexFacePeriodicSet;
     std::map<int,int> _indexFacePeriodicMap;
    
    /* List of boundary faces*/
    std::vector< int > _boundaryFaceIds;
    /* List of boundary nodes*/
    std::vector< int > _boundaryNodeIds;
    /* Boundary mesh */
    MEDCoupling::MEDCouplingUMesh * _boundaryMesh;
    
    double _epsilon;
};

#endif /* MESH_HXX_ */
