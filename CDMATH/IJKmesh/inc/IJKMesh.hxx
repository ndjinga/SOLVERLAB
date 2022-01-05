/*
 * IJKmesh.hxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#ifndef IJKMESH_HXX_
#define IJKMESH_HXX_

/**
 * IJKMesh class is defined by
 * - case 1: file name of mesh med file (MEDCouplingCMesh)
 * - case 2: 1D cartesian, xmin and xmax and number of cells (MEDCouplingIMesh)
 * - case 3: 2D cartesian, xmin, xmax, ymin and ymax and numbers of cells in x direction and y direction (MEDCouplingIMesh)
 * - case 4: 3D cartesian, xmin, xmax, ymin, ymax, zmin and zmax and numbers of cells in x direction, y direction and z direction (MEDCouplingIMesh)
 * 
 *  Mesh cell and node structures are stored in a single MEDCouplingStructuredMesh _mesh
 *  nx    *  ny    *  nz    cell structure
 * (nx+1) * (ny+1) * (nz+1) node structure
 * 
 * The face structure is more tricky because it depends on the dimension
 * if dim=1, a single grid : 
 *                  nx+1 nodes
 * if dim=2, union of two grids : 
 *                  (nx+1)*ny nodes (faces orthogonal to x-axis), origin (0,dy/2) 
 *                  nx*(ny+1) nodes (faces orthogonal to x-axis), origin (dx/2,0)
 * if dim=3, union of three grids : 
 *                  nx*ny*(nz+1) (faces orthogonal to x-axis), origin (0,dy/2,dz/2) 
 *                  nx*(ny+1)*nz (faces orthogonal to y-axis), origin (dx/2,0,dz/2)  
 *                  (nx+1)*ny*nz (faces orthogonal to z-axis), origin (dx/2,dy/2,0)
 * 
 * Mesh face structures are stored in a vector of meshes _faceMeshes of size meshDim
 * The face centers are located on the nodes of the meshes in _faceMeshes
 * The origins of the meshes in _faceMeshes are shifted from the origin of _mesh by dx/2 on x-axis, dy/2 on y-axis and dz/2 on z-axis
 * 
 * - number  of nodes surounding each cell : known constant : 2*_Ndim
 * - number  of faces surounding each cell : known constant : 2 _Ndim
 * - normal vectors surrounding each cell
 * - measure of each cell : known constant : dx*dy*dz
 * - measures of faces : known constants : dx*dy, dx*dz and dy*dz
 */

namespace MEDCoupling
{
class MEDFileCMesh;
class MEDCouplingMesh;
class MEDCouplingIMesh;
class DataArrayIdType;
}
namespace ParaMEDMEM
{
class DataArrayIdType;
}
#include <MCAuto.hxx>
#include "NormalizedGeometricTypes"

class IJKNode;
class IJKCell;
class IJKFace;

enum EntityType
  {
    CELLS = 0,
    NODES = 1,
    FACES = 2,
  };

#include <vector>
#include <string>
#include <map>

class IJKMesh
{

public: //----------------------------------------------------------------
	/**
	 * \brief default constructor
	 */
	Mesh ( void ) ;

	/**
	 * \brief constructor with data to load a structured MEDCouplingIMesh
	 * @param filename name of structured mesh file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	Mesh ( const std::string filename, const std::string & meshName="" , int meshLevel=0 ) ;

	/**
	 * \brief constructor with data for a regular 1D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param nx : Number of cells in x direction
	 */
	Mesh( double xmin, double xmax, int nx, std::string meshName="MESH1D_Regular_Grid" ) ;

	/**
	 * \brief constructor with data for a regular 2D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName="MESH2D_Regular_Rectangle_Grid") ;

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
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName="MESH3D_Regular_Cuboid_Grid") ;

	Mesh( const MEDCoupling::MEDCouplingMesh* mesh ) ;

    /**
     * \brief Computes the global cell number from its IJK position
	 * @param int i : cell index along x-axis
	 * @param int j : cell index along y-axis
	 * @param int k : cell index along z-axis
	 * @return global cell number
     */
    int getCellNumber(int i, int j, int k) const { return _mesh->getCellIdFromPos ( i, j, k); };
    /**
     * \brief Computes the global node number from its IJK position
	 * @param int i : node index along x-axis
	 * @param int j : node index along y-axis
	 * @param int k : node index along z-axis
	 * @return global node number
     */
    int getNodeNumber(int i, int j, int k) const { return _mesh->getNodeIdFromPos ( i, j, k); };
    /**
     * \brief Computes the global face number from its IJK position
     * @param int face grid number corresponding to the direction of the face normal : 0->x, 1->y, 2->z
	 * @param int i : node index along x-axis
	 * @param int j : node index along y-axis
	 * @param int k : node index along z-axis
	 * @return global node number
     */
    int getFaceNumber(int face_grid_number, int i, int j, int k) const { return _faceMeshes[face_grid_number]->getNodeIdFromPos ( i, j, k); };
    /**
     * \brief Computes the IJK position of a cell from its index  number
	 * @param int cellId : global cell index 
	 * @return vector of i,j,k indices of the cell
     */
	std::vector< int > getIJKCellCoordinates(int cellId) const { return _mesh->getLocationFromCellId(cellId); };
    /**
     * \brief Computes the IJK position of a node from its index  number
	 * @param int nodeId : global node index 
	 * @return vector of i,j,k indices of the node
     */
	std::vector< int > getIJKNodeCoordinates(int nodeId) const { return _mesh->getLocationFromNodeId(nodeId); };
    /**
     * \brief Computes the IJK position of a face from its index  number
     * @param int face grid number corresponding to the direction of the face normal : 0->x, 1->y, 2->z
	 * @param int faceId : global face index 
	 * @return vector of i,j,k indices of the node
     */
	std::vector< int > getIJKFaceCoordinates( int face_grid_number, int faceId) const { return _faceMeshes[face_grid_number]->getLocationFromNodeId(faceId); };
    /**
     * \brief Computes the indices of nodes surrounding a given cell
	 * @param int cellId : global cell index 
	 * @return vector of indices of the nodes surrounding the cell 
     */
    std::vector< mcIdType > getNodeIdsOfCell(mcIdType cellId) const { std::vector< mcIdType > conn; _mesh_>getNodeIdsOfCell(mcIdType cellId, conn) ; return conn; };
    /**
     * \brief Computes the coordinates of a node
	 * @param int nodeId : global node index 
	 * @return vector of components of the node coordinates 
     */
    std::vector< double >   getNodeCoordinates (mcIdType nodeId) const { std::vector< double > coo; _mesh_>getCoordinatesOfNode(mcIdType nodeId, coo) ; return coo; };
    /**
     * \brief Computes the coordinates of a face center of mass
     * @param int face grid number corresponding to the direction of the face normal : 0->x, 1->y, 2->z
	 * @param int faceId : global face index 
	 * @return vector of components of the face coordinates 
     */
    std::vector< double >   getFaceCoordinates ( int face_grid_number, mcIdType faceId) const { std::vector< double > coo; _faceMeshes[face_grid_number]>getCoordinatesOfNode(mcIdType faceId, coo) ; return coo; };
    /**
     * \brief Computes the isobarycenter of a cell
	 * @param int cellId : cell number
     */
    std::vector< double >   getCellCenterCoordinates (mcIdType cellId) const ;
    
	 double getMeasureOfAnyCell () const;
	 
	/**
	 * \brief constructor with data
	 * @param filename : file name of structured mesh med file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	void readMeshMed( const std::string filename, const std::string & meshName="" , int meshLevel=0 ) ;

	/**
	 * \brief constructor by copy
	 * @param mesh : The Mesh object to be copied
	 */
	Mesh ( const IJKMesh & mesh ) ;

	/**
	 * \brief destructor
	 */
	~Mesh( void ) ;

	/**
	 * \brief return mesh name
	 * @return _name
	 */
	std::string getName( void ) const { return _mesh->getName (); };

	/**
	 * \brief return Space dimension
	 * @return spaceDim
	 */
	int getSpaceDimension( void ) const { return _mesh->getSpaceDimension (); };

	/**
	 * \brief return Mesh dimension
	 * @return meshDim
	 */
	int getMeshDimension( void ) const { return _mesh->getMeshDimension (); };

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

	std::vector<int> getCellGridStructure() const;

	std::vector<int> getNodeGridStructure() const;

	std::vector< vector< int > > getFaceGridStructures() const;

	/**
	 * \brief overload operator =
	 * @param mesh : The Mesh object to be copied
	 */
	const Mesh& operator= ( const Mesh& mesh ) ;

	/**
	 * \brief return the MEDCouplingStructuredMesh
	 * @return _mesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> getMEDCouplingStructuredMesh ( void )  const ;

	/**
	 * \brief return the MEDCouplingStructuredMesh
	 * @return _faceMesh
	 */
	std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> > getFaceMeshes const;	

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
	 * \brief write the cell mesh in the VTK format
	 */
	void writeVTK ( const std::string fileName ) const ;

	/**
	 * \brief write the cell and face meshes in the VTK format
	 */
	void writeVTKAll meshes ( const std::string fileName ) const ;

	/**
	 * \brief write the cell mesh in the MED format
	 */
	void writeMED ( const std::string fileName, bool fromScratch = true ) const ;

	/**
	 * \brief write the cell and face meshes in the MED format
	 */
	void writeMEDAllMeshes ( const std::string fileName, bool fromScratch = true ) const ;

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
	
	bool isQuadrangular() const ;
	bool isHexahedral() const ;
    bool isStructured() const ;

	// epsilon used in mesh comparisons
	double getComparisonEpsilon() const {return _epsilon;};
	void setComparisonEpsilon(double epsilon){ _epsilon=epsilon;};
    // Quick comparison of two meshes to see if they are identical with high probability (three cells are compared)
    void checkFastEquivalWith( Mesh m) const { return _mesh()->checkFastEquivalWith(m.getMEDCouplingStructuredMesh(),_epsilon);};
    // Deep comparison of two meshes to see if they are identical Except for their names and units
    bool isEqualWithoutConsideringStr( Mesh m) const { return _mesh->isEqualWithoutConsideringStr(m.getMEDCouplingStructuredMesh(),_epsilon);};

    std::vector< std::string > getElementTypesNames() const ;
	/**
	 * \brief Compute the minimum value over all cells of the ratio cell perimeter/cell volume
	 */
    double minRatioVolSurf() const{ return _cellMeasure / *max_element(begin(_faceMeasures), end(_faceMeasures));};
    
	/**
	 * \brief Compute the maximum number of neighbours around an element (cells around a cell or nodes around a node)
	 */
    int getMaxNbNeighbours(EntityType type) const;
    
    /**
	 * return the measure of a cell (length in 1D, surface in 2D or volume in 3D)
	 * @return _cellMeasure
	 */
	double getCellMeasure ( ) const { return _cellMeasure;};

    
    /**
	 * return the measure of a cell (length in 1D, surface in 2D or volume in 3D)
	 * @return _faceMeasures
	 */
	std::vector< double > getFaceMeasures ( ) const { return _faceMeasures;};
	/**
	 * return normal vectors around each cell
	 */
	std::vector< Vector > getNormalVectors (void) const { return _faceNormals;};


private: //----------------------------------------------------------------

    /**
     * \brief The cell mesh
	 * Question : can _mesh be const since no buildUnstructured is applied?
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> _mesh;//This is either a MEDCouplingIMesh (creation from scratch) or a MEDCouplingCMesh (loaded from a file)
    /**
     * \brief The face meshes
	 */
	std::vector< MEDCoupling::MCAuto<MEDCoupling::MEDCouplingStructuredMesh> > _faceMeshes;//These are MEDCouplingIMesh
    /**
     * \brief Generate the face meshes from the cell mesh
	 */
	void setFaceMeshes();

	double _cellMeasure;
	std::vector< double > _faceMeasures;
	std::vector< std::vector< double > > _faceNormals;

	/* Boundary data */
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
     *  \brief List of boundary faces
     */
    std::vector< int > _boundaryFaceIds;
    /**
     *  \brief List of boundary nodes
     */
    std::vector< int > _boundaryNodeIds;
    
    /**
     * \brief Tools to manage periodic boundary conditions in square/cube geometries
     */
     bool _indexFacePeriodicSet;
     std::map<int,int> _indexFacePeriodicMap;
    
    double _epsilon;
};

#endif /* IJKMESH_HXX_ */
