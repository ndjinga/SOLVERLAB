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
 * - case 1: file name of mesh med file (MEDCouplingIMesh)
 * - case 2: 1D cartesian, xmin and xmax and number of cells
 * - case 3: 2D cartesian, xmin, xmax, ymin and ymax and numbers of cells in x direction and y direction
 * - case 4: 3D cartesian, xmin, xmax, ymin, ymax, zmin and zmax and numbers of cells in x direction, y direction and z direction
 * 
 *  Mesh structure
 * nx  , ny  , nz   cell structure
 * nx+1, ny+1, nz+1 node structure
 * 
 * The face structure is more tricky because it depends on the dimension
 * if dim=1, a single grid : nx+1 
 * if dim=2, union of two grids : nx,ny+1  and nx,ny+1
 * if dim=3, union of three grids : nx,ny,nz+1, nx,ny+1,nz  and nx+1,ny,nz
 * 
 * - number  of nodes surounding each cell : known constant : 2*_Ndim
 * - number  of faces surounding each cell : known constant : 2 _Ndim
 * - normal vectors surrounding each cell
 * - measure of each cell : known constant : dx*dy*dz
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

typedef enum
  {
    CELLS = 0,
    NODES = 1,
    FACES = 2,
  } EntityType;

#include <vector>
#include <string>
#include <map>

class IJKMesh
{

public: //----------------------------------------------------------------
	/**
	 * default constructor
	 */
	Mesh ( void ) ;

	/**
	 * constructor with data to load a structured MEDCouplingIMesh
	 * @param filename name of structured mesh file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	Mesh ( const std::string filename, int meshLevel=0 ) ;

	/**
	 * constructor with data for a regular 1D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param nx : Number of cells in x direction
	 */
	Mesh( double xmin, double xmax, int nx, std::string meshName="MESH1D_Regular_Grid" ) ;

	/**
	 * constructor with data for a regular 2D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName="MESH2D_Regular_Rectangle_Grid") ;

	/**
	 * constructor with data for a regular 3D grid 
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

	Mesh( const MEDCoupling::MEDCouplingIMesh* mesh ) ;

	/**
	 * constructor with data
	 * @param filename : file name of structured mesh med file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	void readMeshMed( const std::string filename, int meshLevel=0 ) ;

	/**
	 * constructor by copy
	 * @param mesh : The Mesh object to be copied
	 */
	Mesh ( const IJKMesh & mesh ) ;

	/**
	 * destructor
	 */
	~Mesh( void ) ;

	/**
	 * return mesh name
	 * @return _name
	 */
	std::string getName( void ) const ;

	/**
	 * return Space dimension
	 * @return _spaceDim
	 */
	int getSpaceDimension( void ) const ;

	/**
	 * return Mesh dimension
	 * @return _meshDim
	 */
	int getMeshDimension( void ) const ;

	/**
	 * return the number of nodes in this mesh
	 * @return _numberOfNodes
	 */
	int getNumberOfNodes ( void )  const ;

	/**
	 * return the number of faces in this mesh
	 * @return _numberOfFaces
	 */
	int getNumberOfFaces ( void )  const ;

	/**
	 * return the number of cells in this mesh
	 * @return _numberOfCells
	 */
	int getNumberOfCells ( void )  const ;

	/**
	 * return The cell i in this mesh
	 * @return _cells[i]
	 */
	IJKCell& getCell ( int i )  const ;

	/**
	 * return The face i in this mesh
	 * @return _faces[i]
	 */
	IJKFace& getFace ( int i )  const ;

	/**
	 * return The node i in this mesh
	 * @return _nodes[i]
	 */
	IJKNode& getNode ( int i )  const ;

	/**
	 * return number of cell in x direction (structured mesh)
	 * return _nX
	 */
	int getNx( void )  const ;

	/**
	 * return number of cell in y direction (structured mesh)
	 * return _nY
	 */
	int getNy( void )  const ;

	/**
	 * return number of cell in z direction (structured mesh)
	 * return _nZ
	 */
	int getNz( void )  const ;

	double getXMin( void )  const ;

	double getXMax( void )  const ;

	double getYMin( void )  const ;

	double getYMax( void )  const ;

	double getZMin( void )  const ;

	double getZMax( void )  const ;

	std::vector<double> getDXYZ() const ;

	std::vector<int> getCellGridStructure() const;

	std::vector<int> getNodeGridStructure() const;

	/**
	 * surcharge operator =
	 * @param mesh : The Mesh object to be copied
	 */
	const Mesh& operator= ( const Mesh& mesh ) ;

	/**
	 * return the mesh MEDCoupling
	 * return _mesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingIMesh> getMEDCouplingIMesh ( void )  const ;

	/**
	 * return the list of face group names
	 * return _faaceGroupNames
	 */
	std::vector<std::string> getNameOfFaceGroups( void )  const ;

	/**
	 * return the list of node group names
	 * return _nodeGroupNames
	 */
	std::vector<std::string> getNameOfNodeGroups( void )  const ;

	/**
	 * write mesh in the VTK format
	 */
	void writeVTK ( const std::string fileName ) const ;

	/**
	 * write mesh in the MED format
	 */
	void writeMED ( const std::string fileName ) const ;

	/*
     * Functions to manage periodic boundary condition in square/cubic geometries 
     */
    void setPeriodicFaces() ;
    int getIndexFacePeriodic(int indexFace, bool check_groups= false, bool use_central_inversion=false);
    void setBoundaryNodes();
    bool isIndexFacePeriodicSet() const ;
    
	bool isBorderNode(int nodeid) const ;
	bool isBorderFace(int faceid) const ;
	
	bool isQuadrangular() const ;
	bool isHexahedral() const ;
    bool isStructured() const ;
    std::string getElementTypes() const ;
    
	/**
	 * Compute the minimum value over all cells of the ratio cell perimeter/cell vaolume
	 */
    double minRatioVolSurf();
    
	/**
	 * Return the maximum number of neighbours around an element (cells around a cell or nodes around a node)
	 */
    int getMaxNbNeighbours(EntityType type) const;
    
    /**
	 * return the measure of all cells (length in 1D, surface in 2D or volume in 3D)
	 * @return _measureOfCells
	 */
	double getCellsMeasure ( void ) const ;

	/**
	 * return normal vectors around each cell
	 */
	std::vector< Vector > getNormalVectors (void) const ;


private: //----------------------------------------------------------------

    std::string _name;
    
	/**
	 * Space dimension
	 */
	int _spaceDim ;

	/**
	 * Mesh dimension
	 */
	int _meshDim ;
    
    /*
     * Structured mesh parameters
     */

	double _xMin;

	double _xMax;

	double _yMin;

	double _yMax;

	double _zMin;

	double _zMax;

	std::vector<int> _nxyz;//Number of cells in each direction

	std::vector<double> _dxyz;//lenght depth and height of each cell

	/*
	 * The number of nodes in this mesh.
	 */
	int _numberOfNodes;

	/*
	 * The numbers of faces in this mesh.
	 */
	int _numberOfFaces;

	/*
	 * The number of cells in this mesh.
	 */
	int _numberOfCells;

	/*
	 * The names of face groups.
	 */
	std::vector<std::string> _faceGroupNames;

	/*
	 * The names of node groups.
	 */
	std::vector<std::string> _nodeGroupNames;

	/*
	 * The list of face groups.
	 */
	std::vector<MEDCoupling::MEDCouplingUMesh *> _faceGroups;
	/*
	 * The list of node groups.
	 */
	std::vector<MEDCoupling::DataArrayIdType *> _nodeGroups;
	/*
	 * The mesh MEDCouplingIMesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingIMesh> _mesh;
	std::vector< INTERP_KERNEL::NormalizedCellType > _eltsTypes;//List of cell types contained in the mesh
    
    /*
     * Tools to manage periodic boundary conditions in square/cube geometries
     */
     bool _indexFacePeriodicSet;
     std::map<int,int> _indexFacePeriodicMap;
    
    /* List of boundary faces*/
    std::vector< int > _boundaryFaceIds;
    /* List of boundary nodes*/
    std::vector< int > _boundaryNodeIds;
    
	/*
	 * The number of nodes surrounding each cell.
	 */
   	int _numberOfNodesSurroundingCells ;

	/*
	 * The number of faces surounding each cell.
	 */
	int _numberOfFacesSurroundingCells ;

	/*
	 * The length or surface or volume of each cell.
	 */
	double _measureOfCells ;

	/*
	 * The normal vectors surrounding each cell.
	 */
	std::vector< Vector > _normalVectorsAroundCells ;
};

#endif /* IJKMESH_HXX_ */
