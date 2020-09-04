/*
 * field.hxx
 *
 *  Created on: 07 fevrier. 2012
 *      Authors: CDMAT groups
 */

#ifndef FIELD_HXX_
#define FIELD_HXX_

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
  class DataArrayDouble;
}

#include "DoubleTab.hxx"
#include "Vector.hxx"
#include "Mesh.hxx"

#include <MCAuto.hxx>

/**
 * Field class is defined by
 * - ........
 */

class Field
{
    public: //----------------------------------------------------------------
    /**
     * default constructor
     */
    Field ( EntityType typeField = CELLS ) ;

    /**
    * constructor with data:
    * @param fieldName : name of the field
    * @param type : type of the field
    * @param mesh : mesh of the field
    * @param numberOfComponents : number of the component
    * @param time : time of the field
    */
    Field(const std::string fieldName, EntityType type, const Mesh& mesh, int numberOfComponents=1, double time=0.0) ;

    /**
    * destructor
    */
    ~Field ( void ) ;

    /**
    * constructor by copy
    * @param field : The Field object to be copied
    */
    Field ( const Field & field ) ;

    /**
     * deep copy of a field (values are copied)
     * @param field : The Field object to be copied
     */
    Field deepCopy( ) const;

    /**
     * constructor with data
     * @param filename : file name of field med file
     * @param fieldType: field type
     * @param fieldName: field name
     * @param iteration: iteration number (optional)
     * @param order:     order inside an iteration (optional)
     * @param numberOfComponents:     number of components of the field (optional)
     * @param time:     time index of the field (optional)
     */
    Field( const std::string filename, EntityType fieldType,
           const std::string & fieldName = "",
           int iteration = -1, int order = -1, int meshLevel=0,
           int numberOfComponents=1, double time=0.0);
  
    /**
     * constructor with data
	 * \brief defines a constant field on a mesh stored in a med file
	 * \details
	 * \param [in] string : the mesh file name
     * \param fieldType: field type
	 * \param [in] vector<double> : the value in each cell
     * \param [in] fieldName: field name
	 * \param [in] meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 *  */
	Field(const std::string meshfileName, EntityType fieldType, 
		  const std::vector<double> Vconstant,const std::string & fieldName = "",
		   int meshLevel=0, double time=0.0);

    /**
     * constructor with data
	 * \brief defines a constant field 
	 * \details
	 * \param [in] Mesh 
     * \param [in] fieldType: field type
	 * \param [in] Vector
     * \param [in] fieldName: field name
	 *  */
	Field(const Mesh& M, EntityType fieldType, const Vector Vconstant,
		  const std::string & fieldName = "", double time=0.0);

    /**
     * constructor with data
	 * \brief defines a constant field 
	 * \details
	 * \param [in] Mesh
     * \param [in] fieldType: field type
	 * \param [in] vector<double>
     * \param [in] fieldName: field name
	 *  */
	Field(const Mesh& M, EntityType fieldType, const std::vector<double> Vconstant, const std::string & fieldName = "", double time=0.0);

    /**
     * constructor with data
	 * \brief Builds a rectangular mesh M and defines a constant field on M
	 * \details
	 * \param [in] int the space dimension
	 * \param [in] vector<double> the value in each cell
     * \param [in] fieldType: field type
     * \param [in] fieldName: field name
	 * \param [in] double the lowest value in the x direction
	 * \param [in] double the highest value in the x direction
	 * \param [in] string name of the left boundary
	 * \param [in] string name of the right boundary
	 * \param [in] double the lowest value in the y direction
	 * \param [in] double the highest value in the y direction
	 * \param [in] string name of the back boundary
	 * \param [in] string name of the front boundary
	 * \param [in] double the lowest value in the z direction
	 * \param [in] double the highest value in the z direction
	 * \param [in] string name of the bottom boundary
	 * \param [in] string name of the top boundary
	 *  */
	Field( int nDim, const std::vector<double> Vconstant, EntityType type, 
			double xmin, double xmax,int nx, std::string leftSide, std::string rightSide,
			double ymin=0, double ymax=0, int ny=0, std::string backSide="", std::string frontSide="",
			double zmin=0, double zmax=0, int nz=0, std::string bottomSide="", std::string topSide="",
			const std::string & fieldName="", double time=0.0,double epsilon=1e-6);

    /**
     * constructor with data
	 * \brief Builds a step function field on the mesh M. The direction of the discontinuity is determined by the parameter "direction". The field takes value VV_left for x,y or z<disc_pos and VV_right for x,y or z>disc_pos
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector
	 * \param [in] Vector
	 * \param [in] double position of the discontinuity on one of the three axis
	 * \param [in] int direction (axis carrying the discontinuity) : 0 for x, 1 for y, 2 for z
     * \param [in] fieldType: field type
     * \param [in] fieldName: field name
	 *  */
	Field(const Mesh M, const Vector VV_left, const Vector VV_right, double disc_pos, 
			EntityType type, int direction=0, const std::string & fieldName="", double time=0.0);

    /**
     * constructor with data
	 * \brief Builds a rectangular mesh M and defines a step function field on M that takes values VV_left for x<xstep and VV_right for x>xstep
	 * \param [in] int the space dimension
	 * \param [in] vector<double> the value left of the discontinuity
	 * \param [in] vector<double> the value right of the discontinuity
	 * \param [in] double the position of the discontinuity in the x direction
     * \param [in] fieldType: field type
     * \param [in] fieldName: field name
	 * \param [in] double the lowest value in the x direction
	 * \param [in] double the highest value in the x direction
	 * \param [in] string name of the left boundary
	 * \param [in] string name of the right boundary
	 * \param [in] double the lowest value in the y direction
	 * \param [in] double the highest value in the y direction
	 * \param [in] string name of the back boundary
	 * \param [in] string name of the front boundary
	 * \param [in] double the lowest value in the z direction
	 * \param [in] double the highest value in the z direction
	 * \param [in] string name of the bottom boundary
	 * \param [in] string name of the top boundary
	 * \param [out] void
	 *  */
	Field( int nDim, const std::vector<double> VV_Left, std::vector<double> VV_Right, 
			double xstep, EntityType type,
			double xmin, double xmax,int nx, std::string leftSide, std::string rightSide,
			double ymin=0, double ymax=0, int ny=0, std::string backSide="", std::string frontSide="",
			double zmin=0, double zmax=0, int nz=0, std::string bottomSide="", std::string topSide="",
			int direction=0, const std::string & fieldName="", double time=0.0, double epsilon=1e-6);

    /**
     * constructor with data
	 * \brief builds a step function field on mesh M with values Vin inside the ball with radius Radius and Vout outside
	 * \details
	 * \param [in] Mesh
	 * \param [in] Vector Vin, value inside the ball
	 * \param [in] Vector Vout, value outside the ball
	 * \param [in] double radius of the ball
	 * \param [in] Vector Center, coordinates of the ball center
     * \param [in] fieldType: field type
     * \param [in] fieldName: field name
	 *  */
	Field(const Mesh M, const Vector Vin, const Vector Vout, double Radius, 
			Vector Center, EntityType type, const std::string & fieldName="", double time=0.0);

    void readFieldMed( const std::string & fileNameRadical,
                       EntityType type,
                       const std::string & fieldName = "",
                       int iteration = -1,
                       int order = -1) ;

    void buildFieldMemoryStructure();

    MEDCoupling::DataArrayDouble * getArray();

    double& operator[] ( int ielem ) ;

    double operator[] ( int ielem ) const;

    double& operator() ( int ielem ) ;

    double operator() ( int ielem ) const;

    double& operator() ( int ielem, int jcomp ) ;

    double operator() ( int ielem, int jcomp ) const ;

    int getNumberOfComponents ( void ) const ;

    const double* getValues ( void ) const ;

    const std::string getName ( void ) const;

    const Mesh& getMesh ( void ) const ;

    int getNumberOfElements ( void ) const ;

    EntityType getTypeOfField ( void ) const ;

    /**
     * return the MEDCouplingField pointer
     * return _field
     */
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> getField ( void )  const ;

    void setFieldByMEDCouplingFieldDouble ( const MEDCoupling::MEDCouplingFieldDouble* field );

    void setFieldByDataArrayDouble ( const MEDCoupling::DataArrayDouble* array );

    DoubleTab getNormEuclidean( void ) const ;

    double max( int component=0 ) const ;

    double min( int component=0 ) const ;

    void setTime ( double time, int iter );

    Vector getValuesOnComponent(int compo) const ;

    Vector getValuesOnAllComponents(int elem) const ;

    int getSpaceDimension( void ) const;

    double getTime ( void ) const;

    void setName ( const std::string fieldName ) ;

    void setInfoOnComponent(int icomp, std::string nameCompo) ;

    std::string getInfoOnComponent(int icomp) const;

    /**
     * Computes all the components of the sum of values of the field multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the integral of the field
     */
    Vector integral() const;

    /**
     * Computes the sum of values of a given component of the field multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * @param the index of the component of interest
     * return the numerical value of the integral of the field
     */
    double integral(int compId) const;

    /**
     * Computes for each component the sum of the absolute values of the field components multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh.
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the L1 norm of each component of the field
     */
    Vector normL1() const;

    /**
     * Computes all the components of the sum of squares of the values of the field components multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the L2 norm of each component of the field
     */
    Vector normL2() const;

    /**
     * Computes the maximum of each component of the field
     * The field may be multicomponent so the result of the function is a vector
     * return the vector of numerical value of the Linfinity norm of each component of the field
     */
    Vector normMax() const;

    /**
     * Computes the maximum of each component of the field as well as the index where the maximum was found
     * The field may be multicomponent so the result of the function is a vector of values and a vector of indices
     * return the vector of numerical value of the Linfinity norm of each component of the field AND the corresponding vector of indices
     */
    Vector componentMax(Vector & Indices) const;

    const Field& operator= ( const Field& f ) ;

    Field operator+ ( const Field& f ) const ;

    Field operator- ( const Field& f ) const ;

    const Field& operator+= ( const Field& f ) ;

    const Field& operator-= ( const Field& f ) ;

    const Field& operator*= ( double s ) ;

    const Field& operator/= ( double s ) ;

    const Field& operator-= ( double s ) ;

    const Field& operator+= ( double s ) ;

    void writeVTK ( const std::string fileName, bool fromScratch=true ) const ;

    void writeMED ( const std::string fileName, bool fromScratch=true ) const ;

    void writeCSV ( const std::string fileName ) const ;

    friend Field operator* (double value , const Field& field ) ;

    friend Field operator* (const Field& field, double value ) ;

    friend Field operator/ (const Field& field, double value) ;

    friend std::ostream& operator<<(std::ostream& out, const Field& field ) ;

    protected: //----------------------------------------------------------------

    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> _field;
    Mesh _mesh ;
    EntityType _typeField;
	int _numberOfComponents;
	double _time;
	std::string _fieldName;
	
    private:

};

#endif /* Field_HXX_ */
