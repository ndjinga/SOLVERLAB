/*
 * field.cxx
 *
 *  Created on: 7 fevrier. 2012
 *      Author: CDMAT
 */

#include "Node.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Field.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField1TS.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "CdmathException.hxx"

#include <fstream>
#include <sstream>
#include <cstring>

using namespace MEDCoupling;
using namespace std;


//----------------------------------------------------------------------
Field::Field( EntityType typeField )
//----------------------------------------------------------------------
{
	_field=NULL;
	_typeField=typeField;
	_numberOfComponents=0;
}

//----------------------------------------------------------------------
Field::~Field( void )
//----------------------------------------------------------------------
{
	//std::cerr << "dtor Field, _field = " <<_field << std::endl;
	//if (_field) _field->decrRef();
}


Field::Field(const std::string fieldName, EntityType type, const Mesh& mesh, int numberOfComponents, double time)
{
	_field = NULL;
	_mesh=Mesh(mesh);
	_typeField=type;
	_numberOfComponents=numberOfComponents;
	_time=time;
	_fieldName=fieldName;

	buildFieldMemoryStructure();
}

Field::Field(const MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> field )
{
	_field = field->clone(true);//Use same mesh but new data array for the new field
	_numberOfComponents = _field->getNumberOfComponents() ;
	int iteration, order;
	_time = _field->getTime(iteration, order);

	_mesh=Mesh(_field->getMesh() );
	switch (_field->getTypeOfField()) 
	{
		case ON_CELLS:
			_typeField=CELLS;
			break;
		case ON_NODES:
			_typeField=NODES;
			break;
		default:
			throw CdmathException("Input field type cannot be used in SOLVERLAB. Field values should be on nodes or cells.");
	}
	
	_fieldName=_field->getName();
}

//Allocation d'un tableau pour le stockage des valeurs du champs_field
void Field::buildFieldMemoryStructure()
{
	MEDCouplingUMesh* mu=_mesh.getMEDCouplingMesh()->buildUnstructured();
	DataArrayDouble *array=DataArrayDouble::New();
	if (_typeField==CELLS)
	{
		_field=MEDCouplingFieldDouble::New(ON_CELLS);
		array->alloc(_mesh.getNumberOfCells(),_numberOfComponents);
		_field->setMesh(mu);
	}else if(_typeField==NODES)
	{
		_field=MEDCouplingFieldDouble::New(ON_NODES);
		array->alloc(mu->getNumberOfNodes(),_numberOfComponents);
		_field->setMesh(mu);
	}else if(_typeField==FACES)
	{
		_field=MEDCouplingFieldDouble::New(ON_CELLS);
		array->alloc(_mesh.getNumberOfFaces(),_numberOfComponents);
		DataArrayIdType *desc=DataArrayIdType::New();
		DataArrayIdType *descI=DataArrayIdType::New();
		DataArrayIdType *revDesc=DataArrayIdType::New();
		DataArrayIdType *revDescI=DataArrayIdType::New();
		MEDCouplingUMesh *m3=mu->buildDescendingConnectivity(desc,descI,revDesc,revDescI);
		_field->setMesh(m3);
		desc->decrRef();
		descI->decrRef();
		revDesc->decrRef();
		revDescI->decrRef();
		m3->decrRef();
	}else if(_typeField==GAUSS_PT)
	{
        int spaceDimension = _mesh.getSpaceDimension();
        int meshDimension  = _mesh.getMeshDimension();
        int nbGaussPoints, nbCellNodes;
        std::vector<double> refCoords, gaussCoords, weights;
        INTERP_KERNEL::NormalizedCellType cellType;
        
        if( meshDimension!=1 && !_mesh.isTriangular() && !_mesh.isTetrahedral())
            throw CdmathException("Field construction failed : Support mesh should be composed of segments (1D) triangles (2D) or tetrahedra to allow for the use of Gauss point");

        if(meshDimension==1)//GaussLegendre3 (qf3pE) from FreeFem: 3 quadrature points on each segment (order 3)
        {
            nbGaussPoints=3;
            nbCellNodes=2;
            cellType=INTERP_KERNEL::NORM_SEG2;

            const double gauss_n3_0=  0.5 ;
            const double gauss_n3_1=  (1-sqrt(3./5.)) /2  ;
            const double gauss_n3_2 =  1 - gauss_n3_1 ;
            
            const double pgauss_n3_0=  8./18.;
            const double pgauss_n3_1=  5./18.;
            const double pgauss_n3_2=  5./18.;
                    
            if( meshDimension == spaceDimension )//1D line in 1D space
            {
                refCoords   = std::vector<double>{0.0, 1.0};//nbCellNodes*spaceDimension
                gaussCoords = std::vector<double>{gauss_n3_0, gauss_n3_1, gauss_n3_2};//nbGaussPoints*spaceDimension
            }
            else if( 1+meshDimension == spaceDimension )//1D line embedded in 2D space
            {
                refCoords   = std::vector<double>{ 0.0, 0.0, 1.0 , 0.0 };//nbCellNodes*spaceDimension
                gaussCoords = std::vector<double>{gauss_n3_0, 0, gauss_n3_1, 0, gauss_n3_2, 0};//nbGaussPoints*spaceDimension
            }
            else if( 2+meshDimension == spaceDimension )//1D line embedded in 3D space
            {
                refCoords   = std::vector<double>{ 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0 };//nbCellNodes*spaceDimension
                gaussCoords = std::vector<double>{gauss_n3_0 , 0.0, 0.0, gauss_n3_1 , 0.0, 0.0, gauss_n3_2 , 0.0, 0.0};//nbGaussPoints*spaceDimension
            }
            weights=std::vector<double>{ pgauss_n3_0, pgauss_n3_1, pgauss_n3_2};//nbGaussPoints
        }
        else if(meshDimension==2)//QuadratureFormular_T_5 (qf5pT) from FreeFem: 7 quadrature points on each triangle (order 5)
        {
            nbGaussPoints=7;
            nbCellNodes  =3;
            cellType=INTERP_KERNEL::NORM_TRI3;

            // ----------------------------------------------------------------------
            // STROUD page  314 
            // -----------------------------
            const double sqrt15 = 3.87298334620741688517926539978;
            const double t_T5 =1.E0/3.E0        ,                           A_T5 = 0.225E0;
            const double r_T5 = (6-sqrt15)/21   ,  s_T5 = (9+2*sqrt15)/21 , B_T5 = (155-sqrt15)/1200;
            const double u_T5 = (6+sqrt15)/21   ,  v_T5 = (9-2*sqrt15)/21 , C_T5 = (155+sqrt15)/1200;
            
            if( meshDimension == spaceDimension )//2D surface
            {
                refCoords   = std::vector<double>{0.0, 0.0, 1.0 , 0.0, 0.0, 1.0};//nbCellNodes*spaceDimension
                gaussCoords = std::vector<double>{t_T5, t_T5, r_T5, r_T5, r_T5, s_T5, s_T5, r_T5, u_T5, u_T5, u_T5, v_T5, v_T5, u_T5};//nbGaussPoints*spaceDimension
            }
            else//2D surface embedded in 3D
            {
                refCoords = std::vector<double>{ 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0, 0.0, 1.0, 0.0 };//nbCellNodes*spaceDimension
                gaussCoords = std::vector<double>{t_T5, t_T5, 0, r_T5, r_T5, 0, r_T5, s_T5, 0, s_T5, r_T5, 0, u_T5, u_T5, 0, u_T5, v_T5, 0, v_T5, u_T5, 0};//nbGaussPoints*spaceDimension
            }
            weights=std::vector<double>{A_T5, B_T5, B_T5, B_T5, C_T5, C_T5, C_T5};//nbGaussPoints
        }
        else //QuadratureFormular_Tet_5 (qfV5) from FreeFem: 14 quadrature points on each tetrahedron (order 5)
        {
            nbGaussPoints=14;
            nbCellNodes  =4;
            cellType=INTERP_KERNEL::NORM_TETRA4;

            // 5  14  (formule 1)
            /*
            GM78
            A. Grundmann and H.M. Möller, Invariant integration formulas for the n-simplex by combinatorial methods, SIAM J. Numer. Anal. 15 (1978), 282--290.
             */
            const double w1 = 0.0122488405193936582572850342477212*6;
            const double w5 = 0.0187813209530026417998642753888810*6;
            const double w9 = 7.09100346284691107301157135337624E-3*6;
            const double c1 = 0.7217942490673263207930282587889082;
            const double c2 = 0.0927352503108912264023239137370306;
            const double c3 = 0.067342242210098170607962798709629;
            const double c4 = 0.310885919263300609797345733763457;
            const double c5 = 0.454496295874350350508119473720660;
            const double c6 = 00.045503704125649649491880526279339;
            
            refCoords   = std::vector<double>{ 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };//nbCellNodes*spaceDimension
            gaussCoords = std::vector<double>{ c1, c2, c2, 
                                               c2, c1, c2, 
                                               c2, c2, c1, 
                                               c2, c2, c2, 
                                               c3, c4, c4, 
                                               c4, c3, c4,
                                               c4, c4, c3,
                                               c4, c4, c4,
                                               c5, c5, c6,
                                               c5, c6, c5,
                                               c6, c5, c5,
                                               c6, c6, c5,
                                               c6, c5, c6,
                                               c5, c6, c6 };//nbGaussPoints*spaceDimension
            weights=std::vector<double>{ w1, w1, w1, w1, w5, w5, w5, w5, w9, w9, w9, w9, w9, w9 };//nbGaussPoints
        }
		_field=MEDCouplingFieldDouble::New(ON_GAUSS_PT);
		array->alloc(_mesh.getNumberOfCells(),_numberOfComponents*nbGaussPoints);
		_field->setMesh(mu);
        _field->setGaussLocalizationOnType(cellType,refCoords,gaussCoords,weights);
	}else
		throw CdmathException("Type of Mesh cells is not compatible. Cell types accepted are 1D segments (NORM_SEG2), 2D triangles (NORM_TRI3) and 3D tetrahedra (NORM_TETRA4)");

	_field->setName(_fieldName.c_str()) ;
	_field->setArray(array);
	_field->setTime(_time,0,0);
	array->decrRef();
}

int Field::getNumberOfGaussPtPerCell()
{
    if( _typeField != GAUSS_PT )
        throw CdmathException("Field::getNumberOfGaussPtPerCell() : Field should be of type GAUSS_PT (not CELLS, not NODES not FACES)");
        
    INTERP_KERNEL::NormalizedCellType cellType;

    switch( _mesh.getMeshDimension() )
    {
        case 1: 
            cellType = INTERP_KERNEL::NORM_SEG2;
            break;
        case 2: 
            cellType = INTERP_KERNEL::NORM_TRI3;
            break;
        case 3: 
            cellType = INTERP_KERNEL::NORM_TETRA4;
            break;
        default: 
    		throw CdmathException("Type of Mesh cells is not compatible. Cell types accepted are 1D segments (NORM_SEG2), 2D triangles (NORM_TRI3) and 3D tetrahedra (NORM_TETRA4)");
    }
    
    int locID = _field->getGaussLocalizationIdOfOneType ( cellType );//returns exception if different nb of gauss points in cells

    return _field->getGaussLocalization( locID ).getNumberOfGaussPt( );
}

Field::Field( const std::string filename, EntityType type,
		const std::string & fieldName,
		int iteration, int order, int meshLevel)
{
	_field = NULL;
	_typeField=type;
	_fieldName=fieldName;

	readFieldMed(filename, type, fieldName, iteration, order);
}

Field::Field(const std::string meshFileName, EntityType type, const std::vector<double> Vconstant, 
		const std::string & fieldName, int meshLevel, double time, std::string meshName )
{
	_field = NULL;
	_mesh=Mesh(meshFileName + ".med", meshName, meshLevel);
	_typeField=type;
	_numberOfComponents=Vconstant.size();
	_time=time;
	_fieldName=fieldName;

	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	for (int ielem=0 ; ielem<nbelem; ielem++)
		for (int jcomp=0 ; jcomp<nbcomp ; jcomp++)
			_field->getArray()->getPointer()[jcomp+ielem*nbcomp]=Vconstant[jcomp];
}
Field::Field(const Mesh& M, EntityType type, const Vector Vconstant, const std::string & fieldName, double time)
{
	_field = NULL;
	_mesh=Mesh(M);
	_typeField=type;
	_numberOfComponents=Vconstant.size();
	_time=time;
	_fieldName=fieldName;

	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	for (int ielem=0 ; ielem<nbelem; ielem++)
		for (int jcomp=0 ; jcomp<nbcomp ; jcomp++)
			_field->getArray()->getPointer()[jcomp+ielem*nbcomp]=Vconstant[jcomp];
}
Field::Field(const Mesh& M, EntityType type, const vector<double> Vconstant, const std::string & fieldName, double time) 
{
	_field = NULL;
	_mesh=Mesh(M);
	_typeField=type;
	_numberOfComponents=Vconstant.size();
	_time=time;
	_fieldName=fieldName;

	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	for (int ielem=0 ; ielem<nbelem; ielem++)
		for (int jcomp=0 ; jcomp<nbcomp ; jcomp++)
			_field->getArray()->getPointer()[jcomp+ielem*nbcomp]=Vconstant[jcomp];
}
Field::Field( int nDim, const vector<double> Vconstant, EntityType type, 
		double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide, 
		const std::string & fieldName, double time,double epsilon)
{
	_field = NULL;
	_typeField=type;
	_numberOfComponents=Vconstant.size();
	_time=time;
	_fieldName=fieldName;

	//Build mesh
	if(nDim==1){
		_mesh=Mesh(xmin,xmax,nx);
	}
	else if(nDim==2)
		_mesh=Mesh(xmin,xmax,nx,ymin,ymax,ny);
	else if(nDim==3)
		_mesh=Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
	else{
		cout<<"Field::Field: Space dimension nDim should be between 1 and 3"<<endl;
		throw CdmathException("Space dimension nDim should be between 1 and 3");
	}

	_mesh.setGroupAtPlan(xmax,0,epsilon,rightSide);
	_mesh.setGroupAtPlan(xmin,0,epsilon,leftSide);
	if(nDim>=2){
		_mesh.setGroupAtPlan(ymax,1,epsilon,frontSide);
		_mesh.setGroupAtPlan(ymin,1,epsilon,backSide);
	}
	if(nDim==3){
		_mesh.setGroupAtPlan(zmax,2,epsilon,topSide);
		_mesh.setGroupAtPlan(zmin,2,epsilon,bottomSide);
	}

	// Build field
	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	for (int ielem=0 ; ielem<nbelem; ielem++)
		for (int jcomp=0 ; jcomp<nbcomp ; jcomp++)
			_field->getArray()->getPointer()[jcomp+ielem*nbcomp]=Vconstant[jcomp];
}
Field::Field(const Mesh M, const Vector VV_Left, const Vector VV_Right, double disc_pos,
		EntityType type, int direction, const std::string & fieldName, double time)
{
	if  (VV_Right.getNumberOfRows()!=VV_Left.getNumberOfRows())
		throw CdmathException( "Field::Field: Vectors VV_Left and VV_Right have different sizes");

	_field = NULL;
	_mesh=Mesh(M);
	_typeField=type;
	_numberOfComponents=VV_Left.getNumberOfRows();
	_time=time;
	_fieldName=fieldName;

	// Build field
	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;
	double component_value;

	for (int j = 0; j < nbelem; j++) {
		if(direction==0)
			component_value=M.getCell(j).x();
		else if(direction==1)
			component_value=M.getCell(j).y();
		else if(direction==2)
			component_value=M.getCell(j).z();
		else
			throw CdmathException( "Field::Field: direction should be an integer between 0 and 2");

		for (int i=0; i< nbcomp; i++)
			if (component_value< disc_pos )
				_field->getArray()->getPointer()[j+i*nbcomp] = VV_Left[i];
			else
				_field->getArray()->getPointer()[j+i*nbcomp] = VV_Right[i];
	}
}
Field::Field( int nDim, const vector<double> VV_Left, vector<double> VV_Right, 
		double xstep, EntityType type,
		double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide,
		int direction, const std::string & fieldName, double time, double epsilon)
{
	if  (VV_Right.size()!=VV_Left.size())
		throw CdmathException( "Field::Field: Vectors VV_Left and VV_Right have different sizes");
	Mesh M;
	if(nDim==1)
		M=Mesh(xmin,xmax,nx);
	else if(nDim==2)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny);
	else if(nDim==3)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
	else
		throw CdmathException("Field::Field : Space dimension nDim should be between 1 and 3");

	M.setGroupAtPlan(xmax,0,epsilon,rightSide);
	M.setGroupAtPlan(xmin,0,epsilon,leftSide);
	if(nDim>=2){
		M.setGroupAtPlan(ymax,1,epsilon,frontSide);
		M.setGroupAtPlan(ymin,1,epsilon,backSide);
	}
	if(nDim==3){
		M.setGroupAtPlan(zmax,2,epsilon,topSide);
		M.setGroupAtPlan(zmin,2,epsilon,bottomSide);
	}
	Vector V_Left(VV_Left.size()), V_Right(VV_Right.size());
	for(int i=0;i<VV_Left.size(); i++){
		V_Left(i)=VV_Left[i];
		V_Right(i)=VV_Right[i];
	}
	Field(M, V_Left, V_Right, xstep,  type, direction, fieldName, time);
}

Field::Field(const Mesh M, const Vector Vin, const Vector Vout, double radius, 
		const Vector Center, EntityType type, const std::string & fieldName, double time)
{
	if((Center.size()!=M.getSpaceDimension()) || (Vout.size() != Vin.size()) )
	{
		cout<< "Vout.size()= "<<Vout.size() << ", Vin.size()= "<<Vin.size()<<", Center.size()="<<Center.size()<<", M.getSpaceDim= "<< M.getSpaceDimension()<<endl;
		throw CdmathException("Field::Field : Vector size error");
	}

	_field = NULL;
	_mesh=Mesh(M);
	_typeField=type;
	_numberOfComponents=Vout.size();
	_time=time;
	_fieldName=fieldName;

	// Build field
	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	int spaceDim=M.getSpaceDimension();
	Vector currentPoint(spaceDim);

	for(int i=0;i<nbelem;i++)
	{
		currentPoint(0)=M.getCell(i).x();
		if(spaceDim>1)
		{
			currentPoint(1)=M.getCell(i).y();
			if(spaceDim>2)
				currentPoint(2)=M.getCell(i).z();
		}
		if((currentPoint-Center).norm()<radius)
			for(int j=0;j<nbcomp;j++)
				_field->getArray()->getPointer()[j+i*nbcomp]=Vin[j];
		else
			for(int j=0;j<nbcomp;j++)
				_field->getArray()->getPointer()[j+i*nbcomp]=Vout[j];
	}
}

MEDCoupling::DataArrayDouble * Field::getArray(){
	return _field->getArray();
}

void
Field::readFieldMed( const std::string & fileNameRadical,
		EntityType type,
		const std::string & fieldName,
		int iteration,
		int order)
{
	/**
	 * Reads the file fileNameRadical.med and creates a Field from it.
	 */
	std::string completeFileName = fileNameRadical + ".med";
	std::vector<std::string> fieldNames = MEDCoupling::GetAllFieldNames(completeFileName);
	size_t iField = 0;
	std::string attributedFieldName;
	_field = NULL;

	// Get the name of the right field that we will attribute to the Field.
	if (fieldName == "") {
		if (fieldNames.size() > 0)
		{
			cout<<"Warning : No field name imposed, taking the first field name found : "<< fieldNames[0]<<" in file "<<completeFileName<<endl;;
			attributedFieldName = fieldNames[0];
		}
		else 
		{
			std::ostringstream message;
			message << "Warning : No field found in file " << completeFileName;
			throw CdmathException(message.str().c_str());
		}
	}
	else 
	{
		int mycount = std::count(fieldNames.begin(), fieldNames.end(), fieldName);
		
		if( mycount>0 )
		{
			attributedFieldName = fieldName;

			if( mycount> 1 )
				cout<<"Warning : " << mycount << " fields are associated to the name " << fieldName <<" in file "<<completeFileName<<endl;
		}
		else {
			std::ostringstream message;
			message << "No field named " << fieldName << " in file " << completeFileName;
			throw CdmathException(message.str().c_str());
		}
	}

	// Get the name of the right mesh that we will attribute to the Field.
	std::vector<std::string> meshNames	= MEDCoupling::GetMeshNamesOnField(completeFileName, attributedFieldName);
	if (meshNames.size() == 0) {
		std::ostringstream message;
		message << "No mesh associated to " << fieldName<< " in file " << completeFileName;
		throw CdmathException(message.str().c_str());
	}
	else
		if( meshNames.size() > 1 )
		{
			cout<<"Warning : " << meshNames.size() << " meshes are associated to field named " << fieldName <<" in file "<<completeFileName<<endl;
			cout<<"Mesh names are : ";
			for (int iMesh=0; iMesh < meshNames.size(); iMesh++)
				cout<< meshNames[iMesh]<<", ";
			cout<<endl<<"Taking the mesh with name "<< meshNames[0]<<endl;			
		}

	std::string attributedMeshName = meshNames[0];

	// Create Field.
	switch (type) 
	{
		case CELLS:
			_field = dynamic_cast< MEDCoupling::MEDCouplingFieldDouble * > ( 
						MEDCoupling::ReadFieldCell( completeFileName,
						attributedMeshName, 0,
						attributedFieldName, iteration, order) );
			break;
		case NODES:
			_field = dynamic_cast< MEDCoupling::MEDCouplingFieldDouble * > (
						MEDCoupling::ReadFieldNode( completeFileName,
						attributedMeshName, 0,
						attributedFieldName, iteration, order) );
			break;
		case FACES:
			_field = dynamic_cast< MEDCoupling::MEDCouplingFieldDouble * > ( 
						MEDCoupling::ReadFieldCell( completeFileName,
						attributedMeshName, -1,
						attributedFieldName, iteration, order) );
			break;
	}

	//Read and store the number of components
	_numberOfComponents = _field->getNumberOfComponents() ;
	_time = _field->getTime(iteration, order);

	cout<<"Found field " << fieldName << " in file " << completeFileName << " on mesh named "<< _field->getMesh()->getName()<< endl;
	_mesh=Mesh( completeFileName, _field->getMesh()->getName());
}

void Field::applyFunc(  const std::string& func )
{
    _field->applyFunc( func );
}

void Field::fillFromAnalytic( int nbComp, const std::string & func )
{
    _field->fillFromAnalytic( nbComp, func ); 
}

Vector
Field::getNormEuclidean() const
{
	Vector result(_numberOfComponents);
	DoubleTab norm(_numberOfComponents,_field->magnitude()->getArray()->getConstPointer());
	result.setValues(norm);
	
	return result;
}

double
Field::max(int component) const
{
	if( component >= getNumberOfComponents() || component < 0)
		throw CdmathException("double Field::max() : component number should be between 0 and the field number of components");
		
	double result=-1e100;
	for(int i=0; i<getNumberOfElements() ; i++)
		if( result < (*this)(i,component))
			result = (*this)(i,component);

	return result;
}

double
Field::min(int component) const
{
	if( component >= getNumberOfComponents() || component < 0)
		throw CdmathException("double Field::min() : component number should be between 0 and the field number of components");
		
	double result=1e100;
	for(int i=0; i<getNumberOfElements() ; i++)
		if( result > (*this)(i,component))
			result = (*this)(i,component);

	return result;
}

string
Field::getInfoOnComponent(int icomp) const
{
	return _field->getArray()->getInfoOnComponent(icomp);
}

void
Field::setInfoOnComponent(int icomp, string nameCompo)
{
	_field.retn()->getArray()->setInfoOnComponent(icomp,nameCompo);
}

//----------------------------------------------------------------------
Field::Field( const Field & f )
//----------------------------------------------------------------------
{
	_mesh=f.getMesh() ;
	MCAuto<MEDCouplingFieldDouble> f1=f.getMEDCouplingField()->clone(true);
	_field=f1;
	_typeField=f.getTypeOfField();
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingFieldDouble>
Field::getMEDCouplingField ( void )  const
//----------------------------------------------------------------------
{
	return _field->clone(true);
}

//----------------------------------------------------------------------
void
Field::setFieldByMEDCouplingFieldDouble ( const MEDCouplingFieldDouble* field )
//----------------------------------------------------------------------
{
	MCAuto<MEDCouplingFieldDouble> ff=field->clone(true);
	_field=ff;
}

//----------------------------------------------------------------------
void
Field::setFieldByDataArrayDouble ( const DataArrayDouble* array )
//----------------------------------------------------------------------
{
	_field->setArray(const_cast<DataArrayDouble*>(array));
}

//----------------------------------------------------------------------
Vector
Field::integral ( ) const
//----------------------------------------------------------------------
{
	int nbComp=_field->getNumberOfComponents();
	double res[nbComp];//Pointer containing the inegral of each component
	_field->integral(false,res);
	Vector result(nbComp);//Vector containing the inegral of each component

	for(int i=0; i<nbComp ; i++)
		result(i)=res[i];

	return result;
}

//----------------------------------------------------------------------
double
Field::integral ( int compId ) const
//----------------------------------------------------------------------
{
	return _field->integral(compId, false);
}

//----------------------------------------------------------------------
Vector
Field::normL1 ( ) const
//----------------------------------------------------------------------
{
	int nbComp=_field->getNumberOfComponents();
	double res[nbComp];//Pointer containing the L1 norm of each component
	_field->normL1(res);
	Vector result(nbComp);//Vector containing the L1 norm of each component

	for(int i=0; i<nbComp ; i++)
		result(i)=res[i];

	return result;
}

//----------------------------------------------------------------------
Vector
Field::normL2 ( ) const
//----------------------------------------------------------------------
{
	int nbComp=_field->getNumberOfComponents();
	double res[nbComp];//Pointer containing the L2 norm of each component
	_field->normL2(res);
	Vector result(nbComp);//Vector containing the L2 norm of each component

	for(int i=0; i<nbComp ; i++)
		result(i)=res[i];

	return result;
}

//----------------------------------------------------------------------
Vector
Field::normMax ( ) const
//----------------------------------------------------------------------
{
	int nbComp=_field->getNumberOfComponents();
	double res[nbComp];//Pointer containing the L2 norm of each component
	_field->normMax(res);
	Vector result(nbComp);//Vector containing the L2 norm of each component

	for(int i=0; i<nbComp ; i++)
		result(i)=res[i];

	return result;
}

Vector 
Field::componentMax(Vector & Indices) const
{
	int nbComp=_field->getNumberOfComponents();
   	int nbElems=getNumberOfElements();

	Vector result(nbComp);//Vector containing the Linfinity norm of each component

	for(int i=0; i<nbElems ; i++)
        for(int j=0; j<nbComp ; j++)
            if(fabs((*this)(i,j))>result(j))
            {
                result(j)=fabs((*this)(i,j));
                Indices(j)=i;
            }
	return result;    
}

//----------------------------------------------------------------------
double&
Field::operator() ( int ielem )
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || ielem<0)
		throw CdmathException("double& Field::operator(ielem) : ielem>number of values !");
	return _field->getArray()->getPointer()[ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
double&
Field::operator[] ( int ielem )
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || ielem<0)
		throw CdmathException("double& Field::operator[ielem] : ielem>number of values !");
	return _field->getArray()->getPointer()[ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
double
Field::operator() ( int ielem ) const
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || ielem<0)
		throw CdmathException("double Field::operator(ielem) : ielem>number of values !");
	return _field->getArray()->getConstPointer()[ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
double
Field::operator[] ( int ielem ) const
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || ielem<0)
		throw CdmathException("double Field::operator[ielem] : ielem>number of values !");
	return _field->getArray()->getConstPointer()[ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
double&
Field::operator() ( int ielem, int jcomp )
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || jcomp>_field->getNumberOfComponents() || ielem<0 || jcomp<0)
		throw CdmathException("double& Field::operator( int ielem, int jcomp ) : ielem>number of values or jcomp>number of components !");
	return _field->getArray()->getPointer()[jcomp+ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
double
Field::operator() (  int ielem, int jcomp ) const
//----------------------------------------------------------------------
{
	if(ielem>_field->getNumberOfTuples() || jcomp>_field->getNumberOfComponents() || ielem<0 || jcomp<0)
		throw CdmathException("double Field::operator(  int ielem, int jcomp ) : ielem>number of values or jcomp>number of components !");
	return _field->getArray()->getConstPointer()[jcomp+ielem*_field->getNumberOfComponents()];
}

//----------------------------------------------------------------------
void
Field::setTime ( double time, int iter )
//----------------------------------------------------------------------
{
	_field->setTime(time,iter,0.0);
}
//----------------------------------------------------------------------
void
Field::setTimeIteration ( int iter )
//----------------------------------------------------------------------
{
	 double time=getTime();
	_field->setTime(time,iter,0.0);
}
//----------------------------------------------------------------------
double
Field::getTime ( void ) const
//----------------------------------------------------------------------
{
	int a,b;
	return _field->getTime(a,b);
}

//----------------------------------------------------------------------
int
Field::getTimeIteration ( void ) const
//----------------------------------------------------------------------
{
	int a,b;
	_field->getTime(a,b);
	return a;
}

//----------------------------------------------------------------------
int
Field::getNumberOfElements ( void ) const
//----------------------------------------------------------------------
{
	return _field->getNumberOfTuples() ;
}

int
Field::getSpaceDimension( void ) const
{
	return _mesh.getSpaceDimension() ;
}

//----------------------------------------------------------------------
int
Field::getNumberOfComponents ( void ) const
//----------------------------------------------------------------------
{
	return _field->getNumberOfComponents() ;
}

//----------------------------------------------------------------------
const double*
Field::getValues ( void ) const
//----------------------------------------------------------------------
{
	return _field->getArray()->getConstPointer() ;
}

//----------------------------------------------------------------------
void
Field::getValues ( Vector myVector ) const
//----------------------------------------------------------------------
{
	if( myVector.size() != _field->getNumberOfTuples() * _field->getNumberOfComponents() )
		throw CdmathException("Error : Field::getValues requires a vector having the same number of component as fiedl values");

    const double * fieldValues = _field->getArray()->getConstPointer();
	double * vectorValues = myVector.getValues().getValues();
    
	memcpy(vectorValues, fieldValues,myVector.size()*sizeof(double)) ;	
}

void 
Field::setValues ( Vector myVector )
//----------------------------------------------------------------------
{
	if( myVector.size() != _field->getNumberOfTuples() * _field->getNumberOfComponents() )
		throw CdmathException("Error : Field::setValues requires a vector having the same number of component as fiedl values");
		
	double * vectorValues = myVector.getValues().getValues();

	setValues ( vectorValues, myVector.size() );
}

void 
Field::setValues ( double * values, int nbValues )
{
	double * fieldValues = _field->getArray()->getPointer() ;
	memcpy(fieldValues,values,nbValues*sizeof(double)) ;	
}

//----------------------------------------------------------------------
const string
Field::getName ( void ) const
//----------------------------------------------------------------------
{
	return _field->getName() ;
}

//----------------------------------------------------------------------
const Mesh&
Field::getMesh ( void ) const
//----------------------------------------------------------------------
{
	return _mesh ;
}

//----------------------------------------------------------------------
EntityType
Field::getTypeOfField ( void ) const
//----------------------------------------------------------------------
{
	return _typeField;
}

double 
Field::getElementComponent(int i, int comp) const
{
	switch( _typeField )
	{
		case CELLS:
			switch( comp )
			{
				case 0:
					return _mesh.getCell(i).x();
				case 1:
					return _mesh.getCell(i).y();
				case 2:
					return _mesh.getCell(i).z();
				default:
					cout<<"Wrong component number "<< comp <<" , dimension is "<< _mesh.getSpaceDimension() << ", field values are on CELLS" <<endl;
					throw CdmathException("Field::getElementComponent : Wrong component number");
			}
		case NODES:
			switch( comp )
			{
				case 0:
					return _mesh.getNode(i).x();
				case 1:
					return _mesh.getNode(i).y();
				case 2:
					return _mesh.getNode(i).z();
				default:
					cout<<"Wrong component number "<< comp <<" , dimension is "<< _mesh.getSpaceDimension() << ", field values are on NODES" <<endl;
					throw CdmathException("Field::getElementComponent : Wrong component number");
			}
		case FACES:
			switch( comp )
			{
				case 0:
					return _mesh.getFace(i).x();
				case 1:
					return _mesh.getFace(i).y();
				case 2:
					return _mesh.getFace(i).z();
				default:
					cout<<"Wrong component number "<< comp <<" , dimension is "<< _mesh.getSpaceDimension() << ", field values are on FACES"<<endl;
					throw CdmathException("Field::getElementComponent : Wrong component number");
			}
		default:
			throw CdmathException("Accepted field supports are CELLS, NODES and FACES");
	}
	
}
//----------------------------------------------------------------------
void
Field::setName ( const string fieldName )
//----------------------------------------------------------------------
{
	_field->setName(fieldName.c_str()) ;
}


//----------------------------------------------------------------------
Field
Field::operator+ ( const Field& f ) const
//----------------------------------------------------------------------
{
	//if(f.getMesh().getMEDCouplingMesh() != _mesh.getMEDCouplingMesh())
	//throw CdmathException("Field::operator+ : Field addition requires identical meshes");
	_mesh.getMEDCouplingMesh()->checkFastEquivalWith(f.getMesh().getMEDCouplingMesh(),1e-6);
	if(f.getTypeOfField() != getTypeOfField())
		throw CdmathException("Field::operator+ : Field addition requires identical field types (CELLS, NODES or FACES");
	if(f.getNumberOfComponents() != getNumberOfComponents())
		throw CdmathException("Field::operator+ : Field addition requires identical number of components");

	Field fres(getName(),f.getTypeOfField(),f.getMesh(),f.getNumberOfComponents(),f.getTime());
	int nbComp=f.getNumberOfComponents();
	int nbElem=f.getNumberOfElements();
	for (int ielem=0 ; ielem<nbElem; ielem++)
		for (int jcomp=0 ; jcomp<nbComp ; jcomp++)
			fres(ielem, jcomp)=_field->getArray()->getConstPointer()[jcomp+ielem*nbComp]+f(ielem, jcomp);
	return fres;
}

//----------------------------------------------------------------------
Field
Field::operator- ( const Field& f ) const
//----------------------------------------------------------------------
{
	//if(f.getMesh().getMEDCouplingMesh() != _mesh.getMEDCouplingMesh())
	//throw CdmathException("Field::operator- : Field subtraction requires identical meshes");
	_mesh.getMEDCouplingMesh()->checkFastEquivalWith(f.getMesh().getMEDCouplingMesh(),1e-6);
	if(f.getTypeOfField() != getTypeOfField())
		throw CdmathException("Field::operator- : Field subtraction requires identical field types (CELLS, NODES or FACES");
	if(f.getNumberOfComponents() != getNumberOfComponents())
		throw CdmathException("Field::operator- : Field subtraction requires identical number of components");

	Field fres(getName(),f.getTypeOfField(),f.getMesh(),f.getNumberOfComponents(),f.getTime());
	int nbComp=f.getNumberOfComponents();
	int nbElem=f.getNumberOfElements();
	for (int ielem=0 ; ielem<nbElem; ielem++)
		for (int jcomp=0 ; jcomp<nbComp ; jcomp++)
			fres(ielem, jcomp)=_field->getArray()->getConstPointer()[jcomp+ielem*nbComp]-f(ielem, jcomp);
	return fres;
}

//----------------------------------------------------------------------
const Field&
Field::operator= ( const Field& f )
//----------------------------------------------------------------------
{
	_mesh=f.getMesh() ;
	_typeField=f.getTypeOfField() ;
	_numberOfComponents=f.getNumberOfComponents();
	_time=f.getTime();
	_fieldName=f.getName();
	MCAuto<MEDCouplingFieldDouble> f1=f.getMEDCouplingField()->clone(true);
	_field=f1;
	return *this;
}

Field Field::deepCopy( ) const
{
    Field F(getName(), getTypeOfField(), getMesh(), getNumberOfComponents(), getTime()) ;
	MCAuto<MEDCouplingFieldDouble> f1=getMEDCouplingField()->clone(true);
	F.setFieldByMEDCouplingFieldDouble(f1);
    
    return F;
}


//----------------------------------------------------------------------
const Field&
Field::operator+= ( const Field& f )
//----------------------------------------------------------------------
{
	//if(f.getMesh().getMEDCouplingMesh() != _mesh.getMEDCouplingMesh())
	//throw CdmathException("Field::operator+= : Field addition requires identical meshes");
	_mesh.getMEDCouplingMesh()->checkFastEquivalWith(f.getMesh().getMEDCouplingMesh(),1e-6);
	if(f.getTypeOfField() != getTypeOfField())
		throw CdmathException("Field::operator+= : Field addition requires identical field types (CELLS, NODES or FACES");
	if(f.getNumberOfComponents() != getNumberOfComponents())
		throw CdmathException("Field::operator+= : Field addition requires identical number of components");

	_field->setMesh(f.getMEDCouplingField()->getMesh());
	(*_field)+=(*f.getMEDCouplingField());
	return *this;
}

//----------------------------------------------------------------------
const Field&
Field::operator-= ( const Field& f )
//----------------------------------------------------------------------
{
	//if(f.getMesh().getMEDCouplingMesh() != _mesh.getMEDCouplingMesh())
	//throw CdmathException("Field::operator-= : Field subtraction requires identical meshes");
	_mesh.getMEDCouplingMesh()->checkFastEquivalWith(f.getMesh().getMEDCouplingMesh(),1e-6);
	if(f.getTypeOfField() != getTypeOfField())
		throw CdmathException("Field::operator-= : Field subtraction requires identical field types (CELLS, NODES or FACES");
	if(f.getNumberOfComponents() != getNumberOfComponents())
		throw CdmathException("Field::operator-= : Field subtraction requires identical number of components");

	_field->setMesh(f.getMEDCouplingField()->getMesh());
	(*_field)-=(*f.getMEDCouplingField());
	return *this;
}

//----------------------------------------------------------------------
const Field&
Field::operator*= ( double s )
//----------------------------------------------------------------------
{
	int nbComp=getNumberOfComponents();
	int nbElem=getNumberOfElements();
	for (int i=0 ; i<nbComp ; i++)
		for (int j=0 ; j<nbElem; j++)
			_field->getArray()->getPointer()[i+j*nbComp]*=s;
	return *this;
}

//----------------------------------------------------------------------
const Field&
Field::operator/= ( double s )
//----------------------------------------------------------------------
{
	int nbComp=getNumberOfComponents();
	int nbElem=getNumberOfElements();
	for (int i=0 ; i<nbComp ; i++)
		for (int j=0 ; j<nbElem; j++)
			_field->getArray()->getPointer()[i+j*nbComp]/=s;
	return *this;
}

//----------------------------------------------------------------------
const Field&
Field::operator-= ( double s )
//----------------------------------------------------------------------
{
	int nbComp=getNumberOfComponents();
	int nbElem=getNumberOfElements();
	for (int i=0 ; i<nbComp ; i++)
		for (int j=0 ; j<nbElem; j++)
			_field->getArray()->getPointer()[i+j*nbComp]-=s;
	return *this;
}

//----------------------------------------------------------------------
const Field&
Field::operator+= ( double s )
//----------------------------------------------------------------------
{
	int nbComp=getNumberOfComponents();
	int nbElem=getNumberOfElements();
	for (int i=0 ; i<nbComp ; i++)
		for (int j=0 ; j<nbElem; j++)
			_field->getArray()->getPointer()[i+j*nbComp]+=s;
	return *this;
}

//----------------------------------------------------------------------
void
Field::writeVTK (std::string fileName, bool fromScratch) const
//----------------------------------------------------------------------
{
	string fname=fileName+".pvd";
	int iter,order;
	double time=_field->getTime(iter,order);

	if (fromScratch)
	{
		ofstream file(fname.c_str()) ;
		file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n" ;
		ostringstream numfile;
		numfile << iter ;
		string filetmp=fileName+"_";
		filetmp=filetmp+numfile.str();
		string ret=_field->writeVTK(filetmp.c_str()) ;
		file << "<DataSet timestep=\""<< time << "\" group=\"\" part=\"0\" file=\"" << ret << "\"/>\n" ;
		file << "</Collection></VTKFile>\n" ;
		file.close() ;
	}
	else
	{
		ifstream file1(fname.c_str()) ;
		string contenus;
		getline(file1, contenus, '\0');
		string to_remove="</Collection></VTKFile>";
		size_t m = contenus.find(to_remove);
		size_t n = contenus.find_first_of("\n", m + to_remove.length());
		contenus.erase(m, n - m + 1);
		file1.close() ;
		ofstream file(fname.c_str()) ;
		file << contenus ;
		ostringstream numfile;
		numfile << iter ;
		string filetmp=fileName+"_";
		filetmp=filetmp+numfile.str();
		string ret=_field->writeVTK(filetmp.c_str()) ;
		file << "<DataSet timestep=\""<< time << "\" group=\"\" part=\"0\" file=\"" << ret << "\"/>\n" ;
		file << "</Collection></VTKFile>\n" ;
		file.close() ;
	}
}

//----------------------------------------------------------------------
void
Field::writeCSV ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	int iter,order;
	double time=_field->getTime(iter,order);

	ostringstream numfile;
	numfile << iter ;
	string filetmp=fileName+"_";
	filetmp=filetmp+numfile.str();
	filetmp=filetmp+".csv";
	ofstream file(filetmp.c_str()) ;
	int dim=_mesh.getSpaceDimension();
	int nbElements;
	if (getTypeOfField()==CELLS)
		nbElements=_mesh.getNumberOfCells();
	else
		nbElements=_mesh.getNumberOfNodes();

	if (dim==1)
	{
		int nbCompo=getNumberOfComponents();
		if (nbCompo==1)
			file << "x " << _field->getName() << endl;
		else if (nbCompo>1)
		{
			file << "x";
			for (int i=0;i<nbCompo;i++)
				file << " " << _field->getName() << " Compo " << i+1 << " "<< getInfoOnComponent(i);
			file << endl;
		}
		for (int i=0;i<nbElements;i++)
		{
			if (getTypeOfField()==CELLS)
				file << _mesh.getCell(i).x() ;
			else
				file << _mesh.getNode(i).x() ;
			for (int j=0;j<nbCompo;j++)
				file << " " << getValues()[j+i*nbCompo] ;
			file << endl;
		}
	}else if (dim==2)
	{
		int nbCompo=getNumberOfComponents();
		if (nbCompo==1)
			file << "x y " << _field->getName() << endl;
		else if (nbCompo>1)
		{
			file << "x y";
			for (int i=0;i<nbCompo;i++)
				file << " " << _field->getName() << " Compo " << i+1<< " "<< getInfoOnComponent(i);
			file << endl;
		}
		for (int i=0;i<nbElements;i++)
		{
			if (getTypeOfField()==CELLS)
				file << _mesh.getCell(i).x() << " " << _mesh.getCell(i).y() ;
			else
				file << _mesh.getNode(i).x() << " " << _mesh.getNode(i).y() ;
			for (int j=0;j<nbCompo;j++)
				file << " " << getValues()[j+i*nbCompo] ;
			file << endl;
		}
	}else
	{
		int nbCompo=getNumberOfComponents();
		if (nbCompo==1)
			file << "x y z " << _field->getName() << endl;
		else if (nbCompo>1)
		{
			file << "x y z";
			for (int i=0;i<nbCompo;i++)
				file << " " << _field->getName() << " Compo " << i+1<< " "<< getInfoOnComponent(i);
			file << endl;
		}
		for (int i=0;i<nbElements;i++)
		{
			if (getTypeOfField()==CELLS)
				file << _mesh.getCell(i).x() << " " << _mesh.getCell(i).y() << " " << _mesh.getCell(i).z();
			else
				file << _mesh.getNode(i).x() << " " << _mesh.getNode(i).y() << " " << _mesh.getNode(i).z();
			for (int j=0;j<nbCompo;j++)
				file << " " << getValues()[j+i*nbCompo] ;
			file << endl;
		}
	}
	file.close() ;
}

//----------------------------------------------------------------------
void
Field::writeMED ( const std::string fileName, bool fromScratch)
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	
	if(_mesh.isStructured() || _mesh.meshNotDeleted())
		if (fromScratch)
			MEDCoupling::WriteField(fname.c_str(),_field,fromScratch);
		else
			MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(fname.c_str(),_field);
	else//The mesh has ben deleted, use a MEDFileField1TS instead of _field to save the values
	{
		if ( not fromScratch)
		{
			MEDCoupling::MCAuto<MEDCoupling::MEDFileField1TS> ff=MEDFileField1TS::New();//To save the field when the mesh has been deleted
			ff->setFieldNoProfileSBT(  _field );
			ff->write(fname.c_str(), fromScratch);
		}
		else
			throw CdmathException("Field::writeMED Error !!! The mesh has been deleted, cannot write field from scratch");
	}
}

Field
operator* (double value , const Field& field )
{
	Field fres(field.getName(),field.getTypeOfField(),field.getMesh(),field.getNumberOfComponents(),field.getTime());
	int nbComp=field.getNumberOfComponents();
	int nbElem=field.getNumberOfElements();
	for (int ielem=0 ; ielem<nbElem; ielem++)
		for (int jcomp=0 ; jcomp<nbComp ; jcomp++)
			fres(ielem, jcomp)=value*field(ielem, jcomp);
	return fres;
}

Field
operator* (const Field& field, double value )
{
	Field fres(field.getName(),field.getTypeOfField(),field.getMesh(),field.getNumberOfComponents(),field.getTime());
	int nbComp=field.getNumberOfComponents();
	int nbElem=field.getNumberOfElements();
	for (int ielem=0 ; ielem<nbElem; ielem++)
		for (int jcomp=0 ; jcomp<nbComp ; jcomp++)
			fres(ielem, jcomp)=value*field(ielem, jcomp);
	return fres;
}

Field operator/ (const Field& field, double value)
				{
	Field fres(field.getName(),field.getTypeOfField(),field.getMesh(),field.getNumberOfComponents(),field.getTime());
	int nbComp=field.getNumberOfComponents();
	int nbElem=field.getNumberOfElements();
	for (int ielem=0 ; ielem<nbElem; ielem++)
		for (int jcomp=0 ; jcomp<nbComp ; jcomp++)
			fres(ielem, jcomp)=field(ielem, jcomp)/value;
	return fres;
				}

Vector
Field::getValuesOnAllComponents(int elem) const
{
	Vector v(getNumberOfComponents());
	for(int i=0;i<getNumberOfComponents();i++)
		v[i]=(*this)(elem,i);
	return v;
}

Vector
Field::getValuesOnComponent(int compo) const
{
	Vector v(getNumberOfElements());
	for(int i=0;i<getNumberOfElements();i++)
		v[i]=(*this)(i,compo);
	return v;
}

std::vector< double > 
Field::getFieldValues(int compo) const 
{
	std::vector< double > v(getNumberOfElements());
	for(int i=0;i<getNumberOfElements();i++)
		v[i]=(*this)(i,compo);
	return v;
}

std::ostream& operator<<(std::ostream& out, const Field& field )
{
	cout << "Field " << field.getName() << " : " << endl ;
	out<< field.getMEDCouplingField().retn()->getArray()->repr();
	return out;
}

void Field::deleteMEDCouplingMesh()
{ 
	return _mesh.deleteMEDCouplingMesh();
}
