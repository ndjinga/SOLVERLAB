/*
 * field.cxx
 *
 *  Created on: 7 fevrier. 2012
 *      Author: CDMAT
 */

#include "Node.hxx"
#include "Cell.hxx"
#include "Field.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "CdmathException.hxx"

#include <fstream>
#include <sstream>
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
		array->alloc(_mesh.getNumberOfNodes(),_numberOfComponents);
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
	}else
		throw CdmathException("Type of Field::Field() is not compatible");

	_field->setName(_fieldName.c_str()) ;
	_field->setArray(array);
	_field->setTime(_time,0,0);
	array->decrRef();
	mu->decrRef();
}

Field::Field( const std::string filename, EntityType type,
		const std::string & fieldName,
		int iteration, int order, int meshLevel,
		int numberOfComponents, double time)
{
	_field = NULL;
	_mesh=Mesh(filename + ".med", meshLevel);
	_typeField=type;
	_numberOfComponents=numberOfComponents;
	_time=time;
	_fieldName=fieldName;

	readFieldMed(filename, type, fieldName, iteration, order);
}

Field::Field(const std::string meshFileName, EntityType type, const std::vector<double> Vconstant, 
		const std::string & fieldName, int meshLevel, double time )
{
	_field = NULL;
	_mesh=Mesh(meshFileName + ".med", meshLevel);
	_typeField=type;
	_numberOfComponents=Vconstant.size();
	_time=time;
	_fieldName=fieldName;

	buildFieldMemoryStructure();

	int nbelem=_field->getNumberOfTuples();
	int nbcomp=_field->getNumberOfComponents() ;

	for (int ielem=0 ; ielem<nbelem; ielem++)
		for (int jcomp=0 ; jcomp<nbcomp ; jcomp++)
			_field->getArray()->getPointer()[jcomp+ielem*_field->getNumberOfComponents()]=Vconstant[jcomp];
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
			_field->getArray()->getPointer()[jcomp+ielem*_field->getNumberOfComponents()]=Vconstant[jcomp];
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
			_field->getArray()->getPointer()[jcomp+ielem*_field->getNumberOfComponents()]=Vconstant[jcomp];
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
			_field->getArray()->getPointer()[jcomp+ielem*_field->getNumberOfComponents()]=Vconstant[jcomp];
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
				_field->getArray()->getPointer()[j+i*_field->getNumberOfComponents()] = VV_Left[i];
			else
				_field->getArray()->getPointer()[j+i*_field->getNumberOfComponents()] = VV_Right[i];
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
				_field->getArray()->getPointer()[j+i*_field->getNumberOfComponents()]=Vin[j];
		else
			for(int j=0;j<nbcomp;j++)
				_field->getArray()->getPointer()[j+i*_field->getNumberOfComponents()]=Vout[j];
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
			attributedFieldName = fieldNames[0];
		else {
			std::ostringstream message;
			message << "No field in file " << completeFileName;
			throw CdmathException(message.str().c_str());
		}
	}
	else {
		for (; iField < fieldNames.size(); iField++)
			if (fieldName == fieldNames[iField]) break;

		if (iField < fieldNames.size())
			attributedFieldName = fieldName;
		else {
			std::ostringstream message;
			message << "No field named " << fieldName << " in file " << completeFileName;
			throw CdmathException(message.str().c_str());
		}
	}

	// Get the name of the right mesh that we will attribute to the Field.
	std::vector<std::string> meshNames
	= MEDCoupling::GetMeshNamesOnField(completeFileName, attributedFieldName);
	if (meshNames.size() == 0) {
		std::ostringstream message;
		message << "No mesh associated to " << fieldName
				<< " in file " << completeFileName;
		throw CdmathException(message.str().c_str());
	}
	std::string attributedMeshName = meshNames[0];

	// Create Field.
	MEDCoupling::TypeOfField medFieldType[3] = { ON_CELLS, ON_NODES, ON_CELLS };
	switch (type) {
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
}


DoubleTab
Field::getNormEuclidean() const
{
	DoubleTab norm(getNumberOfElements(),_field->magnitude()->getArray()->getConstPointer());
	return norm;
}

double
Field::max() const
{
	if( getNumberOfComponents() !=1)
		throw CdmathException("double Field::max() : field should have a single component in order to extract maximum value");
		
	double result=0;
	for(int i=0; i<getNumberOfElements() ; i++)
		if( result < (*this)(i,0))
			result = (*this)(i,0);

	return result;
}

double
Field::min() const
{
	if( getNumberOfComponents() !=1)
		throw CdmathException("double Field::min() : field should have a single component in order to extract minimum value");
		
	double result=0;
	for(int i=0; i<getNumberOfElements() ; i++)
		if( result > (*this)(i,0))
			result = (*this)(i,0);

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
	MCAuto<MEDCouplingFieldDouble> f1=f.getField()->deepCopy();
	_field=f1;
	_typeField=f.getTypeOfField();
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingFieldDouble>
Field::getField ( void )  const
//----------------------------------------------------------------------
{
	return _field ;
}

//----------------------------------------------------------------------
void
Field::setFieldByMEDCouplingFieldDouble ( const MEDCouplingFieldDouble* field )
//----------------------------------------------------------------------
{
	MCAuto<MEDCouplingFieldDouble> ff=field->deepCopy();
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
   	int nbElems=getNumberOfElements();

	Vector result(nbComp);//Vector containing the Linfinity norm of each component

	for(int i=0; i<nbElems ; i++)
        for(int j=0; j<nbComp ; j++)
            if(fabs((*this)(i,j))>result(j))
                result(j)=fabs((*this)(i,j));

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
double
Field::getTime ( void ) const
//----------------------------------------------------------------------
{
	int a,b;
	return _field->getTime(a,b);
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
			fres(ielem, jcomp)=_field->getArray()->getConstPointer()[jcomp+ielem*_field->getNumberOfComponents()]+f(ielem, jcomp);
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
			fres(ielem, jcomp)=_field->getArray()->getConstPointer()[jcomp+ielem*_field->getNumberOfComponents()]-f(ielem, jcomp);
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
	MCAuto<MEDCouplingFieldDouble> f1=f.getField()->deepCopy();
	_field=f1;
	return *this;
}

Field Field::deepCopy( ) const
{
    Field F(getName(), getTypeOfField(), getMesh(), getNumberOfComponents(), getTime()) ;
	MCAuto<MEDCouplingFieldDouble> f1=getField()->deepCopy();
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

	_field->setMesh(f.getField()->getMesh());
	(*_field)+=(*f.getField());
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

	_field->setMesh(f.getField()->getMesh());
	(*_field)-=(*f.getField());
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
			_field->getArray()->getPointer()[i+j*_field->getNumberOfComponents()]*=s;
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
			_field->getArray()->getPointer()[i+j*_field->getNumberOfComponents()]/=s;
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
			_field->getArray()->getPointer()[i+j*_field->getNumberOfComponents()]-=s;
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
			_field->getArray()->getPointer()[i+j*_field->getNumberOfComponents()]+=s;
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
Field::writeMED ( const std::string fileName, bool fromScratch) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	if (fromScratch)
		MEDCoupling::WriteField(fname.c_str(),_field,fromScratch);
	else
		MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(fname.c_str(),_field);
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

std::ostream& operator<<(std::ostream& out, const Field& field )
{
	cout << "Field " << field.getName() << " : " << endl ;
	out<< field.getField().retn()->getArray()->repr();
	return out;
}
