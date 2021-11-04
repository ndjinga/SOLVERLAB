//============================================================================
// Name        : Problème CoreFlows
// Author      : M. Ndjinga
// Version     :
// Copyright   : CEA Saclay 2014
// Description : Classe générique des problèmes de Thermohydraulique
//			     Il s'agit d'une classe virtuelle pure (non instanciable qui regroupe
//				les fonctionalités communes à tous les codes thermohydraulique.
//				Le but est d'avoir une interface commune rendant
//				plus lisibles et plus evolutifs les codes dévelopés avec CDMATH
//============================================================================

#include "ProblemCoreFlows.hxx"
#include "SparseMatrixPetsc.hxx"

#include <limits.h>
#include <unistd.h>

using namespace std;

ProblemCoreFlows::ProblemCoreFlows(MPI_Comm comm)
{
	/* Initialisation of PETSC */
	//check if PETSC is already initialised
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
	{//check if MPI is already initialised
		int mpiInitialized;
		MPI_Initialized(&mpiInitialized);
		if(mpiInitialized)
			PETSC_COMM_WORLD = comm;
		PetscInitialize(NULL,NULL,0,0);//Note this is ok if MPI has been been initialised independently from PETSC
	}
	MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&_size);
	PetscPrintf(PETSC_COMM_WORLD,"Simulation on %d processors\n",_size);//Prints to standard out, only from the first processor in the communicator. Calls from other processes are ignored. 
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Processor [%d] ready for action\n",_rank);//Prints synchronized output from several processors. Output of the first processor is followed by that of the second, etc. 
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

	if(_size>1)
	{
		PetscPrintf(PETSC_COMM_WORLD,"---- More than one processor detected : running a parallel simulation ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Limited parallelism : input and output fields remain sequential ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Only the matrixoperations are done in parallel thanks to PETSc----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Processor %d is in charge of building the mesh, saving the results, filling and then distributing the matrix to other processors.\n\n",_rank);
	}
	
	/* Numerical parameter */
	_dt = 0;
	_time = 0;
	_nbTimeStep=0;
	_cfl=0;
	_timeMax =0.;
	_maxNbOfTimeStep = 0;
	_precision=1.e-6;
	_precision_Newton=_precision;
	_erreur_rel= 0;
	_isStationary=false;
	_timeScheme=Explicit;
	_wellBalancedCorrection=false;
    _FECalculation=false;
	_nonLinearSolver = Newton_SOLVERLAB;
	_conditionNumber=false;
	_maxvp=0;
	_runLogFile=new ofstream;

	/* Monitoring of simulation */
	_restartWithNewTimeScheme=false;
	_restartWithNewFileName=false;
	_fileName = "myCoreFlowsProblem";
	_freqSave = 1;
	_verbose = false;
	_system = false;

	/* Mesh parameters */
	_Ndim=0;
	_minl=0;
	_neibMaxNb=0;

	/* Memory and restart */
	_initialDataSet=false;
	_initializedMemory=false;

	/* PETSc and linear solver parameters */
	_MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
	_NEWTON_its=0;
	_maxPetscIts=50;
	_maxNewtonIts=50;
	_PetscIts=0;//the number of iterations of the linear solver
	_ksptype = (char*)&KSPGMRES;
	_pctype = (char*)&PCLU;

	/* Physical parameter */
	_heatPowerFieldSet=false;
	_heatTransfertCoeff=0;
	_rodTemperatureFieldSet=false;
	_heatSource=0;

	//extracting current directory
	char result[ PATH_MAX ];
	getcwd(result, PATH_MAX );
	_path=string( result );
	_saveFormat=VTK;
}

TimeScheme ProblemCoreFlows::getTimeScheme()
{
	return _timeScheme;
}

void ProblemCoreFlows::setTimeScheme(TimeScheme timeScheme)
{
	if( _nbTimeStep>0 && timeScheme!=_timeScheme)//This is a change of time scheme during a simulation
		_restartWithNewTimeScheme=true;
	_timeScheme = timeScheme;
}

bool ProblemCoreFlows::isStationary() const
{
	return _isStationary;
}

double ProblemCoreFlows::presentTime() const
{
	return _time;
}
void ProblemCoreFlows::setTimeMax(double timeMax){
	_timeMax = timeMax;
}
void ProblemCoreFlows::setPresentTime (double time)
{
	_time=time;
}
void ProblemCoreFlows::setMaxNbOfTimeStep(int maxNbOfTimeStep){
	_maxNbOfTimeStep = maxNbOfTimeStep;
}
void ProblemCoreFlows::setCFL(double cfl)
{
	_cfl=cfl;
}
void ProblemCoreFlows::setPrecision(double precision)
{
	_precision=precision;
}
void ProblemCoreFlows::setInitialField(const Field &VV)
{

	if(_Ndim != VV.getSpaceDimension()){
		*_runLogFile<<"ProblemCoreFlows::setInitialField: mesh has incorrect space dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::setInitialField: mesh has incorrect space dimension");
	}
	if(_nVar!=VV.getNumberOfComponents())
	{
		*_runLogFile<<"ProblemCoreFlows::setInitialField: Initial field has incorrect number of components"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::setInitialField: Initial field has incorrect number of components");
	}
	if(_FECalculation && VV.getTypeOfField()!=NODES)
	{
		*_runLogFile<<"ProblemCoreFlows::setInitialField: Initial field has incorrect support : should be on nodes for the Finite Elements method"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::setInitialField: Initial field has incorrect support : should be on nodes for the Finite Elements method");
	}
	else if(!_FECalculation && VV.getTypeOfField()==NODES)
	{
		*_runLogFile<<"ProblemCoreFlows::setInitialField: Initial field has incorrect support : should be on cells or faces for the Finite Volumes method"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::setInitialField: Initial field has incorrect support : should be on cells or faces for the Finite Volumes method");
	}

	_VV=VV;
	_VV.setName("SOLVERLAB results");
	_time=_VV.getTime();
	_mesh=_VV.getMesh();
	_initialDataSet=true;

	//Mesh data
	_Nmailles = _mesh.getNumberOfCells();
	_Nnodes =   _mesh.getNumberOfNodes();
	_Nfaces =   _mesh.getNumberOfFaces();
	_perimeters=Field("Perimeters", CELLS, _mesh,1);

	// find _minl (delta x) and maximum nb of neibourghs
	_minl  = INFINITY;
	int nbNeib,indexFace;
	Cell Ci;
	Face Fk;

	if(_verbose)
		PetscPrintf(PETSC_COMM_WORLD,"Computing cell perimeters and mesh minimal diameter\n");

	//Compute the maximum number of neighbours for nodes or cells
    if(VV.getTypeOfField()==NODES)
        _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
    else
        _neibMaxNb=_mesh.getMaxNbNeighbours(CELLS);
        
	//Compute Delta x and the cell perimeters
	for (int i=0; i<_mesh.getNumberOfCells(); i++){
		Ci = _mesh.getCell(i);
		if (_Ndim > 1){
			_perimeters(i)=0;
			for (int k=0 ; k<Ci.getNumberOfFaces() ; k++){
				indexFace=Ci.getFacesId()[k];
				Fk = _mesh.getFace(indexFace);
				_minl = min(_minl,Ci.getMeasure()/Fk.getMeasure());
				_perimeters(i)+=Fk.getMeasure();
			}
		}else{
			_minl = min(_minl,Ci.getMeasure());
			_perimeters(i)=Ci.getNumberOfFaces();
		}
	}
	if(_verbose)
		cout<<_perimeters<<endl;
}
//Function needed because swig of enum EntityType fails
void ProblemCoreFlows::setInitialField(string fileName, string fieldName, int timeStepNumber, int field_support_type)
{
	if(_FECalculation && field_support_type!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	Field VV;
	
	switch(field_support_type)
	{
	case CELLS:
		VV = Field(fileName, CELLS, fieldName, timeStepNumber, 0);
		break;
	case NODES:
		VV = Field(fileName, NODES, fieldName, timeStepNumber, 0);
		break;
	case FACES:
		VV = Field(fileName, FACES, fieldName, timeStepNumber, 0);
		break;
	default:
		std::ostringstream message;
		message << "Error ProblemCoreFlows::setInitialField \n Accepted field support integers are "<< CELLS <<" (for CELLS), "<<NODES<<" (for NODES), and "<< FACES <<" (for FACES)" ;
		throw CdmathException(message.str().c_str());
	}	

	setInitialField(VV);
}
//Function needed because swig of enum EntityType fails
void ProblemCoreFlows::setInitialFieldConstant( int nDim, const vector<double> Vconstant, double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide, int field_support_type)
{	
	if(_FECalculation && field_support_type!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	EntityType typeField;
	
	switch(field_support_type)
	{
	case CELLS:
		typeField=CELLS;
		break;
	case NODES:
		typeField=NODES;
		break;
	case FACES:
		typeField=FACES;
		break;
	default:
		std::ostringstream message;
		message << "Error ProblemCoreFlows::setInitialField \n Accepted field support integers are "<< CELLS <<" (for CELLS), "<<NODES<<" (for NODES), and "<< FACES <<" (for FACES)" ;
		throw CdmathException(message.str().c_str());
	}	
	
	setInitialFieldConstant( nDim, Vconstant, xmin, xmax, nx, leftSide, rightSide, ymin, ymax, ny, backSide, frontSide, zmin, zmax, nz, bottomSide, topSide, typeField);
}
void ProblemCoreFlows::setInitialField(string fileName, string fieldName, int timeStepNumber, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	Field VV(fileName, typeField, fieldName, timeStepNumber, 0);
	
	setInitialField(VV);
}
void ProblemCoreFlows::setInitialFieldConstant(string fileName, const vector<double> Vconstant, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;
	Mesh M(fileName);
	Field VV("SOLVERLAB results", typeField, M, Vconstant.size());

	for (int j = 0; j < VV.getNumberOfElements(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant[i];
	}

	setInitialField(VV);
}
void ProblemCoreFlows::	setInitialFieldConstant(const Mesh& M, const Vector Vconstant, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	Field VV("SOLVERLAB results", typeField, M, Vconstant.getNumberOfRows());

	for (int j = 0; j < VV.getNumberOfElements(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant(i);
	}
	setInitialField(VV);
}
void ProblemCoreFlows::	setInitialFieldConstant(const Mesh& M, const vector<double> Vconstant, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	Field VV("SOLVERLAB results", typeField, M, Vconstant.size());

	for (int j = 0; j < VV.getNumberOfElements(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant[i];
	}
	setInitialField(VV);
}
void ProblemCoreFlows::setInitialFieldConstant( int nDim, const vector<double> Vconstant, double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide, EntityType typeField)
{
	Mesh M;
	if(nDim==1)
		M=Mesh(xmin,xmax,nx);
	else if(nDim==2)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny);
	else if(nDim==3)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
	else{
		PetscPrintf(PETSC_COMM_WORLD,"ProblemCoreFlows::setInitialFieldConstant: Space dimension nDim should be between 1 and 3\n");
		*_runLogFile<<"ProblemCoreFlows::setInitialFieldConstant: Space dimension nDim should be between 1 and 3"<<endl;
		_runLogFile->close();
		throw CdmathException("Space dimension nDim should be between 1 and 3");
	}

	M.setGroupAtPlan(xmax,0,_precision,rightSide);
	M.setGroupAtPlan(xmin,0,_precision,leftSide);
	if(nDim>=2){
		M.setGroupAtPlan(ymax,1,_precision,frontSide);
		M.setGroupAtPlan(ymin,1,_precision,backSide);
	}
	if(nDim==3){
		M.setGroupAtPlan(zmax,2,_precision,topSide);
		M.setGroupAtPlan(zmin,2,_precision,bottomSide);
	}

	setInitialFieldConstant(M, Vconstant, typeField);
}
void ProblemCoreFlows::setInitialFieldStepFunction(const Mesh M, const Vector VV_Left, const Vector VV_Right, double disc_pos, int direction, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	if  (VV_Right.getNumberOfRows()!=VV_Left.getNumberOfRows())
	{
		*_runLogFile<<"ProblemCoreFlows::setStepFunctionInitialField: Vectors VV_Left and VV_Right have different sizes"<<endl;
		_runLogFile->close();
		throw CdmathException( "ProblemCoreFlows::setStepFunctionInitialField: Vectors VV_Left and VV_Right have different sizes");
	}
	Field VV("SOLVERLAB results", typeField, M, VV_Left.getNumberOfRows());

	double component_value;

	for (int j = 0; j < VV.getNumberOfElements(); j++) 
	{
		if(direction<=2)
			component_value=VV.getElementComponent(j, direction);
		else
		{
			PetscPrintf(PETSC_COMM_WORLD,"Error : space dimension is %d,  direction asked is \%d \n",M.getSpaceDimension(),direction);
			_runLogFile->close();
			throw CdmathException( "ProblemCoreFlows::setStepFunctionInitialField: direction should be an integer between 0 and 2");
		}

		for (int i=0; i< VV.getNumberOfComponents(); i++)
			if (component_value< disc_pos )
				VV(j,i) = VV_Left[i];
			else
				VV(j,i) = VV_Right[i];
	}
	setInitialField(VV);
}
void ProblemCoreFlows::setInitialFieldStepFunction( int nDim, const vector<double> VV_Left, vector<double> VV_Right, double xstep,
		double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide, EntityType typeField)
{
	Mesh M;
	if(nDim==1)
		M=Mesh(xmin,xmax,nx);
	else if(nDim==2)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny);
	else if(nDim==3)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
	else
	{
		_runLogFile->close();
		throw CdmathException("Space dimension nDim should be between 1 and 3");
	}

	M.setGroupAtPlan(xmax,0,_precision,rightSide);
	M.setGroupAtPlan(xmin,0,_precision,leftSide);
	if(nDim>=2){
		M.setGroupAtPlan(ymax,1,_precision,frontSide);
		M.setGroupAtPlan(ymin,1,_precision,backSide);
	}
	if(nDim==3){
		M.setGroupAtPlan(zmax,2,_precision,topSide);
		M.setGroupAtPlan(zmin,2,_precision,bottomSide);
	}
	Vector V_Left(VV_Left.size()), V_Right(VV_Right.size());
	for(int i=0;i<VV_Left.size(); i++){
		V_Left(i)=VV_Left[i];
		V_Right(i)=VV_Right[i];
	}
	setInitialFieldStepFunction(M, V_Left, V_Right, xstep, typeField);
}

void ProblemCoreFlows::setInitialFieldSphericalStepFunction(const Mesh M, const Vector Vin, const Vector Vout, double radius, const Vector Center, EntityType typeField)
{
	if(_FECalculation && typeField!= NODES)
		cout<<"Warning : finite element simulation should have initial field on nodes!!!"<<endl;

	if((Center.size()!=M.getSpaceDimension()) || (Vout.size() != Vin.size()) )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Vout.size() = %d, Vin.size()= %d, Center.size() = %d, M.getSpaceDim = %d \n",Vout.size(),Vin.size(),Center.size(), M.getSpaceDimension());
		throw CdmathException("ProblemCoreFlows::setInitialFieldSphericalStepFunction : Vector size error");
	}
	int nVar=Vout.size();
	int spaceDim=M.getSpaceDimension();
	Field VV("Primitive variables for spherical step function", typeField, M, nVar);

	Vector currentPoint(spaceDim);
	for(int i=0;i<VV.getNumberOfElements();i++)
	{
		currentPoint(0)=VV.getElementComponent(i,0);
		if(spaceDim>1)
		{
			currentPoint(1)=VV.getElementComponent(i,1);
			if(spaceDim>2)
				currentPoint(2)=VV.getElementComponent(i,2);
		}
		if((currentPoint-Center).norm()<radius)
			for(int j=0;j<nVar;j++)
				VV(i,j)=Vin[j];
		else
			for(int j=0;j<nVar;j++)
				VV(i,j)=Vout[j];
	}
	setInitialField(VV);
}

double ProblemCoreFlows::getTime()
{
	return _time;
}
unsigned ProblemCoreFlows::getNbTimeStep()
{
	return _nbTimeStep;
}
double ProblemCoreFlows::getCFL()
{
	return _cfl;
}
double ProblemCoreFlows::getPrecision()
{
	return _precision;
}
Mesh ProblemCoreFlows::getMesh()
{
	return _mesh;
}

void ProblemCoreFlows::setLinearSolver(linearSolver kspType, preconditioner pcType, double maxIts)
{
	_maxPetscIts=maxIts;
	// set linear solver algorithm
	if (kspType==GMRES)
		_ksptype = (char*)&KSPGMRES;
	else if (kspType==CG)
		_ksptype = (char*)&KSPCG;
	else if (kspType==BCGS)
		_ksptype = (char*)&KSPBCGS;
	else {
		PetscPrintf(PETSC_COMM_WORLD,"!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!\n");
		*_runLogFile << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'GMRES', 'CG' or 'BCGS' algorithm is acceptable !!!");
	}
	// set preconditioner
	if (pcType == NOPC)
		_pctype = (char*)&PCNONE;
	else if (pcType ==LU)
		_pctype = (char*)&PCLU;
	else if (pcType == ILU)
		_pctype = (char*)&PCILU;
	else if (pcType ==CHOLESKY)
		_pctype = (char*)&PCCHOLESKY;
	else if (pcType == ICC)
		_pctype = (char*)&PCICC;
	else {
		PetscPrintf(PETSC_COMM_WORLD,"!!! Error : only 'NOPC', 'LU', 'ILU', 'CHOLESKY' or 'ICC' preconditioners are acceptable !!!\n");
		*_runLogFile << "!!! Error : only 'NOPC' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'NOPC' or 'LU' or 'ILU' preconditioners are acceptable !!!" );
	}
}

void ProblemCoreFlows::setNewtonSolver(double precision, int iterations, nonLinearSolver solverName)
{
	_maxNewtonIts=iterations;
	_precision_Newton=precision;
	_nonLinearSolver=solverName;
}


// Description:
// Cette methode lance une execution du ProblemCoreFlows
// Elle peut etre utilisee si le probleme n'est couple a aucun autre.
// (s'il n'a besoin d'aucun champ d'entree).
// Precondition: initialize
// Seule la methode terminate peut etre appelée apres
bool ProblemCoreFlows::run()
{
	if(!_initializedMemory)
	{
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::run() call initialize() first");
	}
	bool stop=false; // Does the Problem want to stop (error) ?
	bool ok; // Is the time interval successfully solved ?
	_isStationary=false;//in case of a second run with a different physics or cfl

	PetscPrintf(PETSC_COMM_WORLD,"Running test case %d\n",_fileName);

	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation
	*_runLogFile<< "Running test case "<< _fileName<<endl;

	// Time step loop
	while(!stop && !_isStationary &&_time<_timeMax && _nbTimeStep<_maxNbOfTimeStep)
	{
		ok=false; // Is the time interval successfully solved ?

		// Guess the next time step length
		_dt=computeTimeStep(stop);
		if (stop){
			PetscPrintf(PETSC_COMM_WORLD,"Failed computing time step %d, time = %.2f, dt= %.2f, stopping calculation",_nbTimeStep,_time,_dt);
			*_runLogFile << "Failed computing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
			break;
		}
		// Loop on the time interval tries
		while (!ok && !stop )
		{
			stop=!initTimeStep(_dt);
			// Prepare the next time step
			if (stop){
				PetscPrintf(PETSC_COMM_WORLD,"Failed initializing time step %d, time = %.2f, dt= %.2f, stopping calculation",_nbTimeStep,_time,_dt);
				*_runLogFile << "Failed initializing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
				break;
			}
			// Solve the next time step
			ok=solveTimeStep();

			if (!ok)   // The resolution failed, try with a new time interval.
			{
				/* if(_dt>_precision){
					cout<<"ComputeTimeStep returned _dt="<<_dt<<endl;
					cout << "Failed solving time step "<<_nbTimeStep<<", time = " << _time <<" _dt= "<<_dt<<", cfl= "<<_cfl<<", trying again with dt/2"<< endl;
					*_runLogFile << "Failed solving time step "<<_nbTimeStep<<", time = " << _time <<" _dt= "<<_dt<<", cfl= "<<_cfl<<", trying again with dt/2"<< endl;
					double dt=_dt/2;//We chose to divide the time step by 2
					abortTimeStep();//Cancel the initTimeStep
					_dt=dt;//new value of time step is previous time step divided by 2 (we do not call computeTimeStep
					//_cfl*=0.5;//If we change the cfl, we must compute the new time step with computeTimeStep
					//_dt=computeTimeStep(stop);
				}
				else{*/
					PetscPrintf(PETSC_COMM_WORLD,"Failed solving time step %d, time = %.2f, dt= %.2f, cfl = %.2f, stopping calculation \n",_nbTimeStep,_time,_dt,_cfl);
					*_runLogFile << "Failed solving time step "<<_nbTimeStep<<", _time = " << _time<<" _dt= "<<_dt<<", cfl= "<<_cfl <<", stopping calculation"<< endl;
					stop=true; // Impossible to solve the next time step, the Problem has given up
					break;
				//}
			}
			else // The resolution was successful, validate and go to the next time step.
			{
				validateTimeStep();
				if ((_nbTimeStep-1)%_freqSave ==0){
					PetscPrintf(PETSC_COMM_WORLD,"Time step = %d, dt = %.2f, time = %.2f, ||Un+1-Un||= %.2f\n\n",_nbTimeStep,_dt,_time,_erreur_rel);
					*_runLogFile << "Time step = "<< _nbTimeStep << ", dt = "<< _dt <<", time = "<<_time << ", ||Un+1-Un||= "<<_erreur_rel<<endl<<endl;
				}
			}
		}
	}
	if(_isStationary){
		PetscPrintf(PETSC_COMM_WORLD,"Stationary state reached\n");
		*_runLogFile << "Stationary state reached" <<endl;
	}
	else if(_time>=_timeMax){
		PetscPrintf(PETSC_COMM_WORLD,"Maximum time %.2f reached\n",_timeMax);
		*_runLogFile<<"Maximum time "<<_timeMax<<" reached"<<endl;
	}
	else if(_nbTimeStep>=_maxNbOfTimeStep){
		PetscPrintf(PETSC_COMM_WORLD,"Maximum number of time steps %d reached\n",_maxNbOfTimeStep);
		*_runLogFile<<"Maximum number of time steps "<<_maxNbOfTimeStep<<" reached"<<endl;
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD,"Error problem wants to stop!\n");
		*_runLogFile<<"Error problem wants to stop!"<<endl;
	}
	PetscPrintf(PETSC_COMM_WORLD,"End of calculation at time t = %.2f and time step number %d\n",_time,_nbTimeStep);
	*_runLogFile << "End of calculation time t= " << _time << " at time step number "<< _nbTimeStep << endl;

	_runLogFile->close();
	return !stop;
}

void ProblemCoreFlows::displayMatrix(double *matrix, int size, string name)
{
	cout<<name<<endl;
	for(int p=0; p<size; p++)
	{
		for(int q=0; q<size; q++)
		{
			cout << matrix[p*size+q] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void ProblemCoreFlows::displayVector(double *vector, int size, string name)
{
	cout<<name<<endl;
	for(int q=0; q<size; q++)
	{
		cout << vector[q] << "\t";
	}
	cout << endl;
}
void ProblemCoreFlows::setFileName(string fileName){
	if( _nbTimeStep>0 && fileName!=_fileName)//This is a change of file name during a simulation
		_restartWithNewFileName=true;
	_fileName = fileName;
}

void ProblemCoreFlows::setFreqSave(int freqSave){
	_freqSave = freqSave;
}

bool ProblemCoreFlows::solveTimeStep(){
	_NEWTON_its=0;
	bool converged=false, ok=true;
	while(!converged && ok && _NEWTON_its < _maxNewtonIts){
		ok=iterateTimeStep(converged);//resolution du systeme lineaire si schema implicite

		if(_timeScheme == Implicit && (_nbTimeStep-1)%_freqSave ==0)//To monitor the convergence of the newton scheme
		{
			PetscPrintf(PETSC_COMM_WORLD," Newton iteration %d, %s iterations : %d maximum variation ||Uk+1-Uk||: %.2f\n",_NEWTON_its,_ksptype,_PetscIts,_erreur_rel);
			*_runLogFile<< " Newton iteration " << _NEWTON_its<< ", "<< _ksptype << " iterations : " << _PetscIts<< " maximum variation ||Uk+1-Uk||: " << _erreur_rel << endl;

			if(_conditionNumber)
			{
				PetscReal sv_max, sv_min;
				KSPComputeExtremeSingularValues(_ksp, &sv_max, &sv_min);
				PetscPrintf(PETSC_COMM_WORLD," Singular value max = %.2f, singular value min = %.2f, condition number = %.2f\n",sv_max,sv_min,sv_max/sv_min);
				*_runLogFile<<" Singular value max = " << sv_max <<", singular value min = " << sv_min <<", condition number = " << sv_max/sv_min <<endl;
			}
		}
		_NEWTON_its++;
	}
	if(!converged){
		if(_NEWTON_its >= _maxNewtonIts){
			PetscPrintf(PETSC_COMM_WORLD,"Maximum number of Newton iterations %d reached\n",_maxNewtonIts);
			*_runLogFile << "Maximum number of Newton iterations "<<_maxNewtonIts<<" reached"<< endl;
		}
		else if(!ok){
			PetscPrintf(PETSC_COMM_WORLD,"iterateTimeStep: solving Newton iteration %d Failed\n",_NEWTON_its);
			*_runLogFile<<"iterateTimeStep: solving Newton iteration "<<_NEWTON_its<<" Failed"<<endl;
		}
	}
	else if(_timeScheme == Implicit && (_nbTimeStep-1)%_freqSave ==0)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Nombre d'iterations de Newton %d, Nombre max d'iterations %s : %d\n\n",_NEWTON_its, _ksptype, _MaxIterLinearSolver);
		*_runLogFile <<endl;
		*_runLogFile << "Nombre d'iterations de Newton "<< _NEWTON_its << "Nombre max d'iterations "<< _ksptype << " : " << _MaxIterLinearSolver << endl << endl;
		_MaxIterLinearSolver = 0;
	}

	return converged;
}

ProblemCoreFlows::~ProblemCoreFlows()
{
	/*
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(petscInitialized)
		PetscFinalize();
	 */
	delete _runLogFile;
}

double 
ProblemCoreFlows::getConditionNumber(bool isSingular, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getConditionNumber( isSingular, tol);
}
std::vector< double > 
ProblemCoreFlows::getEigenvalues(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getEigenvalues( nev, which, tol);
}
std::vector< Vector > 
ProblemCoreFlows::getEigenvectors(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getEigenvectors( nev, which, tol);
}
Field 
ProblemCoreFlows::getEigenvectorsField(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  MEDCoupling::DataArrayDouble * d = A.getEigenvectorsDataArrayDouble( nev, which, tol);
  Field my_eigenfield;
  
  my_eigenfield = Field("Eigenvectors field", _VV.getTypeOfField(), _mesh, nev);

  my_eigenfield.setFieldByDataArrayDouble(d);
  
  return my_eigenfield;
}

std::vector< double > 
ProblemCoreFlows::getSingularValues( int nsv, SVDWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getSingularValues( nsv, which, tol);
}
std::vector< Vector > 
ProblemCoreFlows::getSingularVectors(int nsv, SVDWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getSingularVectors( nsv, which, tol);
}

Field 
ProblemCoreFlows::getUnknownField() const
{
	if(!_initializedMemory)
	{
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::getUnknownField() call initialize() first");
	}
	return _VV;
}
