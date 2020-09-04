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

ProblemCoreFlows::ProblemCoreFlows()
{
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
		PetscInitialize(NULL,NULL,0,0);
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
	_Ndim=0;
	_minl=0;
	_neibMaxNb=0;
	_fileName = "myCoreFlowsProblem";
	_freqSave = 1;
	_initialDataSet=false;
	_initializedMemory=false;
	_restartWithNewTimeScheme=false;
	_restartWithNewFileName=false;
	_timeScheme=Explicit;
	_wellBalancedCorrection=false;
    _FECalculation=false;
	_maxPetscIts=50;
	_MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
	_maxNewtonIts=50;
	_NEWTON_its=0;
	int _PetscIts=0;//the number of iterations of the linear solver
	_ksptype = (char*)&KSPGMRES;
	_pctype = (char*)&PCLU;
	_heatPowerFieldSet=false;
	_heatTransfertCoeff=0;
	_rodTemperatureFieldSet=false;
	_heatSource=0;
	_verbose = false;
	_system = false;
	_conditionNumber=false;
	_maxvp=0;
	_runLogFile=new ofstream;

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

	_VV=VV;
	_time=_VV.getTime();
	_mesh=_VV.getMesh();
	_Nmailles = _mesh.getNumberOfCells();
	_Nnodes =   _mesh.getNumberOfNodes();
	_Nfaces =   _mesh.getNumberOfFaces();
	_perimeters=Field("Perimeters", CELLS, _mesh,1);

	// find _minl and maximum nb of neibourghs
	_minl  = INFINITY;
	int nbNeib,indexFace;
	Cell Ci;
	Face Fk;

	if(_verbose)
		cout<<"Computing cell perimeters and mesh minimal diameter"<<endl;

    if(VV.getTypeOfField()==NODES)
    {
        _minl = _mesh.getMaxNbNeighbours(NODES);
        _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
    }
    else
        for (int i=0; i<_mesh.getNumberOfCells(); i++){
            Ci = _mesh.getCell(i);
            //Detect mesh with junction
            nbNeib=0;
            for(int j=0; j<Ci.getNumberOfFaces(); j++){
                Fk=_mesh.getFace(Ci.getFacesId()[j]);
                nbNeib+=Fk.getNumberOfCells()-1;
            }
            if(nbNeib>_neibMaxNb)
                _neibMaxNb=nbNeib;
            //Compute mesh data
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
	_initialDataSet=true;

	if(_verbose)
		cout<<_perimeters<<endl;
}
void ProblemCoreFlows::setInitialField(string fileName, string fieldName, int timeStepNumber)
{
	Field VV(fileName,CELLS,fieldName,timeStepNumber,0);
	setInitialField(VV);
}
void ProblemCoreFlows::setInitialFieldConstant(string fileName, const vector<double> Vconstant)
{
	Mesh M(fileName);
	Field VV("Primitive", CELLS, M, Vconstant.size());

	for (int j = 0; j < M.getNumberOfCells(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant[i];
	}

	setInitialField(VV);
}
void ProblemCoreFlows::	setInitialFieldConstant(const Mesh& M, const Vector Vconstant)
{
	Field VV("Primitive", CELLS, M, Vconstant.getNumberOfRows());

	for (int j = 0; j < M.getNumberOfCells(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant(i);
	}
	setInitialField(VV);
}
void ProblemCoreFlows::	setInitialFieldConstant(const Mesh& M, const vector<double> Vconstant)
{
	Field VV("Primitive", CELLS, M, Vconstant.size());

	for (int j = 0; j < M.getNumberOfCells(); j++) {
		for (int i=0; i< VV.getNumberOfComponents(); i++)
			VV(j,i) = Vconstant[i];
	}
	setInitialField(VV);
}
void ProblemCoreFlows::setInitialFieldConstant( int nDim, const vector<double> Vconstant, double xmin, double xmax, int nx, string leftSide, string rightSide,
		double ymin, double ymax, int ny, string backSide, string frontSide,
		double zmin, double zmax, int nz, string bottomSide, string topSide)
{
	Mesh M;
	if(nDim==1){
		//cout<<"coucou1 xmin="<<xmin<<", xmax= "<<xmax<< ", nx= "<<nx<<endl;
		M=Mesh(xmin,xmax,nx);
		//cout<<"coucou2"<<endl;
	}
	else if(nDim==2)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny);
	else if(nDim==3)
		M=Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
	else{
		cout<<"ProblemCoreFlows::setInitialFieldConstant: Space dimension nDim should be between 1 and 3"<<endl;
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

	setInitialFieldConstant(M, Vconstant);
}
void ProblemCoreFlows::setInitialFieldStepFunction(const Mesh M, const Vector VV_Left, const Vector VV_Right, double disc_pos, int direction)
{
	if  (VV_Right.getNumberOfRows()!=VV_Left.getNumberOfRows())
	{
		*_runLogFile<<"ProblemCoreFlows::setStepFunctionInitialField: Vectors VV_Left and VV_Right have different sizes"<<endl;
		_runLogFile->close();
		throw CdmathException( "ProblemCoreFlows::setStepFunctionInitialField: Vectors VV_Left and VV_Right have different sizes");
	}
	Field VV("Primitive", CELLS, M, VV_Left.getNumberOfRows());

	double component_value;

	for (int j = 0; j < M.getNumberOfCells(); j++) {
		if(direction==0)
			component_value=M.getCell(j).x();
		else if(direction==1)
			component_value=M.getCell(j).y();
		else if(direction==2)
			component_value=M.getCell(j).z();
		else{
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
		double zmin, double zmax, int nz, string bottomSide, string topSide)
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
	setInitialFieldStepFunction(M, V_Left, V_Right, xstep);
}

void ProblemCoreFlows::setInitialFieldSphericalStepFunction(const Mesh M, const Vector Vin, const Vector Vout, double radius, const Vector Center)
{
	if((Center.size()!=M.getSpaceDimension()) || (Vout.size() != Vin.size()) )
	{
		cout<< "Vout.size()= "<<Vout.size() << ", Vin.size()= "<<Vin.size()<<", Center.size()="<<Center.size()<<", M.getSpaceDim= "<< M.getSpaceDimension()<<endl;
		throw CdmathException("ProblemCoreFlows::setInitialFieldSphericalStepFunction : Vector size error");
	}
	int nVar=Vout.size();
	int spaceDim=M.getSpaceDimension();
	Field VV("Primitive variables for spherical step function", CELLS, M, nVar);

	Vector currentPoint(spaceDim);
	for(int i=0;i<M.getNumberOfCells();i++)
	{
		currentPoint(0)=M.getCell(i).x();
		if(spaceDim>1)
		{
			currentPoint(1)=M.getCell(i).y();
			if(spaceDim>2)
				currentPoint(2)=M.getCell(i).z();
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

void ProblemCoreFlows::setLinearSolver(linearSolver kspType, preconditioner pcType)
{
	//_maxPetscIts=maxIterationsPetsc;
	// set linear solver algorithm
	if (kspType==GMRES)
		_ksptype = (char*)&KSPGMRES;
	else if (kspType==CG)
		_ksptype = (char*)&KSPCG;
	else if (kspType==BCGS)
		_ksptype = (char*)&KSPBCGS;
	else {
		cout << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
		*_runLogFile << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'GMRES', 'CG' or 'BCGS' algorithm is acceptable !!!");
	}
	// set preconditioner
	if (pcType == NONE)
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
		cout << "!!! Error : only 'NONE', 'LU', 'ILU', 'CHOLESKY' or 'ICC' preconditioners are acceptable !!!" << endl;
		*_runLogFile << "!!! Error : only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" );
	}
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

	cout<< "Running test case "<< _fileName<<endl;

	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation
	*_runLogFile<< "Running test case "<< _fileName<<endl;

	// Time step loop
	while(!stop && !_isStationary &&_time<_timeMax && _nbTimeStep<_maxNbOfTimeStep)
	{
		ok=false; // Is the time interval successfully solved ?

		// Guess the next time step length
		_dt=computeTimeStep(stop);
		if (stop){
			cout << "Failed computing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
			*_runLogFile << "Failed computing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
			break;
		}
		// Loop on the time interval tries
		while (!ok && !stop )
		{
			stop=!initTimeStep(_dt);
			// Prepare the next time step
			if (stop){
				cout << "Failed initializing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
				*_runLogFile << "Failed initializing time step "<<_nbTimeStep<<", time = " << _time <<", dt= "<<_dt<<", stopping calculation"<< endl;
				break;
			}
			// Solve the next time step
			ok=solveTimeStep();

			if (!ok)   // The resolution failed, try with a new time interval.
			{
				if(_dt>_precision){
					cout<<"ComputeTimeStep returned _dt="<<_dt<<endl;
					cout << "Failed solving time step "<<_nbTimeStep<<", time = " << _time <<" _dt= "<<_dt<<", cfl= "<<_cfl<<", trying again with dt/2"<< endl;
					*_runLogFile << "Failed solving time step "<<_nbTimeStep<<", time = " << _time <<" _dt= "<<_dt<<", cfl= "<<_cfl<<", trying again with dt/2"<< endl;
					double dt=_dt/2;//We chose to divide the time step by 2
					abortTimeStep();//Cancel the initTimeStep
					_dt=dt;//new value of time step is previous time step divided by 2 (we do not call computeTimeStep
					//_cfl*=0.5;//If we change the cfl, we must compute the new time step with computeTimeStep
					//_dt=computeTimeStep(stop);
				}
				else{
					cout << "Failed solving time step "<<_nbTimeStep<<", _time = " << _time<<" _dt= "<<_dt<<", cfl= "<<_cfl <<", stopping calculation"<< endl;
					*_runLogFile << "Failed solving time step "<<_nbTimeStep<<", _time = " << _time<<" _dt= "<<_dt<<", cfl= "<<_cfl <<", stopping calculation"<< endl;
					stop=true; // Impossible to solve the next time step, the Problem has given up
					break;
				}
			}
			else // The resolution was successful, validate and go to the next time step.
			{
				validateTimeStep();
				if (_nbTimeStep%_freqSave ==0){
					cout << "Time step = "<< _nbTimeStep << ", dt = "<< _dt <<", time = "<<_time << ", ||Un+1-Un||= "<<_erreur_rel<<endl;
					*_runLogFile << "Time step = "<< _nbTimeStep << ", dt = "<< _dt <<", time = "<<_time << ", ||Un+1-Un||= "<<_erreur_rel<<endl;
				}
			}
		}
	}
	if(_isStationary){
		cout << "Stationary state reached" <<endl;
		*_runLogFile << "Stationary state reached" <<endl;
	}
	else if(_time>=_timeMax){
		cout<<"Maximum time "<<_timeMax<<" reached"<<endl;
		*_runLogFile<<"Maximum time "<<_timeMax<<" reached"<<endl;
	}
	else if(_nbTimeStep>=_maxNbOfTimeStep){
		cout<<"Maximum number of time steps "<<_maxNbOfTimeStep<<" reached"<<endl;
		*_runLogFile<<"Maximum number of time steps "<<_maxNbOfTimeStep<<" reached"<<endl;
	}
	else{
		cout<<"Error problem wants to stop!"<<endl;
		*_runLogFile<<"Error problem wants to stop!"<<endl;
	}
	cout << "End of calculation time t= " << _time << " at time step number "<< _nbTimeStep << endl;
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

		if(_timeScheme == Implicit && _nbTimeStep%_freqSave ==0)//To monitor the convergence of the newton scheme
		{
			cout << "\n Newton iteration " << _NEWTON_its<< ", "<< _ksptype << " iterations : " << " : " << _PetscIts<< " maximum variation ||Uk+1-Uk||: " << _erreur_rel << endl;
			*_runLogFile<< "\n Newton iteration " << _NEWTON_its<< ", "<< _ksptype << " iterations : " << " : " << _PetscIts<< " maximum variation ||Uk+1-Uk||: " << _erreur_rel << endl;

			if(_conditionNumber)
			{
				PetscReal sv_max, sv_min;
				KSPComputeExtremeSingularValues(_ksp, &sv_max, &sv_min);
				cout<<" singular value max = " << sv_max <<" singular value min = " << sv_min <<" condition number = " << sv_max/sv_min <<endl;
				*_runLogFile<<" singular value max = " << sv_max <<" singular value min = " << sv_min <<" condition number = " << sv_max/sv_min <<endl;
			}
		}
		_NEWTON_its++;
	}
	if(!converged){
		if(_NEWTON_its >= _maxNewtonIts){
			cout << "Maximum number of Newton iterations "<<_maxNewtonIts<<" reached"<< endl;
			*_runLogFile << "Maximum number of Newton iterations "<<_maxNewtonIts<<" reached"<< endl;
		}
		else if(!ok){
			cout<<"iterateTimeStep: solving Newton iteration "<<_NEWTON_its<<" Failed"<<endl;
			*_runLogFile<<"iterateTimeStep: solving Newton iteration "<<_NEWTON_its<<" Failed"<<endl;
		}
	}
	else if(_timeScheme == Implicit && _nbTimeStep%_freqSave ==0)
	{
		cout<<endl;
		cout << "Nombre d'iterations de Newton "<< _NEWTON_its << ", Nombre max d'iterations "<< _ksptype << " : " << _MaxIterLinearSolver << endl;
		*_runLogFile <<endl;
		*_runLogFile << "Nombre d'iterations de Newton "<< _NEWTON_its << " Variation ||Un+1-Un||= "<< _erreur_rel<<endl;
		*_runLogFile << "Nombre max d'iterations "<< _ksptype << " : " << _MaxIterLinearSolver << endl;
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
  
  if(_FECalculation)
    my_eigenfield = Field("Eigenvectors field", NODES, _mesh, nev);
  else
    my_eigenfield = Field("Eigenvectors field", CELLS, _mesh, nev);

  my_eigenfield.setFieldByDataArrayDouble(d);
  
  return my_eigenfield;
}
