#include "TransportEquation.hxx"
#include "math.h"
#include <fstream>
#include <sstream>

using namespace std;

TransportEquation::TransportEquation(phase fluid, pressureMagnitude pEstimate,vector<double> vitesseTransport, MPI_Comm comm):ProblemCoreFlows(comm)
{
	if(pEstimate==around1bar300KTransport){
		_Tref=300;
		if(fluid==GasPhase){//Nitrogen pressure 1 bar and temperature 27°C
			_href=3.11e5; //nitrogen enthalpy at 1 bar and 300K
			_cpref=1041;//nitrogen specific heat at constant pressure 1 bar and 300K
			//saturation data for nitrogen at 1 bar and 77K
			_hsatv=0.77e5;//nitrogen vapour enthalpy at saturation at 1 bar
			_hsatl=-1.22e5;//nitrogen liquid enthalpy at saturation at 1 bar
			_rhosatv=4.556;//nitrogen vapour density at saturation at 1 bar
			_rhosatl=806.6;//nitrogen liquid density at saturation at 1 bar
		}
		else{
			//Water at pressure 1 bar and temperature 27°C
			_href=1.127e5; //water enthalpy at 1 bar and 300K
			_cpref=4181;//water specific heat at 1 bar and 300K
			//saturation data for water at 1 bar and 373K
			_hsatv=2.675e6;//Gas enthalpy at saturation at 1 bar
			_hsatl=4.175e5;//water enthalpy at saturation at 1 bar
			_rhosatv=0.6;//Gas density at saturation at 1 bar
			_rhosatl=958;//water density at saturation at 1 bar
		}
	}
	else{//around155bars600K
		_Tref=618;//=Tsat
		if(fluid==GasPhase){
			_href=2.675e6; //Gas enthalpy at 155 bars and 618K
			_cpref=14001;//Gas specific heat at 155 bar and 618K
		}
		else{//Liquid
			_href=4.175e5;//water enthalpy at 155 bars and 618K
			_cpref=8950;//water specific heat at 155 bar and 618K
		}
		//saturation data for water at 155 bars and 618K
		_hsatv=2.6e6;//Gas enthalpy at saturation at 155 bars
		_hsatl=1.63e6;//water enthalpy at saturation at 155 bars
		_rhosatv=101.9;//Gas density at saturation at 155 bars
		_rhosatl=594.4;//water density at saturation at 155 bars
	}
	_Ndim=vitesseTransport.size();
	_vitesseTransport=Vector(_Ndim);
	for(int i=0;i<_Ndim;i++)
		_vitesseTransport[i]=vitesseTransport[i];
	_nVar=1;
	_dt_transport=0;
	_dt_src=0;
	_transportMatrixSet=false;
	_FECalculation=false;//Only finite volumes available
	_rodTemperatureFieldSet=false;
	_rodTemperature=0;
}

void TransportEquation::initialize()
{
	if(_mpi_rank==0)
	{
		if(!_initialDataSet)
			throw CdmathException("TransportEquation::initialize() set initial data first");
		else if (_VV.getTypeOfField() != CELLS)
			throw CdmathException("TransportEquation::initialize() Initial data should be a field on CELLS, not NODES, neither FACES");
		else
			PetscPrintf(PETSC_COMM_SELF,"Initialising the transport of a fluid enthalpy\n");
	
		/**************** Field creation *********************/
	
		//post processing fields used only for saving results
		_TT=Field ("Temperature", CELLS, _mesh, 1);
		_Alpha=Field ("Void fraction", CELLS, _mesh, 1);
		_Rho=Field ("Mixture density", CELLS, _mesh, 1);
		//Construction des champs de post-traitement
		for(int i =0; i<_Nmailles;i++){
			_TT(i)=temperature(_VV(i));
			_Alpha(i)=voidFraction(_VV(i));
			_Rho(i)=density(_Alpha(i));
		}
		if(!_heatPowerFieldSet){
			_heatPowerField=Field("Heat Power",CELLS,_mesh,1);
			for(int i =0; i<_Nmailles; i++)
				_heatPowerField(i) = _heatSource;
		}
		if(!_rodTemperatureFieldSet){
			_rodTemperatureField=Field("Rod temperature",CELLS,_mesh,1);
			for(int i =0; i<_Nmailles; i++)
				_rodTemperatureField(i) = _rodTemperature;
		}
	}
	
	_globalNbUnknowns = _Nmailles;
	
	/* Vectors creations */
	VecCreate(PETSC_COMM_WORLD, &_Hn);
	VecSetSizes(_Hn,PETSC_DECIDE,_Nmailles);
	VecSetFromOptions(_Hn);
	VecGetLocalSize(_Hn, &_localNbUnknowns);

	VecDuplicate(_Hn, &_Hk);
	VecDuplicate(_Hn, &_Hkm1);
	VecDuplicate(_Hn, &_deltaH);
	VecDuplicate(_Hn, &_b);//RHS of the linear system: _b=Hn/dt + _b0 + puisance
	VecDuplicate(_Hn, &_b0);//part of the RHS that comes from the boundary conditions

	if(_mpi_rank == 0)//Process 0 reads and distributes initial data
		for(int i =0; i<_Nmailles;i++)
			VecSetValue(_Hn,i,_VV(i), INSERT_VALUES);
	VecAssemblyBegin(_Hn);
	VecAssemblyEnd(_Hn);

	//creation de la matrice
   	MatCreateAIJ(PETSC_COMM_WORLD, _localNbUnknowns, _localNbUnknowns, _globalNbUnknowns, _globalNbUnknowns, _d_nnz, PETSC_NULL, _o_nnz, PETSC_NULL, &_A);

	/* Local sequential vector creation */
	if(_mpi_size>1 && _mpi_rank == 0)
		VecCreateSeq(PETSC_COMM_SELF,_globalNbUnknowns,&_Hn_seq);//For saving results on proc 0
	VecScatterCreateToZero(_Hn,&_scat,&_Hn_seq);

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

	_initializedMemory=true;
	save();//save initial data
}

double TransportEquation::computeTimeStep(bool & stop){
	if(!_transportMatrixSet)
		_dt_transport=computeTransportMatrix();

	_dt_src=computeRHS();

    if(_verbose or _system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Right hand side of the linear system\n");
        VecView(_b,PETSC_VIEWER_STDOUT_WORLD);
	}

	stop=false;
	return min(_dt_transport,_dt_src);
}
double TransportEquation::computeTransportMatrix(){
	MatZeroEntries(_A);
	VecZeroEntries(_b0);

	if(_mpi_rank == 0)
	{
		long nbFaces = _mesh.getNumberOfFaces();
		Face Fj;
		Cell Cell1,Cell2;
		string nameOfGroup;
		double inv_dxi, inv_dxj;
		Vector normale(_Ndim);
		double un, hk;
		PetscInt idm, idn;
		std::vector< int > idCells;
		for (int j=0; j<nbFaces;j++){
			Fj = _mesh.getFace(j);
	
			// compute the normal vector corresponding to face j : from idCells[0] to idCells[1]
			idCells = Fj.getCellsId();
			Cell1 = _mesh.getCell(idCells[0]);
			idm = idCells[0];
			if (_Ndim >1){
				for(int l=0; l<Cell1.getNumberOfFaces(); l++){
					if (j == Cell1.getFacesId()[l]){
						for (int idim = 0; idim < _Ndim; ++idim)
							normale[idim] = Cell1.getNormalVector(l,idim);
						break;
					}
				}
			}else{ // _Ndim = 1 : assume that this is normal mesh : the face index increases in positive direction
				if (Fj.getNumberOfCells()<2) {
					if (j==0)
						normale[0] = -1;
					else if (j==nbFaces-1)
						normale[0] = 1;
					else
						throw CdmathException("TransportEquation::ComputeTimeStep(): computation of normal vector failed");
				} else if(Fj.getNumberOfCells()==2){
					if (idCells[0] < idCells[1])
						normale[0] = 1;
					else
						normale[0] = -1;
				}
			}
			//Compute velocity at the face Fj
			un=normale*_vitesseTransport;
			if(abs(un)>_maxvp)
				_maxvp=abs(un);
	
			// compute 1/dxi = volume of Ci/area of Fj
			if (_Ndim > 1)
				inv_dxi = Fj.getMeasure()/Cell1.getMeasure();
			else
				inv_dxi = 1/Cell1.getMeasure();
	
			// If Fj is on the boundary
			if (Fj.getNumberOfCells()==1) {
				if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
				{
					cout << "face numero " << j << " cellule frontiere " << idCells[0] << " ; vecteur normal=(";
					for(int p=0; p<_Ndim; p++)
						cout << normale[p] << ",";
					cout << ") "<<endl;
				}
				nameOfGroup = Fj.getGroupName();
	
				if     (_limitField[nameOfGroup].bcType==NeumannTransport || _limitField[nameOfGroup].bcType==OutletTransport ){
					MatSetValue(_A,idm,idm,inv_dxi*un, ADD_VALUES);
				}
				else if(_limitField[nameOfGroup].bcType==InletTransport   || _limitField[nameOfGroup].bcType==DirichletTransport){
					if(un>0){
						MatSetValue(_A,idm,idm,inv_dxi*un, ADD_VALUES);
					}
					else{
						hk=_limitField[nameOfGroup].h;
						VecSetValue(_b0,idm,-inv_dxi*un*hk, ADD_VALUES);
					}
				}
				else {//bcType=NoneBCTransport
					cout<<"!!!!!!!!!!!!!!! Error TransportEquation::computeTransportMatrix() !!!!!!!!!!"<<endl;
					cout<<"!!!!!!!!! Boundary condition not set for boundary named "<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<" !!!!!!!!!!!!!! "<<endl;
					cout<<"Accepted boundary conditions are NeumannTransport "<<NeumannTransport<< " and InletTransport "<< InletTransport <<endl;
					throw CdmathException("Boundary condition not accepted");
				}
				// if Fj is inside the domain
			} else 	if (Fj.getNumberOfCells()==2 ){
				if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
				{
					cout << "face numero " << j << " cellule gauche " << idCells[0] << " cellule droite " << idCells[1];
					cout << " ; vecteur normal=(";
					for(int p=0; p<_Ndim; p++)
						cout << normale[p] << ",";
					cout << "). "<<endl;
				}
				Cell2 = _mesh.getCell(idCells[1]);
				idn = idCells[1];
				if (_Ndim > 1)
					inv_dxj = Fj.getMeasure()/Cell2.getMeasure();
				else
					inv_dxj = 1/Cell2.getMeasure();
	
				if(un>0){
					MatSetValue(_A,idm,idm,inv_dxi*un, ADD_VALUES);
					MatSetValue(_A,idn,idm,-inv_dxj*un, ADD_VALUES);
				}
				else{
					MatSetValue(_A,idm,idn,inv_dxi*un, ADD_VALUES);
					MatSetValue(_A,idn,idn,-inv_dxj*un, ADD_VALUES);
				}
			}
			else
				throw CdmathException("TransportEquation::ComputeTimeStep(): incompatible number of cells around a face");
		}
	}

	_transportMatrixSet=true;

	MPI_Bcast(&_maxvp, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Maximum speed is %.2f, CFL = %.2f, Delta x = %.2f\n",_maxvp,_cfl,_minl);

    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);          
	VecAssemblyEnd(  _b0);
	
	if(abs(_maxvp)<_precision)
		throw CdmathException("TransportEquation::computeTransportMatrix(): maximum eigenvalue for time step is zero");
	else
		return _cfl*_minl/_maxvp;
}
double TransportEquation::computeRHS(){
	double rhomin=INFINITY;

	VecCopy(_b0,_b);

	if(_mpi_rank == 0)
	{
		if(_system)
			cout<<"Second membre of transport problem"<<endl;
		for (int i=0; i<_Nmailles;i++){
			VecSetValue(_b,i,_heatTransfertCoeff*(_rodTemperatureField(i)-_TT(i))/_Rho(i),ADD_VALUES);
			VecSetValue(_b,i,_heatPowerField(i)/_Rho(i),ADD_VALUES);
			if(_system)
				cout<<_heatPowerField(i)/_Rho(i)<<endl;
			if(_Rho(i)<rhomin)
				rhomin=_Rho(i);
		}
	}
	VecAssemblyBegin(_b);
	VecAssemblyEnd(  _b);

    if(_verbose or _system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Right hand side of the linear system\n");
        VecView(_b,PETSC_VIEWER_STDOUT_WORLD);
	}

	return rhomin*_cpref/_heatTransfertCoeff;
}
void TransportEquation::updatePrimitives()
{
	if(_mpi_size>1){
		VecScatterBegin(_scat,_Hn,_Hn_seq,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(  _scat,_Hn,_Hn_seq,INSERT_VALUES,SCATTER_FORWARD);
	}
	
    if(_verbose or _system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Unknown of the linear system :\n");
        VecView(_Hn,PETSC_VIEWER_STDOUT_WORLD);
	}

	if(_mpi_rank == 0)
	{
		double hi;
		for(int i=0; i<_Nmailles; i++)
		{
			if(_mpi_size>1)
				VecGetValues(_Hn_seq, 1, &i, &hi);
			else
				VecGetValues(_Hn    , 1, &i, &hi);
			_VV(i)=hi;
			_TT(i)=temperature(hi);
			_Alpha(i)=voidFraction(hi);
			_Rho(i)=density(_Alpha(i));
		}
	}
}

bool TransportEquation::initTimeStep(double dt){
	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Matrix of the linear system\n");
		MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
	}

    if(_dt>0 and dt>0)
    {
        //Remove the contribution from dt to prepare for new time step. The diffusion matrix is not recomputed
        if(_timeScheme == Implicit)
            MatShift(_A,-1/_dt+1/dt);
        //No need to remove the contribution to the right hand side since it is recomputed from scratch at each time step
    }
    else if(dt>0)//_dt==0, first time step
    {
		if(_timeScheme == Implicit)
			MatShift(_A,1/dt);
	}
    else//dt<=0
    {
        PetscPrintf(PETSC_COMM_WORLD,"TransportEquation::initTimeStep %.2f = \n",dt);
        throw CdmathException("Error TransportEquation::initTimeStep : cannot set time step to zero");        
    }
    //At this stage _b contains _b0 + power + heat exchange
    VecAXPY(_b, 1/dt, _Hn);        
	
	_dt=dt;
	
	return _dt>0;
}

void TransportEquation::abortTimeStep(){
    //Remove contribution of dt to the RHS
	VecAXPY(_b,  -1/_dt, _Hn);

    //Remove contribution of dt to the matrix
	if(_timeScheme == Implicit)
		MatShift(_A,-1/_dt);

	_dt = 0;
}

bool TransportEquation::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		_maxvp=0;
		computeTimeStep(stop);
	}
	if(stop){
		converged=false;
		return false;
	}
	VecAXPY(_b, 1/_dt, _Hn);
	if(_system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Vecteur Hn : \n");
		VecView(_Hn,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
		PetscPrintf(PETSC_COMM_WORLD,"Vecteur _b : \n");
		VecView(_b,  PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"Matrice A : \n");
		MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
	}

	if(_timeScheme == Explicit)
	{
		MatMult(_A, _Hn, _Hk);
		if(_system)
		{
			PetscPrintf(PETSC_COMM_WORLD,"Nouveau vecteur Hk: \n");
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		VecAXPY(_Hk, -1, _b);
		if(_system)
		{
			PetscPrintf(PETSC_COMM_WORLD,"Nouveau vecteur Hk-b: \n");
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		VecScale(_Hk, -_dt);
		if(_system)
		{
			PetscPrintf(PETSC_COMM_WORLD,"Nouveau vecteur dt*(Hk-b): \n");
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		converged = true;
	}
	else
	{

#if PETSC_VERSION_GREATER_3_5
		KSPSetOperators(_ksp, _A, _A);
#else
		KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

		if(_conditionNumber)
			KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);

		KSPSolve(_ksp, _b, _Hk);
		MatShift(_A,-1/_dt);

		KSPGetIterationNumber(_ksp, &_PetscIts);
		if( _MaxIterLinearSolver < _PetscIts)
			_MaxIterLinearSolver = _PetscIts;
		if(_PetscIts>=_maxPetscIts)
		{
			PetscPrintf(PETSC_COMM_WORLD,"Systeme lineaire : pas de convergence de Petsc. Itérations maximales %d atteintes", _maxPetscIts);
			converged=false;
			return false;
		}
		else
		{
			VecCopy(_Hk, _deltaH);//ici on a deltaH=Hk
			VecAXPY(_deltaH,  -1, _Hkm1);//On obtient deltaH=Hk-Hkm1
			VecNorm(_deltaH,NORM_INFINITY,&_erreur_rel);
			converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
		}
	}

	VecCopy(_Hk, _Hkm1);

	return true;
}
void TransportEquation::validateTimeStep()
{
	VecCopy(_Hk, _deltaH);//ici Hk=Hnp1 donc on a deltaH=Hnp1
	VecAXPY(_deltaH,  -1, _Hn);//On obtient deltaH=Hnp1-Hn
	VecNorm(_deltaH,NORM_INFINITY,&_erreur_rel);

	_isStationary =(_erreur_rel <_precision);

	VecCopy(_Hk, _Hn);
	VecCopy(_Hk, _Hkm1);

	updatePrimitives();

	if(_mpi_rank == 0)
		if((_nbTimeStep-1)%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
		{
			cout <<"Valeur propre locale max: " << _maxvp << endl;
			//Find minimum and maximum void fractions
			double alphamin=INFINITY;
			double alphamax=-INFINITY;
			for(int i=0; i<_Nmailles; i++)
			{
				if(_Alpha(i)>alphamax)
					alphamax=_Alpha(i);
				if(_Alpha(i)<alphamin)
					alphamin=_Alpha(i);
			}
			cout<<"Alpha min = " << alphamin << " Alpha max = " << alphamax<<endl;
		}

	_time+=_dt;
	_nbTimeStep++;
	save();
}

void TransportEquation::terminate(){
	VecDestroy(&_Hn);
	VecDestroy(&_Hk);
	VecDestroy(&_Hkm1);
	VecDestroy(&_deltaH);
	VecDestroy(&_b0);
	VecDestroy(&_b);
	MatDestroy(&_A);
	if(_mpi_size>1 && _mpi_rank == 0)
		VecDestroy(&_Hn_seq);
}

void TransportEquation::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results\n\n");

	string resultFile(_path+"/TransportEquation_");///Results
	resultFile+=_fileName;

	if(_mpi_rank==0){
		_VV.setTime(_time,_nbTimeStep);
		_TT.setTime(_time,_nbTimeStep);
		_Alpha.setTime(_time,_nbTimeStep);
		_Rho.setTime(_time,_nbTimeStep);
	
		// create mesh and component info
		if (_nbTimeStep ==0 || _restartWithNewFileName){
			if (_restartWithNewFileName)
				_restartWithNewFileName=false;
			string suppress ="rm -rf "+resultFile+"_*";
			system(suppress.c_str());//Nettoyage des précédents calculs identiques
	
			switch(_saveFormat)
			{
			case VTK :
				_VV.writeVTK(resultFile+"Enthalpy");
				_TT.writeVTK(resultFile+"Temperature");
				_Alpha.writeVTK(resultFile+"GasFraction");
				_Rho.writeVTK(resultFile+"MixtureDensity");
				break;
			case MED :
				_VV.writeMED(resultFile+"Enthalpy");
				_TT.writeMED(resultFile+"Temperature");
				_Alpha.writeMED(resultFile+"GasFraction");
				_Rho.writeMED(resultFile+"MixtureDensity");
				break;
			case CSV :
				_VV.writeCSV(resultFile+"Enthalpy");
				_TT.writeCSV(resultFile+"Temperature");
				_Alpha.writeCSV(resultFile+"GasFraction");
				_Rho.writeCSV(resultFile+"MixtureDensity");
				break;
			}
		}
		// do not create mesh
		else{
			switch(_saveFormat)
			{
			case VTK :
				_VV.writeVTK(resultFile+"Enthalpy",false);
				_TT.writeVTK(resultFile+"Temperature",false);
				_Alpha.writeVTK(resultFile+"GasFraction",false);
				_Rho.writeVTK(resultFile+"MixtureDensity",false);
				break;
			case MED :
				_VV.writeMED(resultFile+"Enthalpy",false);
				_TT.writeMED(resultFile+"Temperature",false);
				_Alpha.writeMED(resultFile+"GasFraction",false);
				_Rho.writeMED(resultFile+"MixtureDensity",false);
				break;
			case CSV :
				_VV.writeCSV(resultFile+"Enthalpy");
				_TT.writeCSV(resultFile+"Temperature");
				_Alpha.writeCSV(resultFile+"GasFraction");
				_Rho.writeCSV(resultFile+"MixtureDensity");
				break;
			}
		}
	}
}

vector<string> TransportEquation::getInputFieldsNames()
{
	vector<string> result(2);
	
	result[0]="HeatPower";
	result[1]="RodTemperature";
	
	return result;
}

vector<string> TransportEquation::getOutputFieldsNames()
{
	vector<string> result(4);
	
	result[0]="Enthalpy";
	result[1]="FluidTemperature";
	result[2]="VoidFraction";
	result[3]="Density";
	
	return result;
}

Field& TransportEquation::getOutputField(const string& nameField )
{
	if(nameField=="FluidTemperature" || nameField=="FLUIDTEMPERATURE" || nameField=="TemperatureFluide" || nameField=="TEMPERATUREFLUIDE")
		return getFluidTemperatureField();
	else if(nameField=="Enthalpy" || nameField=="ENTHALPY" || nameField=="Enthalpie" || nameField=="ENTHALPIE" )
		return getEnthalpyField();
    else 	if(nameField=="VoidFraction" || nameField=="VOIDFRACTION" || nameField=="TauxDeVide" || nameField=="TAUXDEVIDE")
		return getVoidFractionField();
	else if(nameField=="Density" || nameField=="DENSITY" || nameField=="Densité" || nameField=="DENSITE" )
		return getDensityField();
	else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first" << endl;
        throw CdmathException("TransportEquation::getOutputField error : Unknown Field name");
    }
}

void
TransportEquation::setInputField(const string& nameField, Field& inputField )
{
	if(nameField=="RodTemperature" || nameField=="RODTEMPERATURE" || nameField=="TemperatureCombustible" || nameField=="TEMPERATURECOMBUSTIBLE")
		return setRodTemperatureField( inputField) ;
	else if(nameField=="HeatPower" || nameField=="HEATPOWER" || nameField=="PuissanceThermique" || nameField=="PUISSANCETHERMIQUE" )
		return setHeatPowerField( inputField );
	else
    {
        cout<<"Error : Field name "<< nameField << " is not an input field name, call getInputFieldsNames first" << endl;
        throw CdmathException("TransportEquation::setInputField error : Unknown Field name");
    }
}
