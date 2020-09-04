#include "TransportEquation.hxx"
#include "math.h"
#include <fstream>
#include <sstream>

using namespace std;

TransportEquation::TransportEquation(phase fluid, pressureMagnitude pEstimate,vector<double> vitesseTransport){
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
}

void TransportEquation::initialize()
{
	if(!_initialDataSet)
		throw CdmathException("TransportEquation::initialize() set initial data first");
	else
		cout<<"Initialising the transport of a fluid enthalpy"<<endl;
	/**************** Field creation *********************/

	//post processing fields used only for saving results
	_TT=Field ("Temperature", CELLS, _mesh, 1);
	_Alpha=Field ("Void fraction", CELLS, _mesh, 1);
	_Rho=Field ("Mixture density", CELLS, _mesh, 1);
	//Construction des champs de post-traitement
	VecCreate(PETSC_COMM_SELF, &_Hn);
	VecSetSizes(_Hn,PETSC_DECIDE,_Nmailles);
	VecSetFromOptions(_Hn);
	for(int i =0; i<_Nmailles;i++){
		_TT(i)=temperature(_VV(i));
		_Alpha(i)=voidFraction(_VV(i));
		_Rho(i)=density(_Alpha(i));
		VecSetValue(_Hn,i,_VV(i), INSERT_VALUES);
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

	//creation de la matrice
	MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles, _Nmailles, (1+_neibMaxNb), PETSC_NULL, &_A);
	VecDuplicate(_Hn, &_Hk);
	VecDuplicate(_Hn, &_Hkm1);
	VecDuplicate(_Hn, &_deltaH);
	VecDuplicate(_Hn, &_b);//RHS of the linear system: _b=Hn/dt + _b0 + puisance
	VecDuplicate(_Hn, &_b0);//part of the RHS that comes from the boundary conditions

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

	stop=false;
	return min(_dt_transport,_dt_src);
}
double TransportEquation::computeTransportMatrix(){
	long nbFaces = _mesh.getNumberOfFaces();
	Face Fj;
	Cell Cell1,Cell2;
	string nameOfGroup;
	double inv_dxi, inv_dxj;
	Vector normale(_Ndim);
	double un, hk;
	PetscInt idm, idn;
	std::vector< int > idCells;
	MatZeroEntries(_A);
	VecZeroEntries(_b0);
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
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "face numero " << j << " cellule frontiere " << idCells[0] << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << normale[p] << ",";
				cout << ") "<<endl;
			}
			nameOfGroup = Fj.getGroupName();

			if (_limitField[nameOfGroup].bcType==NeumannTransport){
				MatSetValue(_A,idm,idm,inv_dxi*un, ADD_VALUES);
			}
			else if(_limitField[nameOfGroup].bcType==InletTransport){
				if(un>0){
					MatSetValue(_A,idm,idm,inv_dxi*un, ADD_VALUES);
				}
				else{
					hk=_limitField[nameOfGroup].h;
					VecSetValue(_b0,idm,-inv_dxi*un*hk, ADD_VALUES);
				}
			}
			else {
				cout<<"!!!!!!!!!!!!!!! Error TransportEquation::computeTransportMatrix() !!!!!!!!!!"<<endl;
				cout<<"!!!!!!!!! Boundary condition not treated for boundary named "<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<" !!!!!!!!!!!!!! "<<endl;
				cout<<"Accepted boundary conditions are NeumannTransport "<<NeumannTransport<< " and InletTransport "<< InletTransport <<endl;
				throw CdmathException("Boundary condition not accepted");
			}
			// if Fj is inside the domain
		} else 	if (Fj.getNumberOfCells()==2 ){
			if(_verbose && _nbTimeStep%_freqSave ==0)
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
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);
	VecAssemblyEnd(_b0);
	_transportMatrixSet=true;
	if(abs(_maxvp)<_precision)
		throw CdmathException("TransportEquation::computeTransportMatrix(): maximum eigenvalue for time step is zero");
	else
		return _cfl*_minl/_maxvp;
}
double TransportEquation::computeRHS(){
	double rhomin=INFINITY;
	VecCopy(_b0,_b);
	VecAssemblyBegin(_b);
	if(_system)
		cout<<"second membre of transport problem"<<endl;
	for (int i=0; i<_Nmailles;i++){
		VecSetValue(_b,i,_heatTransfertCoeff*(_rodTemperatureField(i)-_TT(i))/_Rho(i),ADD_VALUES);
		VecSetValue(_b,i,_heatPowerField(i)/_Rho(i),ADD_VALUES);
		if(_system)
			cout<<_heatPowerField(i)/_Rho(i)<<endl;
		if(_Rho(i)<rhomin)
			rhomin=_Rho(i);
	}
	VecAssemblyEnd(_b);
	if(_system)
		VecView(_b,  PETSC_VIEWER_STDOUT_WORLD);

	return rhomin*_cpref/_heatTransfertCoeff;
}
void TransportEquation::updatePrimitives()
{
	double hi;
	for(int i=0; i<_Nmailles; i++)
	{
		VecGetValues(_Hk, 1, &i, &hi);
		_VV(i)=hi;
		_TT(i)=temperature(hi);
		_Alpha(i)=voidFraction(hi);
		_Rho(i)=density(_Alpha(i));
	}
}

bool TransportEquation::initTimeStep(double dt){
	_dt = dt;
	if(_verbose && _nbTimeStep%_freqSave ==0)
		MatView(_A,PETSC_VIEWER_STDOUT_SELF);
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

	if(_timeScheme == Implicit)
		MatShift(_A,1/_dt);

	return _dt>0;
}

void TransportEquation::abortTimeStep(){
	VecAXPY(_b,  -1/_dt, _Hn);
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

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
		cout << "Vecteur Hn: " << endl;
		VecView(_Hn,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
		cout<<"Vecteur _b "<<endl;
		VecView(_b,  PETSC_VIEWER_STDOUT_SELF);
		cout << "Matrice A "<<endl;
		MatView(_A,PETSC_VIEWER_STDOUT_SELF);
	}

	if(_timeScheme == Explicit)
	{
		MatMult(_A, _Hn, _Hk);
		if(_system)
		{
			cout << "Nouveau vecteur Hk: " << endl;
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		VecAXPY(_Hk, -1, _b);
		if(_system)
		{
			cout << "Nouveau vecteur Hk-b: " << endl;
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		VecScale(_Hk, -_dt);
		if(_system)
		{
			cout << "Nouveau vecteur dt*(Hk-b): " << endl;
			VecView(_Hk,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
		converged = true;
	}
	else
	{
		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

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
			cout<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
			converged=false;
			return false;
		}
		else
		{
			VecCopy(_Hk, _deltaH);//ici on a deltaH=Hk
			VecAXPY(_deltaH,  -1, _Hkm1);//On obtient deltaH=Hk-Hkm1
			_erreur_rel= 0;
			double hi, dhi;

			for(int i=0; i<_Nmailles; i++)
			{
				VecGetValues(_deltaH, 1, &i, &dhi);
				VecGetValues(_Hk, 1, &i, &hi);
				if(_erreur_rel < fabs(dhi/hi))
					_erreur_rel = fabs(dhi/hi);
			}
		}

		converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
	}

	updatePrimitives();

	VecCopy(_Hk, _Hkm1);


	return true;
}
void TransportEquation::validateTimeStep()
{
	VecCopy(_Hk, _deltaH);//ici Hk=Hnp1 donc on a deltaH=Hnp1
	VecAXPY(_deltaH,  -1, _Hn);//On obtient deltaH=Hnp1-Hn

	_erreur_rel= 0;
	double hi, dhi;

	for(int i=0; i<_Nmailles; i++)
	{
		VecGetValues(_deltaH, 1, &i, &dhi);
		VecGetValues(_Hk, 1, &i, &hi);
		if(_erreur_rel < fabs(dhi/hi))
			_erreur_rel = fabs(dhi/hi);
	}
	_isStationary =(_erreur_rel <_precision);

	VecCopy(_Hk, _Hn);
	VecCopy(_Hk, _Hkm1);

	if(_nbTimeStep%_freqSave ==0)
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
}

void TransportEquation::save(){
	string resultFile(_path+"/TransportEquation_");///Results
	resultFile+=_fileName;

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

vector<string> TransportEquation::getOutputFieldsNames()
{
	vector<string> result(2);
	
	result[0]="Enthalpy";
	result[1]="FluidTemperature";
	
	return result;
}

Field& TransportEquation::getOutputField(const string& nameField )
{
	if(nameField=="FluidTemperature" || nameField=="FLUIDTEMPERATURE" )
		return getFluidTemperatureField();
	else if(nameField=="Enthalpy" || nameField=="ENTHALPY" || nameField=="Enthalpie" || nameField=="ENTHALPY" )
		return getEnthalpyField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first" << endl;
        throw CdmathException("TransportEquation::getOutputField error : Unknown Field name");
    }
}

