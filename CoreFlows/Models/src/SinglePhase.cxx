/*
 * SinglePhase.cxx
 *
 *  Created on: Sep 15, 2014
 *      Author: tn236279
 */

#include "SinglePhase.hxx"

using namespace std;

SinglePhase::SinglePhase(phaseType fluid, pressureEstimate pEstimate, int dim, bool useDellacherieEOS){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1;
	_dragCoeffs=vector<double>(1,0);
	_fluides.resize(1);
	_useDellacherieEOS=useDellacherieEOS;
	_saveAllFields=false;

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		_Tref=300;
		_Pref=1e5;
		if(fluid==Gas){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(1.4,743,_Tref,2.22e5);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, e=2.22e5, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(996,_Pref,_Tref,1.12e5,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C, e=1.12e5, c_v=4130
		}
	}
	else{//EOS at 155 bars and 618K 
		_Tref=618;
		_Pref=155e5;
		if(fluid==Gas){
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_fluides[0] = new StiffenedGas(102,_Pref,_Tref,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
		}
		else{//To do : change to normal regime: 155 bars and 573K
			cout<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			if(_useDellacherieEOS)
				_fluides[0]= new StiffenedGasDellacherie(2.35,1e9,-1.167e6,1816); //stiffened gas law for water from S. Dellacherie
			else
				_fluides[0]= new StiffenedGas(594.,_Pref,_Tref,1.6e6, 621.,3100.);  //stiffened gas law for water at pressure 155 bar, and temperature 345°C
		}
	}
}
void SinglePhase::initialize(){
	cout<<"Initialising the Navier-Stokes model"<<endl;
	*_runLogFile<<"Initialising the Navier-Stokes model"<<endl;

	_Uroe = new double[_nVar];
	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+1]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];
	if(_saveVelocity || _saveAllFields)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant

	if(_saveAllFields)
	{
		_Enthalpy=Field("Enthalpy",CELLS,_mesh,1);
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		_Density=Field("Density",CELLS,_mesh,1);
		_Temperature=Field("Temperature",CELLS,_mesh,1);
		_VitesseX=Field("Velocity x",CELLS,_mesh,1);
		if(_Ndim>1)
		{
			_VitesseY=Field("Velocity y",CELLS,_mesh,1);
			if(_Ndim>2)
				_VitesseZ=Field("Velocity z",CELLS,_mesh,1);
		}
	}

	if(_entropicCorrection)
		_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

	ProblemFluid::initialize();
}

bool SinglePhase::iterateTimeStep(bool &converged)
{
	if(_timeScheme == Explicit || !_usePrimitiveVarsInNewton)
		return ProblemFluid::iterateTimeStep(converged);
	else
	{
		bool stop=false;

		if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
			_maxvp=0;
			computeTimeStep(stop);
		}
		if(stop){//Le compute time step ne s'est pas bien passé
			cout<<"ComputeTimeStep failed"<<endl;
			converged=false;
			return false;
		}

		computeNewtonVariation();

		//converged=convergence des iterations
		if(_timeScheme == Explicit)
			converged=true;
		else{//Implicit scheme

			KSPGetIterationNumber(_ksp, &_PetscIts);
			if( _MaxIterLinearSolver < _PetscIts)//save the maximum number of iterations needed during the newton scheme
				_MaxIterLinearSolver = _PetscIts;
			if(_PetscIts>=_maxPetscIts)//solving the linear system failed
			{
				cout<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
				*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
				converged=false;
				return false;
			}
			else{//solving the linear system succeeded
				//Calcul de la variation relative Uk+1-Uk
				_erreur_rel = 0.;
				double x, dx;
				int I;
				for(int j=1; j<=_Nmailles; j++)
				{
					for(int k=0; k<_nVar; k++)
					{
						I = (j-1)*_nVar + k;
						VecGetValues(_newtonVariation, 1, &I, &dx);
						VecGetValues(_primitiveVars, 1, &I, &x);
						if (fabs(x)*fabs(x)< _precision)
						{
							if(_erreur_rel < fabs(dx))
								_erreur_rel = fabs(dx);
						}
						else if(_erreur_rel < fabs(dx/x))
							_erreur_rel = fabs(dx/x);
					}
				}
			}
			converged = _erreur_rel <= _precision_Newton;
		}

		double relaxation=1;//Vk+1=Vk+relaxation*deltaV

		VecAXPY(_primitiveVars,  relaxation, _newtonVariation);

		//mise a jour du champ primitif
		updateConservatives();

		if(_nbPhases==2 && fabs(_err_press_max) > _precision)//la pression n'a pu être calculée en diphasique à partir des variables conservatives
		{
			cout<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
			*_runLogFile<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
			converged=false;
			return false;
		}
		if(_system)
		{
			cout<<"Vecteur Vkp1-Vk "<<endl;
			VecView(_newtonVariation,  PETSC_VIEWER_STDOUT_SELF);
			cout << "Nouvel etat courant Vk de l'iteration Newton: " << endl;
			VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_SELF);
		}

		if(_nbPhases==2 && _nbTimeStep%_freqSave ==0){
			if(_minm1<-_precision || _minm2<-_precision)
			{
				cout<<"!!!!!!!!! WARNING masse partielle negative sur " << _nbMaillesNeg << " faces, min m1= "<< _minm1 << " , minm2= "<< _minm2<< " precision "<<_precision<<endl;
				*_runLogFile<<"!!!!!!!!! WARNING masse partielle negative sur " << _nbMaillesNeg << " faces, min m1= "<< _minm1 << " , minm2= "<< _minm2<< " precision "<<_precision<<endl;
			}

			if (_nbVpCplx>0){
				cout << "!!!!!!!!!!!!!!!!!!!!!!!! Complex eigenvalues on " << _nbVpCplx << " cells, max imag= " << _part_imag_max << endl;
				*_runLogFile << "!!!!!!!!!!!!!!!!!!!!!!!! Complex eigenvalues on " << _nbVpCplx << " cells, max imag= " << _part_imag_max << endl;
			}
		}
		_minm1=1e30;
		_minm2=1e30;
		_nbMaillesNeg=0;
		_nbVpCplx =0;
		_part_imag_max=0;

		return true;
	}
}
void SinglePhase::computeNewtonVariation()
{
	if(!_usePrimitiveVarsInNewton)
		ProblemFluid::computeNewtonVariation();
	else
	{
		if(_verbose)
		{
			cout<<"Vecteur courant Vk "<<endl;
			VecView(_primitiveVars,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Matrice du système linéaire avant contribution delta t" << endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Second membre du système linéaire avant contribution delta t" << endl;
			VecView(_b, PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
		if(_timeScheme == Explicit)
		{
			VecCopy(_b,_newtonVariation);
			VecScale(_newtonVariation, _dt);
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout<<"Vecteur _newtonVariation =_b*dt"<<endl;
				VecView(_newtonVariation,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
		}
		else
		{
			MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);

			VecAXPY(_b, 1/_dt, _old);
			VecAXPY(_b, -1/_dt, _conservativeVars);

			for(int imaille = 0; imaille<_Nmailles; imaille++){
				_idm[0] = _nVar*imaille;
				for(int k=1; k<_nVar; k++)
					_idm[k] = _idm[k-1] + 1;
				VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
				primToConsJacobianMatrix(_Vi);
				for(int k=0; k<_nVar*_nVar; k++)
					_primToConsJacoMat[k]*=1/_dt;
				MatSetValuesBlocked(_A, 1, &imaille, 1, &imaille, _primToConsJacoMat, ADD_VALUES);
			}
			MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

#if PETSC_VERSION_GREATER_3_5
			KSPSetOperators(_ksp, _A, _A);
#else
			KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

			if(_verbose)
			{
				cout << "Matrice du système linéaire" << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
				cout << "Second membre du système linéaire" << endl;
				VecView(_b, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}

			if(_conditionNumber)
				KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
			if(!_isScaling)
			{
				KSPSolve(_ksp, _b, _newtonVariation);
			}
			else
			{
				computeScaling(_maxvp);
				int indice;
				VecAssemblyBegin(_vecScaling);
				VecAssemblyBegin(_invVecScaling);
				for(int imaille = 0; imaille<_Nmailles; imaille++)
				{
					indice = imaille;
					VecSetValuesBlocked(_vecScaling,1 , &indice, _blockDiag, INSERT_VALUES);
					VecSetValuesBlocked(_invVecScaling,1,&indice,_invBlockDiag, INSERT_VALUES);
				}
				VecAssemblyEnd(_vecScaling);
				VecAssemblyEnd(_invVecScaling);

				if(_system)
				{
					cout << "Matrice avant le preconditionneur des vecteurs propres " << endl;
					MatView(_A,PETSC_VIEWER_STDOUT_SELF);
					cout << endl;
					cout << "Second membre avant le preconditionneur des vecteurs propres " << endl;
					VecView(_b, PETSC_VIEWER_STDOUT_SELF);
					cout << endl;
				}
				MatDiagonalScale(_A,_vecScaling,_invVecScaling);
				if(_system)
				{
					cout << "Matrice apres le preconditionneur des vecteurs propres " << endl;
					MatView(_A,PETSC_VIEWER_STDOUT_SELF);
					cout << endl;
				}
				VecCopy(_b,_bScaling);
				VecPointwiseMult(_b,_vecScaling,_bScaling);
				if(_system)
				{
					cout << "Produit du second membre par le preconditionneur bloc diagonal  a gauche" << endl;
					VecView(_b, PETSC_VIEWER_STDOUT_SELF);
					cout << endl;
				}

				KSPSolve(_ksp,_b, _bScaling);
				VecPointwiseMult(_newtonVariation,_invVecScaling,_bScaling);
			}
			if(_system)
			{
				cout << "solution du systeme lineaire local:" << endl;
				VecView(_newtonVariation, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
		}
	}
}
void SinglePhase::convectionState( const long &i, const long &j, const bool &IsBord){

	_idm[0] = _nVar*i; 
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);

	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	if(IsBord)
		VecGetValues(_Uext, _nVar, _idm, _Uj);
	else
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Ui[k]<<endl;
		cout<<"Convection Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Uj[k]<<endl;
	}
	if(_Ui[0]<0||_Uj[0]<0)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!densite negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!densite negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("densite negative, arret de calcul");
	}
	PetscScalar ri, rj, xi, xj, pi, pj;
	PetscInt Ii;
	ri = sqrt(_Ui[0]);//racine carre de phi_i rho_i
	rj = sqrt(_Uj[0]);
	_Uroe[0] = ri*rj;	//moyenne geometrique des densites
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Densité moyenne Roe  gauche " << i << ": " << ri*ri << ", droite " << j << ": " << rj*rj << "->" << _Uroe[0] << endl;
	for(int k=0;k<_Ndim;k++){
		xi = _Ui[k+1];
		xj = _Uj[k+1];
		_Uroe[1+k] = (xi/ri + xj/rj)/(ri + rj);
		//"moyenne" des vitesses
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout << "Vitesse de Roe composante "<< k<<"  gauche " << i << ": " << xi/(ri*ri) << ", droite " << j << ": " << xj/(rj*rj) << "->" << _Uroe[k+1] << endl;
	}
	// H = (rho E + p)/rho
	xi = _Ui[_nVar-1];//phi rho E
	xj = _Uj[_nVar-1];
	Ii = i*_nVar; // correct Kieu
	VecGetValues(_primitiveVars, 1, &Ii, &pi);// _primitiveVars pour p
	if(IsBord)
	{
		double q_2 = 0;
		for(int k=1;k<=_Ndim;k++)
			q_2 += _Uj[k]*_Uj[k];
		q_2 /= _Uj[0];	//phi rho u²
		pj =  _fluides[0]->getPressure((_Uj[(_Ndim+2)-1]-q_2/2)/_porosityj,_Uj[0]/_porosityj);
	}
	else
	{
		Ii = j*_nVar; // correct Kieu
		VecGetValues(_primitiveVars, 1, &Ii, &pj);
	}
	xi = (xi + pi)/(ri*ri);
	xj = (xj + pj)/(rj*rj);
	_Uroe[_nVar-1] = (ri*xi + rj*xj)/(ri + rj);
	//on se donne l enthalpie ici
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Enthalpie totale de Roe H  gauche " << i << ": " << xi << ", droite " << j << ": " << xj << "->" << _Uroe[_nVar-1] << endl;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
}

void SinglePhase::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat de diffusion (rho, q, T)
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	if(IsBord)
		VecGetValues(_Uextdiff, _nVar, _idm, _Uj);
	else
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "SinglePhase::diffusionStateAndMatrices cellule gauche" << i << endl;
		cout << "Ui = ";
		for(int q=0; q<_nVar; q++)
			cout << _Ui[q]  << "\t";
		cout << endl;
		cout << "SinglePhase::diffusionStateAndMatrices cellule droite" << j << endl;
		cout << "Uj = ";
		for(int q=0; q<_nVar; q++)
			cout << _Uj[q]  << "\t";
		cout << endl;
	}

	for(int k=0; k<_nVar; k++)
		_Udiff[k] = (_Ui[k]/_porosityi+_Uj[k]/_porosityj)/2;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "SinglePhase::diffusionStateAndMatrices conservative diffusion state" << endl;
		cout << "_Udiff = ";
		for(int q=0; q<_nVar; q++)
			cout << _Udiff[q]  << "\t";
		cout << endl;
		cout << "porosite gauche= "<<_porosityi<< ", porosite droite= "<<_porosityj<<endl;
	}
	consToPrim(_Udiff,_phi,1);
	_Udiff[_nVar-1]=_phi[_nVar-1];

	if(_timeScheme==Implicit)
	{
		double q_2=0;
		for (int i = 0; i<_Ndim;i++)
			q_2+=_Udiff[i+1]*_Udiff[i+1];
		double mu = _fluides[0]->getViscosity(_Udiff[_nVar-1]);
		double lambda = _fluides[0]->getConductivity(_Udiff[_nVar-1]);
		double Cv= _fluides[0]->constante("Cv");
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		for(int i=1;i<(_nVar-1);i++)
		{
			_Diffusion[i*_nVar] =  mu*_Udiff[i]/(_Udiff[0]*_Udiff[0]);
			_Diffusion[i*_nVar+i] = -mu/_Udiff[0];
		}
		int i = (_nVar-1)*_nVar;
		_Diffusion[i]=lambda*(_Udiff[_nVar-1]/_Udiff[0]-q_2/(2*Cv*_Udiff[0]*_Udiff[0]*_Udiff[0]));
		for(int k=1;k<(_nVar-1);k++)
		{
			_Diffusion[i+k]= lambda*_Udiff[k]/(_Udiff[0]*_Udiff[0]*Cv);
		}
		_Diffusion[_nVar*_nVar-1]=-lambda/(_Udiff[0]*Cv);
	}

}
void SinglePhase::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	double porosityj=_porosityField(j);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< " face unit normal vector "<<endl;
		for(int k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne conservatif
		double q_n=0;//q_n=quantité de mouvement normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_n+=_externalStates[(k+1)]*normale[k];
			
		//Pour la convection, inversion du sens de la vitesse normale
		for(int k=0; k<_Ndim; k++)
			_externalStates[(k+1)]-= 2*q_n*normale[k];

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		//Pour la diffusion, paroi à vitesse et temperature imposees
		_idm[0] = _nVar*j;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//L'état fantome contient à présent les variables primitives internes
		double pression=_externalStates[0];
		double T=_limitField[nameOfGroup].T;
		double rho=_fluides[0]->getDensity(pression,T);

		_externalStates[0]=porosityj*rho;
		_externalStates[1]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
		double v2=0;
		v2 +=_limitField[nameOfGroup].v_x[0]*_limitField[nameOfGroup].v_x[0];
		if(_Ndim>1)
		{
			v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
			if(_Ndim==3)
			{
				_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
				v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
			}
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On prend l'état fantôme égal à l'état interne (conditions limites de Neumann)

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne (conditions limites de Neumann)
		double q_int_n=0;//q_int_n=composante normale de la quantité de mouvement à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_int_n+=_externalStates[(k+1)]*normale[k];//On calcule la vitesse normale sortante

		double q_ext_n=_limitField[nameOfGroup].v_x[0]*normale[0];
		if(_Ndim>1)
			{
				q_ext_n+=_limitField[nameOfGroup].v_y[0]*normale[1];
				if(_Ndim>2)
						q_ext_n+=_limitField[nameOfGroup].v_z[0]*normale[2];
			}

		if(q_int_n+q_ext_n<=0){//Interfacial velocity goes inward
			VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On met à jour l'état fantome avec les variables primitives internes
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_fluides[0]->getDensity(pression,T);

			_externalStates[0]=porosityj*rho;//Composante fantome de masse
			_externalStates[1]=_externalStates[0]*(_limitField[nameOfGroup].v_x[0]);//Composante fantome de qdm x
			double v2=0;
			v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
			if(_Ndim>1)
			{
				v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
				_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];//Composante fantome de qdm y
				if(_Ndim==3)
				{
					_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];//Composante fantome de qdm z
					v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);//Composante fantome de l'nrj
		}
		else if(_nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid possibly going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletRotationVelocity){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double u_int_n=0;//u_int_n=composante normale de la vitesse intérieure à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_int_n+=_externalStates[(k+1)]*normale[k];

		double u_ext_n=0;
        //ghost velocity
        if(_Ndim>1)
        {
            Vector omega(3);
            omega[0]=_limitField[nameOfGroup].v_x[0];
            omega[1]=_limitField[nameOfGroup].v_y[0];
            
            Vector Normale(3);
            Normale[0]=normale[0];
            Normale[1]=normale[1];

            if(_Ndim==3)
            {
                omega[2]=_limitField[nameOfGroup].v_z[0];
                Normale[2]=normale[2];
            }
            
            Vector tangent_vel=omega%Normale;
			u_ext_n=-0.01*tangent_vel.norm();
			//Changing external state velocity
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]=u_ext_n*normale[k] + tangent_vel[k];
        }

		if(u_ext_n + u_int_n <= 0){
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_fluides[0]->getDensity(pression,T);

			double v2=0;
			v2 +=_externalStates[1]*_externalStates[1];//v_x*v_x
			_externalStates[0]=porosityj*rho;
			_externalStates[1]*=_externalStates[0];
			if(_Ndim>1)
			{
				v2 +=_externalStates[2]*_externalStates[2];//+v_y*v_y
				_externalStates[2]*=_externalStates[0];
				if(_Ndim==3)
				{
					v2 +=_externalStates[3]*_externalStates[3];//+v_z*v_z
					_externalStates[3]*=_externalStates[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		}
		else if(_nbTimeStep%_freqSave ==0)
		{
			/*
			 * cout<< "Warning : fluid going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;
			 */ 
			VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On définit l'état fantôme avec l'état interne
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Wall boundary condition."<<endl;
			
			//Changing external state momentum
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]-=2*_externalStates[0]*u_int_n*normale[k];
		}

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		//Building the primitive external state
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double u_n=0;//u_n=vitesse normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_n+=_externalStates[(k+1)]*normale[k];
        
		if(u_n<=0)
		{
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress,_limitField[nameOfGroup].T);
			
	        //Contribution from the tangential velocity
	        if(_Ndim>1)
	        {
	            Vector omega(3);
	            omega[0]=_limitField[nameOfGroup].v_x[0];
	            omega[1]=_limitField[nameOfGroup].v_y[0];
	            
	            Vector Normale(3);
	            Normale[0]=normale[0];
	            Normale[1]=normale[1];
	
	            if(_Ndim==3)
	            {
	                omega[2]=_limitField[nameOfGroup].v_z[0];
	                Normale[2]=normale[2];
	            }
	            
	            Vector tangent_vel=omega%Normale;
	
				//Changing external state velocity
	            for(int k=0; k<_Ndim; k++)
	                _externalStates[(k+1)]=u_n*normale[k] + abs(u_n)*tangent_vel[k];
	        }
		}
		else{
			/*
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature (only pressure value is imposed as in outlet BC)."<<endl;
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
			*/
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Wall boundary condition."<<endl;
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_externalStates[0]+hydroPress, _externalStates[_nVar-1]);
			//Changing external state velocity
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]-=2*u_n*normale[k];
		}

		double v2=0;
		for(int k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;//qdm component
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);//nrj component


		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne conservatif
		double q_n=0;//q_n=quantité de mouvement normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_n+=_externalStates[(k+1)]*normale[k];

		if(q_n < -_precision &&  _nbTimeStep%_freqSave ==0)
        {
			cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<" with flow rate "<< q_n<<endl;
            cout<< "Applying Neumann boundary condition for velocity and temperature"<<endl;
        }
		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"Cond lim outlet densite= "<<_externalStates[0]<<" gravite= "<<_GravityField3d[0]<<" Cj.x()= "<<Cj.x()<<endl;
			cout<<"Cond lim outlet pression ref= "<<_limitField[nameOfGroup].p<<" pression hydro= "<<hydroPress<<" total= "<<_limitField[nameOfGroup].p+hydroPress<<endl;
		}
		//Building the external state
		_idm[0] = _nVar*j;// Kieu
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		double v2=0;
		for(int k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);
		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<< " _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}
}

void SinglePhase::convectionMatrices()
{
	//entree: URoe = rho, u, H
	//sortie: matrices Roe+  et Roe-

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhase::convectionMatrices()"<<endl;

	double u_n=0, u_2=0;//vitesse normale et carré du module

	for(int i=0;i<_Ndim;i++)
	{
		u_2 += _Uroe[1+i]*_Uroe[1+i];
		u_n += _Uroe[1+i]*_vec_normal[i];
	}

	vector<complex<double> > vp_dist(3,0);

	if(_spaceScheme==staggered && _nonLinearFormulation==VFFC)//special case
	{
		staggeredVFFCMatricesConservativeVariables(u_n);//Computation of classical upwinding matrices
		if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)//For use in implicit matrix
			staggeredVFFCMatricesPrimitiveVariables(u_n);
	}
	else
	{
		Vector vitesse(_Ndim);
		for(int idim=0;idim<_Ndim;idim++)
			vitesse[idim]=_Uroe[1+idim];

		double  c, H, K, k;
		/***********Calcul des valeurs propres ********/
		H = _Uroe[_nVar-1];
		c = _fluides[0]->vitesseSonEnthalpie(H-u_2/2);//vitesse du son a l'interface
		k = _fluides[0]->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
		K = u_2*k/2; //g-1/2 *|u|²

		vp_dist[0]=u_n-c;vp_dist[1]=u_n;vp_dist[2]=u_n+c;

		_maxvploc=fabs(u_n)+c;
		if(_maxvploc>_maxvp)
			_maxvp=_maxvploc;

		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<"SinglePhase::convectionMatrices Eigenvalues "<<u_n-c<<" , "<<u_n<<" , "<<u_n+c<<endl;

		RoeMatrixConservativeVariables( u_n, H,vitesse,k,K);

		/******** Construction des matrices de decentrement ********/
		if( _spaceScheme ==centered){
			if(_entropicCorrection)
			{
				*_runLogFile<<"SinglePhase::convectionMatrices: entropy scheme not available for centered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhase::convectionMatrices: entropy scheme not available for centered scheme");
			}

			for(int i=0; i<_nVar*_nVar;i++)
				_absAroe[i] = 0;
		}
		else if(_spaceScheme == upwind || _spaceScheme ==pressureCorrection || _spaceScheme ==lowMach){
			if(_entropicCorrection)
				entropicShift(_vec_normal);
			else
				_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

			vector< complex< double > > y (3,0);
			Polynoms Poly;
			for( int i=0 ; i<3 ; i++)
				y[i] = Poly.abs_generalise(vp_dist[i])+_entropicShift[i];
			Poly.abs_par_interp_directe(3,vp_dist, _Aroe, _nVar,_precision, _absAroe,y);

			if( _spaceScheme ==pressureCorrection)
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(1+i)*_nVar+1+j]-=(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			else if( _spaceScheme ==lowMach){
				double M=sqrt(u_2)/c;
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(1+i)*_nVar+1+j]-=(1-M)*(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			}
		}
		else if( _spaceScheme ==staggered ){
			if(_entropicCorrection)//To do: study entropic correction for staggered
			{
				*_runLogFile<<"SinglePhase::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhase::convectionMatrices: entropy scheme not available for staggered scheme");
			}

			staggeredRoeUpwindingMatrixConservativeVariables( u_n, H, vitesse, k, K);
		}
		else
		{
			*_runLogFile<<"SinglePhase::convectionMatrices: scheme not treated"<<endl;
			_runLogFile->close();
			throw CdmathException("SinglePhase::convectionMatrices: scheme not treated");
		}

		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinus[i] = (_Aroe[i]-_absAroe[i])/2;
			_AroePlus[i]  = (_Aroe[i]+_absAroe[i])/2;
		}
		if(_timeScheme==Implicit)
		{
			if(_usePrimitiveVarsInNewton)//Implicitation using primitive variables
			{
				_Vij[0]=_fluides[0]->getPressureFromEnthalpy(_Uroe[_nVar-1]-u_2/2, _Uroe[0]);//pressure
				_Vij[_nVar-1]=_fluides[0]->getTemperatureFromPressure( _Vij[0], _Uroe[0]);//Temperature
				for(int idim=0;idim<_Ndim; idim++)
					_Vij[1+idim]=_Uroe[1+idim];
				primToConsJacobianMatrix(_Vij);
				Polynoms Poly;
				Poly.matrixProduct(_AroeMinus, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroeMinusImplicit);
				Poly.matrixProduct(_AroePlus,  _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroePlusImplicit);
			}
			else
				for(int i=0; i<_nVar*_nVar;i++)
				{
					_AroeMinusImplicit[i] = _AroeMinus[i];
					_AroePlusImplicit[i]  = _AroePlus[i];
				}
		}
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			displayMatrix(_Aroe, _nVar,"Matrice de Roe");
			displayMatrix(_absAroe, _nVar,"Valeur absolue matrice de Roe");
			displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus");
			displayMatrix(_AroePlus, _nVar,"Matrice _AroePlus");
		}
	}

	if(_verbose && _nbTimeStep%_freqSave ==0 && _timeScheme==Implicit)
	{
		displayMatrix(_AroeMinusImplicit, _nVar,"Matrice _AroeMinusImplicit");
		displayMatrix(_AroePlusImplicit,  _nVar,"Matrice _AroePlusImplicit");
	}

	/*********Calcul de la matrice signe pour VFFC, VFRoe et décentrement des termes source*****/
	if(_entropicCorrection)
	{
		InvMatriceRoe( vp_dist);
		Polynoms Poly;
		Poly.matrixProduct(_absAroe, _nVar, _nVar, _invAroe, _nVar, _nVar, _signAroe);
	}
	else if (_spaceScheme==upwind || (_spaceScheme ==pressureCorrection ) || (_spaceScheme ==lowMach ))//upwind sans entropic
		SigneMatriceRoe( vp_dist);
	else if (_spaceScheme==centered)//centre  sans entropic
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
	else if( _spaceScheme ==staggered )//à tester
	{
		double signu;
		if(u_n>0)
			signu=1;
		else if (u_n<0)
			signu=-1;
		else
			signu=0;
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
		_signAroe[0] = signu;
		for(int i=1; i<_nVar-1;i++)
			_signAroe[i*_nVar+i] = -signu;
		_signAroe[_nVar*(_nVar-1)+_nVar-1] = signu;
	}
	else
	{
		*_runLogFile<<"SinglePhase::convectionMatrices: well balanced option not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::convectionMatrices: well balanced option not treated");
	}
}
void SinglePhase::computeScaling(double maxvp)
{
	_blockDiag[0]=1;
	_invBlockDiag[0]=1;
	for(int q=1; q<_nVar-1; q++)
	{
		_blockDiag[q]=1./maxvp;//
		_invBlockDiag[q]= maxvp;//1.;//
	}
	_blockDiag[_nVar - 1]=(_fluides[0]->constante("gamma")-1)/(maxvp*maxvp);//1
	_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
}

void SinglePhase::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)

{
	double lambda=_fluides[0]->getConductivity(_Udiff[_nVar-1]);
	double mu = _fluides[0]->getViscosity(_Udiff[_nVar-1]);

	if(lambda==0 && mu ==0 && _heatTransfertCoeff==0)
		return;

	//extraction des valeurs
	_idm[0] = _nVar*i; // Kieu
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Calcul diffusion: variables primitives maille " << i<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q] << endl;
		cout << endl;
	}

	if(!isBord ){
		for(int k=0; k<_nVar; k++)
			_idn[k] = _nVar*j + k;

		VecGetValues(_primitiveVars, _nVar, _idn, _Vj);
	}
	else
	{
		lambda=max(lambda,_heatTransfertCoeff);//wall nucleate boing -> larger heat transfer
		for(int k=0; k<_nVar; k++)
			_idn[k] = k;

		VecGetValues(_Uextdiff, _nVar, _idn, _phi);
		consToPrim(_phi,_Vj,1);
	}

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Calcul diffusion: variables primitives maille " <<j <<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q] << endl;
		cout << endl;
	}
	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_nVar-1; k++)
		_phi[k] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*mu*(_porosityj*_Vj[k] - _porosityi*_Vi[k]);

	//contribution visqueuse sur l'energie
	_phi[_nVar-1] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*lambda*(_porosityj*_Vj[_nVar-1] - _porosityi*_Vi[_nVar-1]);

	_idm[0] = i;
	VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion au 2nd membre pour la maille " << i << ": "<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _phi[q] << endl;
		cout << endl;
	}

	if(!isBord)
	{
		//On change de signe pour l'autre contribution
		for(int k=0; k<_nVar; k++)
			_phi[k] *= -_inv_dxj/_inv_dxi;
		_idn[0] = j;

		VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Contribution diffusion au 2nd membre pour la maille  " << j << ": "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}
	}

	if(_verbose && _nbTimeStep%_freqSave ==0 && _timeScheme==Implicit)
	{
		cout << "Matrice de diffusion D, pour le couple (" << i << "," << j<< "):" << endl;
		for(int i=0; i<_nVar; i++)
		{
			for(int j=0; j<_nVar; j++)
				cout << _Diffusion[i*_nVar+j]<<", ";
			cout << endl;
		}
		cout << endl;
	}
}

void SinglePhase::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{
	double phirho=Ui[0], T=Vi[_nVar-1];
	double norm_u=0;
	for(int k=0; k<_Ndim; k++)
		norm_u+=Vi[1+k]*Vi[1+k];
	norm_u=sqrt(norm_u);
	if(T>_Tsat)
		Si[0]=_heatPowerField(i)/_latentHeat;
	else
		Si[0]=0;
	for(int k=1; k<_nVar-1; k++)
		Si[k]  =(_gravite[k]-_dragCoeffs[0]*norm_u*Vi[1+k])*phirho;

	Si[_nVar-1]=_heatPowerField(i);

	for(int k=0; k<_Ndim; k++)
		Si[_nVar-1] +=(_GravityField3d[k]-_dragCoeffs[0]*norm_u*Vi[1+k])*Vi[1+k]*phirho;

	if(_timeScheme==Implicit)
	{
		for(int k=0; k<_nVar*_nVar;k++)
			_GravityImplicitationMatrix[k] = 0;
		if(!_usePrimitiveVarsInNewton)
			for(int k=0; k<_nVar;k++)
				_GravityImplicitationMatrix[k*_nVar]=-_gravite[k];
		else
		{
			double pression=Vi[0];
			getDensityDerivatives( pression, T, norm_u*norm_u);
			for(int k=0; k<_nVar;k++)
			{
				_GravityImplicitationMatrix[k*_nVar+0]      =-_gravite[k]*_drho_sur_dp;
				_GravityImplicitationMatrix[k*_nVar+_nVar-1]=-_gravite[k]*_drho_sur_dT;
			}
		}
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::sourceVector"<<endl;
		cout<<"Ui="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Ui[k]<<", ";
		cout<<endl;
		cout<<"Vi="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vi[k]<<", ";
		cout<<endl;
		cout<<"Si="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Si[k]<<", ";
		cout<<endl;
		if(_timeScheme==Implicit)
			displayMatrix(_GravityImplicitationMatrix, _nVar, "Gravity implicitation matrix");
	}
}

void SinglePhase::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{
	double norm_u=0, u_n=0, rho;
	for(int i=0;i<_Ndim;i++)
		u_n += _Uroe[1+i]*_vec_normal[i];

	pressureLoss[0]=0;
	if(u_n>0){
		for(int i=0;i<_Ndim;i++)
			norm_u += Vi[1+i]*Vi[1+i];
		norm_u=sqrt(norm_u);
		rho=Ui[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-1/2*K*rho*norm_u*Vi[1+i];
	}
	else{
		for(int i=0;i<_Ndim;i++)
			norm_u += Vj[1+i]*Vj[1+i];
		norm_u=sqrt(norm_u);
		rho=Uj[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-1/2*K*rho*norm_u*Vj[1+i];
	}
	pressureLoss[_nVar-1]=-1/2*K*rho*norm_u*norm_u*norm_u;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::pressureLossVector K= "<<K<<endl;
		cout<<"Ui="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Ui[k]<<", ";
		cout<<endl;
		cout<<"Vi="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vi[k]<<", ";
		cout<<endl;
		cout<<"Uj="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Uj[k]<<", ";
		cout<<endl;
		cout<<"Vj="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vj[k]<<", ";
		cout<<endl;
		cout<<"pressureLoss="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<pressureLoss[k]<<", ";
		cout<<endl;
	}
}

void SinglePhase::porosityGradientSourceVector()
{
	double u_ni=0, u_nj=0, rhoi,rhoj, pi=_Vi[0], pj=_Vj[0],pij;
	for(int i=0;i<_Ndim;i++) {
		u_ni += _Vi[1+i]*_vec_normal[i];
		u_nj += _Vj[1+i]*_vec_normal[i];
	}
	_porosityGradientSourceVector[0]=0;
	rhoj=_Uj[0]/_porosityj;
	rhoi=_Ui[0]/_porosityi;
	pij=(pi+pj)/2+rhoi*rhoj/2/(rhoi+rhoj)*(u_ni-u_nj)*(u_ni-u_nj);
	for(int i=0;i<_Ndim;i++)
		_porosityGradientSourceVector[1+i]=pij*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
	_porosityGradientSourceVector[_nVar-1]=0;
}

void SinglePhase::jacobian(const int &j, string nameOfGroup,double * normale)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite convection bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_Jcb[k] = 0;//No implicitation at this stage

	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);
	double q_n=0;//quantité de mouvement normale à la paroi
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+1)]*normale[k];

	// loop of boundary types
	if (_limitField[nameOfGroup].bcType==Wall)
	{
		for(k=0; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
		for(k=1; k<_nVar-1;k++)
			for(int l=1; l<_nVar-1;l++)
				_Jcb[k*_nVar + l] -= 2*normale[k-1]*normale[l-1];
	}
	else if (_limitField[nameOfGroup].bcType==Inlet)
	{
		return;
		if(q_n<0){
			double v[_Ndim], ve[_Ndim], v2, ve2;

			_idm[0] = _nVar*j;
			for(k=1; k<_nVar; k++)
				_idm[k] = _idm[k-1] + 1;
			VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
			VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

			ve[0] = _limitField[nameOfGroup].v_x[0];
			v[0]=_Vj[1];
			ve2 = ve[0]*ve[0];
			v2 = v[0]*v[0];
			if (_Ndim >1){
				ve[1] = _limitField[nameOfGroup].v_y[0];
				v[1]=_Vj[2];
				ve2 += ve[1]*ve[1];
				v2 = v[1]*v[1];
			}
			if (_Ndim >2){
				ve[2] = _limitField[nameOfGroup].v_z[0];
				v[2]=_Vj[3];
				ve2 += ve[2]*ve[2];
				v2 = v[2]*v[2];
			}
			double internal_energy=_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,_Uj[0]);
			double total_energy=internal_energy+ve2/2;

			//Mass line
			_Jcb[0]=v2/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[1+k]=-v[k]/internal_energy;
			_Jcb[_nVar-1]=1/internal_energy;
			//Momentum lines
			for(int l =1;l<1+_Ndim;l++){
				_Jcb[l*_nVar]=v2*ve[l-1]/(2*internal_energy);
				for(k=0; k<_Ndim;k++)
					_Jcb[l*_nVar+1+k]=-v[k]*ve[l-1]/internal_energy;
				_Jcb[l*_nVar+_nVar-1]=ve[l-1]/internal_energy;
			}
			//Energy line
			_Jcb[(_nVar-1)*_nVar]=v2*total_energy/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[(_nVar-1)*_nVar+1+k]=-v[k]*total_energy/internal_energy;
			_Jcb[(_nVar-1)*_nVar+_nVar-1]=total_energy/internal_energy;
		}
		else
			for(k=0;k<_nVar;k++)
				_Jcb[k*_nVar+k]=1;
		//Old jacobian
		/*
		 _Jcb[0] = 1;
		_Jcb[_nVar]=_limitField[nameOfGroup].v_x[0];//Kieu
		v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
		if(_Ndim>1)
		{
			_Jcb[2*_nVar]= _limitField[nameOfGroup].v_y[0];
			v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			if(_Ndim==3){
				_Jcb[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
			}
		}
		_Jcb[(_nVar-1)*_nVar]=_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2;
		 */
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure && q_n<0){
		return;
		double v[_Ndim], v2=0;
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

		for(k=0; k<_Ndim;k++){
			v[k]=_Vj[1+k];
			v2+=v[k]*v[k];
		}

		double rho_ext=_fluides[0]->getDensity(_limitField[nameOfGroup].p, _limitField[nameOfGroup].T);
		double rho_int = _externalStates[0];
		double density_ratio=rho_ext/rho_int;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_Jcb[l*_nVar]=-density_ratio*v[l-1];
			_Jcb[l*_nVar+l]=density_ratio;
		}
		//Energy lines
		_Jcb[(_nVar-1)*_nVar]=-v2*density_ratio;
		for(k=0; k<_Ndim;k++)
			_Jcb[(_nVar-1)*_nVar+1+k]=density_ratio*v[k];
	}
	// not wall, not inlet, not inletPressure
	else if(_limitField[nameOfGroup].bcType==Outlet || (_limitField[nameOfGroup].bcType==InletPressure && q_n>=0))
	{
		return;
		double v[_Ndim], v2=0;
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

		for(k=0; k<_Ndim;k++){
			v[k]=_Vj[1+k];
			v2+=v[k]*v[k];
		}

		double rho_ext=_fluides[0]->getDensity(_limitField[nameOfGroup].p, _externalStates[_nVar-1]);
		double rho_int = _externalStates[0];
		double density_ratio=rho_ext/rho_int;
		double internal_energy=_fluides[0]->getInternalEnergy(_externalStates[_nVar-1],rho_int);
		double total_energy=internal_energy+v2/2;

		//Mass line
		_Jcb[0]=density_ratio*(1-v2/(2*internal_energy));
		for(k=0; k<_Ndim;k++)
			_Jcb[1+k]=density_ratio*v[k]/internal_energy;
		_Jcb[_nVar-1]=-density_ratio/internal_energy;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_Jcb[l*_nVar]=density_ratio*v2*v[l-1]/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[l*_nVar+1+k]=density_ratio*v[k]*v[l-1]/internal_energy;
			_Jcb[l*_nVar+1+k]-=density_ratio;
			_Jcb[l*_nVar+_nVar-1]=-density_ratio*v[l-1]/internal_energy;
		}
		//Energy line
		_Jcb[(_nVar-1)*_nVar]=density_ratio*v2*total_energy/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_Jcb[(_nVar-1)*_nVar+1+k]=density_ratio*v[k]*total_energy/internal_energy;
		_Jcb[(_nVar-1)*_nVar+_nVar-1]=density_ratio*(1-total_energy/internal_energy);
		//Old jacobian
		/*
		int idim,jdim;
		double cd = 1,cn=0,p0, gamma;
		_idm[0] = j*_nVar;// Kieu
		for(k=1; k<_nVar;k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_conservativeVars, _nVar, _idm, _phi);
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		// compute the common numerator and common denominator
		p0=_fluides[0]->constante("p0");
		gamma =_fluides[0]->constante("gamma");
		cn =_limitField[nameOfGroup].p +p0;
		cd = _phi[0]*_fluides[0]->getInternalEnergy(_externalStates[_nVar-1],rho)-p0;
		cd*=cd;
		cd*=(gamma-1);
		//compute the v2
		for(k=1; k<_nVar-1;k++)
			v2+=_externalStates[k]*_externalStates[k];
		// drho_ext/dU
		_JcbDiff[0] = cn*(_phi[_nVar-1] -v2 -p0)/cd;
		for(k=1; k<_nVar-1;k++)
			_JcbDiff[k]=cn*_phi[k]/cd;
		_JcbDiff[_nVar-1]= -cn*_phi[0]/cd;
		//dq_ext/dU
		for(idim=0; idim<_Ndim;idim++)
		{
			//premiere colonne
			_JcbDiff[(1+idim)*_nVar]=-(v2*cn*_phi[idim+1])/(2*cd);
			//colonnes intermediaire
			for(jdim=0; jdim<_Ndim;jdim++)
			{
				_JcbDiff[(1+idim)*_nVar + jdim + 1] =_externalStates[idim+1]*_phi[jdim+1];
				_JcbDiff[(1+idim)*_nVar + jdim + 1]*=cn/cd;
			}
			//matrice identite*cn*(rhoe- p0)
			_JcbDiff[(1+idim)*_nVar + idim + 1] +=( cn*(_phi[0]*_fluides[0]->getInternalEnergy(_externalStates[_nVar-1],rho)-p0))/cd;

			//derniere colonne
			_JcbDiff[(1+idim)*_nVar + _nVar-1]=-_phi[idim+1]*cn/cd;
		}
		//drhoE/dU
		_JcbDiff[_nVar*(_nVar-1)] = -(v2*_phi[_nVar -1]*cn)/(2*cd);
		for(int idim=0; idim<_Ndim;idim++)
			_JcbDiff[_nVar*(_nVar-1)+idim+1]=_externalStates[idim+1]*_phi[_nVar -1]*cn/cd;
		_JcbDiff[_nVar*_nVar -1] = -(v2/2+p0)*cn/cd;
		 */
	}
	else  if (_limitField[nameOfGroup].bcType!=Neumann)// not wall, not inlet, not outlet
	{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::jacobian: This boundary condition is not treated");
	}
}

//calcule la jacobienne pour une CL de diffusion
void  SinglePhase::jacobianDiff(const int &j, string nameOfGroup)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite diffusion bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;

	if (_limitField[nameOfGroup].bcType==Wall){
		double v[_Ndim], ve[_Ndim], v2, ve2;

		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

		ve[0] = _limitField[nameOfGroup].v_x[0];
		v[0]=_Vj[1];
		ve2 = ve[0]*ve[0];
		v2 = v[0]*v[0];
		if (_Ndim >1){
			ve[1] = _limitField[nameOfGroup].v_y[0];
			v[1]=_Vj[2];
			ve2 += ve[1]*ve[1];
			v2 = v[1]*v[1];
		}
		if (_Ndim >2){
			ve[2] = _limitField[nameOfGroup].v_z[0];
			v[2]=_Vj[3];
			ve2 += ve[2]*ve[2];
			v2 = v[2]*v[2];
		}
		double rho=_Uj[0];
		double internal_energy=_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho);
		double total_energy=internal_energy+ve2/2;

		//Mass line
		_JcbDiff[0]=v2/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_JcbDiff[1+k]=-v[k]/internal_energy;
		_JcbDiff[_nVar-1]=1/internal_energy;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_JcbDiff[l*_nVar]=v2*ve[l-1]/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_JcbDiff[l*_nVar+1+k]=-v[k]*ve[l-1]/internal_energy;
			_JcbDiff[l*_nVar+_nVar-1]=ve[l-1]/internal_energy;
		}
		//Energy line
		_JcbDiff[(_nVar-1)*_nVar]=v2*total_energy/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_JcbDiff[(_nVar-1)*_nVar+1+k]=-v[k]*total_energy/internal_energy;
		_JcbDiff[(_nVar-1)*_nVar+_nVar-1]=total_energy/internal_energy;
	}
	else if (_limitField[nameOfGroup].bcType==Outlet || _limitField[nameOfGroup].bcType==Neumann
			||_limitField[nameOfGroup].bcType==Inlet || _limitField[nameOfGroup].bcType==InletPressure)
	{
		for(k=0;k<_nVar;k++)
			_JcbDiff[k*_nVar+k]=1;
	}
	else{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::jacobianDiff: This boundary condition is not recognised");
	}
}

void SinglePhase::primToCons(const double *P, const int &i, double *W, const int &j){
	//cout<<"SinglePhase::primToCons i="<<i<<", j="<<j<<", P[i*(_Ndim+2)]="<<P[i*(_Ndim+2)]<<", P[i*(_Ndim+2)+(_Ndim+1)]="<<P[i*(_Ndim+2)+(_Ndim+1)]<<endl;

	double rho=_fluides[0]->getDensity(P[i*(_Ndim+2)], P[i*(_Ndim+2)+(_Ndim+1)]);
	W[j*(_Ndim+2)] =  _porosityField(j)*rho;//phi*rho
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(k+1)] = W[j*(_Ndim+2)]*P[i*(_Ndim+2)+(k+1)];//phi*rho*u

	W[j*(_Ndim+2)+(_Ndim+1)] = W[j*(_Ndim+2)]*_fluides[0]->getInternalEnergy(P[i*(_Ndim+2)+ (_Ndim+1)],rho);//rho*e
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(_Ndim+1)] += W[j*(_Ndim+2)]*P[i*(_Ndim+2)+(k+1)]*P[i*(_Ndim+2)+(k+1)]*0.5;//phi*rho*e+0.5*phi*rho*u^2
}

void SinglePhase::primToConsJacobianMatrix(double *V)
{
	double pression=V[0];
	double temperature=V[_nVar-1];
	double vitesse[_Ndim];
	for(int idim=0;idim<_Ndim;idim++)
		vitesse[idim]=V[1+idim];
	double v2=0;
	for(int idim=0;idim<_Ndim;idim++)
		v2+=vitesse[idim]*vitesse[idim];

	double rho=_fluides[0]->getDensity(pression,temperature);
	double gamma=_fluides[0]->constante("gamma");
	double Pinf=_fluides[0]->constante("p0");
	double q=_fluides[0]->constante("q");
	double cv=_fluides[0]->constante("cv");

	for(int k=0;k<_nVar*_nVar; k++)
		_primToConsJacoMat[k]=0;

	if(		!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		double e=fluide0->getInternalEnergy(temperature);
		double E=e+0.5*v2;

		_primToConsJacoMat[0]=1/((gamma-1)*(e-q));
		_primToConsJacoMat[_nVar-1]=-rho*cv/(e-q);

		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[_nVar+idim*_nVar]=vitesse[idim]/((gamma-1)*(e-q));
			_primToConsJacoMat[_nVar+idim*_nVar+1+idim]=rho;
			_primToConsJacoMat[_nVar+idim*_nVar+_nVar-1]=-rho*vitesse[idim]*cv/(e-q);
		}
		_primToConsJacoMat[(_nVar-1)*_nVar]=E/((gamma-1)*(e-q));
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+1+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*cv*(1-E/(e-q));
	}
	else if(	_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		double h=fluide0->getEnthalpy(temperature);
		double H=h+0.5*v2;
		double cp=_fluides[0]->constante("cp");

		_primToConsJacoMat[0]=gamma/((gamma-1)*(h-q));
		_primToConsJacoMat[_nVar-1]=-rho*cp/(h-q);

		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[_nVar+idim*_nVar]=gamma*vitesse[idim]/((gamma-1)*(h-q));
			_primToConsJacoMat[_nVar+idim*_nVar+1+idim]=rho;
			_primToConsJacoMat[_nVar+idim*_nVar+_nVar-1]=-rho*vitesse[idim]*cp/(h-q);
		}
		_primToConsJacoMat[(_nVar-1)*_nVar]=gamma*H/((gamma-1)*(h-q))-1;
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+1+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*cp*(1-H/(h-q));
	}
	else
	{
		*_runLogFile<<"SinglePhase::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie"<<endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" SinglePhase::primToConsJacobianMatrix" << endl;
		displayVector(_Vi,_nVar," _Vi " );
		cout<<" Jacobienne primToCons: " << endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobienne primToCons: ");
	}
}

void SinglePhase::consToPrim(const double *Wcons, double* Wprim,double porosity)
{
	double q_2 = 0;
	for(int k=1;k<=_Ndim;k++)
		q_2 += Wcons[k]*Wcons[k];
	q_2 /= Wcons[0];	//phi rho u²
	double rhoe=(Wcons[(_Ndim+2)-1]-q_2/2)/porosity;
	double rho=Wcons[0]/porosity;
	Wprim[0] =  _fluides[0]->getPressure(rhoe,rho);//pressure p
	if (Wprim[0]<0){
		cout << "pressure = "<< Wprim[0] << " < 0 " << endl;
		*_runLogFile<< "pressure = "<< Wprim[0] << " < 0 " << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::consToPrim: negative pressure");
	}
	for(int k=1;k<=_Ndim;k++)
		Wprim[k] = Wcons[k]/Wcons[0];//velocity u
	Wprim[(_Ndim+2)-1] =  _fluides[0]->getTemperatureFromPressure(Wprim[0],Wcons[0]/porosity);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"ConsToPrim Vecteur conservatif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Wcons[k]<<endl;
		cout<<"ConsToPrim Vecteur primitif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Wprim[k]<<endl;
	}
}

void SinglePhase::entropicShift(double* n)//TO do: make sure _Vi and _Vj are well set
{
	/*Left and right values */
	double ul_n = 0, ul_2=0, ur_n=0,	ur_2 = 0; //valeurs de vitesse normale et |u|² a droite et a gauche
	for(int i=0;i<_Ndim;i++)
	{
		ul_n += _Vi[1+i]*n[i];
		ul_2 += _Vi[1+i]*_Vi[1+i];
		ur_n += _Vj[1+i]*n[i];
		ur_2 += _Vj[1+i]*_Vj[1+i];
	}


	double cl = _fluides[0]->vitesseSonEnthalpie(_Vi[_Ndim+1]-ul_2/2);//vitesse du son a l'interface
	double cr = _fluides[0]->vitesseSonEnthalpie(_Vj[_Ndim+1]-ur_2/2);//vitesse du son a l'interface

	_entropicShift[0]=fabs(ul_n-cl - (ur_n-cr));
	_entropicShift[1]=fabs(ul_n     - ur_n);
	_entropicShift[2]=fabs(ul_n+cl - (ur_n+cr));

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Entropic shift left eigenvalues: "<<endl;
		cout<<"("<< ul_n-cl<< ", " <<ul_n<<", "<<ul_n+cl << ")";
		cout<<endl;
		cout<<"Entropic shift right eigenvalues: "<<endl;
		cout<<"("<< ur_n-cr<< ", " <<ur_n<<", "<<ur_n+cr << ")";
		cout<< endl;
		cout<<"eigenvalue jumps "<<endl;
		cout<< _entropicShift[0] << ", " << _entropicShift[1] << ", "<< _entropicShift[2] <<endl;
	}
}

Vector SinglePhase::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::convectionFlux start"<<endl;
		cout<<"Ucons"<<endl;
		cout<<U<<endl;
		cout<<"Vprim"<<endl;
		cout<<V<<endl;
	}

	double phirho=U(0);//phi rho
	Vector phiq(_Ndim);//phi rho u
	for(int i=0;i<_Ndim;i++)
		phiq(i)=U(1+i);

	double pression=V(0);
	Vector vitesse(_Ndim);
	for(int i=0;i<_Ndim;i++)
		vitesse(i)=V(1+i);
	double Temperature= V(1+_Ndim);

	double vitessen=vitesse*normale;
	double rho=phirho/porosity;
	double e_int=_fluides[0]->getInternalEnergy(Temperature,rho);

	Vector F(_nVar);
	F(0)=phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F(1+i)=phirho*vitessen*vitesse(i)+pression*normale(i)*porosity;
	F(1+_Ndim)=phirho*(e_int+0.5*vitesse*vitesse+pression/rho)*vitessen;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void SinglePhase::RoeMatrixConservativeVariables(double u_n, double H,Vector velocity, double k, double K)
{
	/******** Construction de la matrice de Roe *********/
	//premiere ligne (masse)
	_Aroe[0]=0;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[1+idim]=_vec_normal[idim];
	_Aroe[_nVar-1]=0;

	//lignes intermadiaires(qdm)
	for(int idim=0; idim<_Ndim;idim++)
	{
		//premiere colonne
		_Aroe[(1+idim)*_nVar]=K*_vec_normal[idim] - u_n*_Uroe[1+idim];
		//colonnes intermediaires
		for(int jdim=0; jdim<_Ndim;jdim++)
			_Aroe[(1+idim)*_nVar + jdim + 1] = _Uroe[1+idim]*_vec_normal[jdim]-k*_vec_normal[idim]*_Uroe[1+jdim];
		//matrice identite
		_Aroe[(1+idim)*_nVar + idim + 1] += u_n;
		//derniere colonne
		_Aroe[(1+idim)*_nVar + _nVar-1]=k*_vec_normal[idim];
	}

	//derniere ligne (energie)
	_Aroe[_nVar*(_nVar-1)] = (K - H)*u_n;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] - k*u_n*_Uroe[idim+1];
	_Aroe[_nVar*_nVar -1] = (1 + k)*u_n;
}
void SinglePhase::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//On remplit la matrice de Roe en variables primitives : F(V_L)-F(V_R)=Aroe (V_L-V_R)
	//EOS is more involved with primitive variables
	// call to getDensityDerivatives(double concentration, double pression, double temperature,double v2) needed
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim]=_drho_sur_dT*u_n;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*vitesse[i]+_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*vitesse[i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim]=_drho_sur_dT*u_n*vitesse[i];
	}
	_AroeImplicit[(1+_Ndim)*_nVar+0]=(_drhoE_sur_dp+1)*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(1+_Ndim)*_nVar+1+i]=rho*(H*_vec_normal[i]+u_n*vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim]=_drhoE_sur_dT*u_n;
}
void SinglePhase::staggeredRoeUpwindingMatrixConservativeVariables( double u_n, double H,Vector velocity, double k, double K)
{
	//Calcul de décentrement de type décalé pour formulation de Roe
	if(fabs(u_n)>_precision)
	{
		//premiere ligne (masse)
		_absAroe[0]=0;
		for(int idim=0; idim<_Ndim;idim++)
			_absAroe[1+idim]=_vec_normal[idim];
		_absAroe[_nVar-1]=0;

		//lignes intermadiaires(qdm)
		for(int idim=0; idim<_Ndim;idim++)
		{
			//premiere colonne
			_absAroe[(1+idim)*_nVar]=-K*_vec_normal[idim] - u_n*_Uroe[1+idim];
			//colonnes intermediaires
			for(int jdim=0; jdim<_Ndim;jdim++)
				_absAroe[(1+idim)*_nVar + jdim + 1] = _Uroe[1+idim]*_vec_normal[jdim]+k*_vec_normal[idim]*_Uroe[1+jdim];
			//matrice identite
			_absAroe[(1+idim)*_nVar + idim + 1] += u_n;
			//derniere colonne
			_absAroe[(1+idim)*_nVar + _nVar-1]=-k*_vec_normal[idim];
		}

		//derniere ligne (energie)
		_absAroe[_nVar*(_nVar-1)] = (-K - H)*u_n;
		for(int idim=0; idim<_Ndim;idim++)
			_absAroe[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] + k*u_n*_Uroe[idim+1];
		_absAroe[_nVar*_nVar -1] = (1 - k)*u_n;

		double signu;
		if(u_n>0)
			signu=1;
		else if (u_n<0)
			signu=-1;

		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] *= signu;
	}
	else//umn=0 ->centered scheme
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] =0;
	}
}
void SinglePhase::staggeredRoeUpwindingMatrixPrimitiveVariables(double rho, double u_n,double H, Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//Calcul de décentrement de type décalé pour formulation Roe
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim]=_drho_sur_dT*u_n;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*vitesse[i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*vitesse[i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim]=_drho_sur_dT*u_n*vitesse[i];
	}
	_AroeImplicit[(1+_Ndim)*_nVar+0]=(_drhoE_sur_dp+1)*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(1+_Ndim)*_nVar+1+i]=rho*(H*_vec_normal[i]+u_n*vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim]=_drhoE_sur_dT*u_n;
}

Vector SinglePhase::staggeredVFFCFlux()
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhase::staggeredVFFCFlux() start"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<< "SinglePhase::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding, pressure = "<<  endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		Vector Fij(_nVar);

		double uijn=0, phiqn=0, uin=0, ujn=0;
		for(int idim=0;idim<_Ndim;idim++)
		{
			uijn+=_vec_normal[idim]*_Uroe[1+idim];//URoe = rho, u, H
			uin +=_vec_normal[idim]*_Ui[1+idim];
			ujn +=_vec_normal[idim]*_Uj[1+idim];
		}
		
		if( (uin>0 && ujn >0) || (uin>=0 && ujn <=0 && uijn>0) ) // formerly (uijn>_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Ui[1+idim];//phi rho u n
			Fij(0)=phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phiqn*_Vi[1+idim]+_Vj[0]*_vec_normal[idim]*_porosityj;
			Fij(_nVar-1)=phiqn/_Ui[0]*(_Ui[_nVar-1]+_Vj[0]*sqrt(_porosityj/_porosityi));
		}
		else if( (uin<0 && ujn <0) || (uin>=0 && ujn <=0 && uijn<0) ) // formerly (uijn<-_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Uj[1+idim];//phi rho u n
			Fij(0)=phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phiqn*_Vj[1+idim]+_Vi[0]*_vec_normal[idim]*_porosityi;
			Fij(_nVar-1)=phiqn/_Uj[0]*(_Uj[_nVar-1]+_Vi[0]*sqrt(_porosityi/_porosityj));
		}
		else//case (uin<=0 && ujn >=0) or (uin>=0 && ujn <=0 && uijn==0), apply centered scheme
		{
			Vector Ui(_nVar), Uj(_nVar), Vi(_nVar), Vj(_nVar), Fi(_nVar), Fj(_nVar);
			Vector normale(_Ndim);
			for(int i1=0;i1<_Ndim;i1++)
				normale(i1)=_vec_normal[i1];
			for(int i1=0;i1<_nVar;i1++)
			{
				Ui(i1)=_Ui[i1];
				Uj(i1)=_Uj[i1];
				Vi(i1)=_Vi[i1];
				Vj(i1)=_Vj[i1];
			}
			Fi=convectionFlux(Ui,Vi,normale,_porosityi);
			Fj=convectionFlux(Uj,Vj,normale,_porosityj);
			Fij=(Fi+Fj)/2;//+_maxvploc*(Ui-Uj)/2;
		}
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"SinglePhase::staggeredVFFCFlux() endf uijn="<<uijn<<endl;
			cout<<Fij<<endl;
		}
		return Fij;
	}
}

void SinglePhase::staggeredVFFCMatricesConservativeVariables(double un)//vitesse normale de Roe en entrée
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhase::staggeredVFFCMatrices()"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<< "SinglePhase::staggeredVFFCMatrices: staggeredVFFCMatrices method should be called only for VFFC formulation and staggered upwinding"<< endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::staggeredVFFCMatrices: staggeredVFFCMatrices method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		double ui_n=0, ui_2=0, uj_n=0, uj_2=0;//vitesse normale et carré du module
		double H;//enthalpie totale (expression particulière)
		consToPrim(_Ui,_Vi,_porosityi);
		consToPrim(_Uj,_Vj,_porosityj);
		for(int i=0;i<_Ndim;i++)
		{
			ui_2 += _Vi[1+i]*_Vi[1+i];
			ui_n += _Vi[1+i]*_vec_normal[i];
			uj_2 += _Vj[1+i]*_Vj[1+i];
			uj_n += _Vj[1+i]*_vec_normal[i];
		}

		double rhoi,pj, Ei, rhoj;
		double cj, Kj, kj;//dérivées de la pression
		rhoi=_Ui[0]/_porosityi;
		Ei= _Ui[_Ndim+1]/(rhoi*_porosityi);
		pj=_Vj[0];
		rhoj=_Uj[0]/_porosityj;
		H =Ei+pj/rhoi;
		cj = _fluides[0]->vitesseSonTemperature(_Vj[_Ndim+1],rhoj);
		kj = _fluides[0]->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
		Kj = uj_2*kj/2; //g-1/2 *|u|²

		double pi, Ej;
		double ci, Ki, ki;//dérivées de la pression
		Ej= _Uj[_Ndim+1]/rhoj;
		pi=_Vi[0];
		H =Ej+pi/rhoj;
		ci = _fluides[0]->vitesseSonTemperature(_Vi[_Ndim+1],rhoi);
		ki = _fluides[0]->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
		Ki = ui_2*ki/2; //g-1/2 *|u|²

		if(un>_precision)
		{
			/***********Calcul des valeurs propres ********/
			vector<complex<double> > vp_dist(3,0);
			vp_dist[0]=ui_n-cj;vp_dist[1]=ui_n;vp_dist[2]=ui_n+cj;

			_maxvploc=fabs(ui_n)+cj;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"Valeurs propres "<<ui_n-cj<<" , "<<ui_n<<" , "<<ui_n+cj<<endl;

			/******** Construction de la matrice A^+ *********/
			//premiere ligne (masse)
			_AroePlus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[1+idim]=_vec_normal[idim];
			_AroePlus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroePlus[(1+idim)*_nVar]=- ui_n*_Vi[1+idim];
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroePlus[(1+idim)*_nVar + jdim + 1] = _Vi[1+idim]*_vec_normal[jdim];
				//matrice identite
				_AroePlus[(1+idim)*_nVar + idim + 1] += ui_n;
				//derniere colonne
				_AroePlus[(1+idim)*_nVar + _nVar-1]=0;
			}

			//derniere ligne (energie)
			_AroePlus[_nVar*(_nVar-1)] = - H*ui_n;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] ;
			_AroePlus[_nVar*_nVar -1] = ui_n;

			/******** Construction de la matrice A^- *********/
			//premiere ligne (masse)
			_AroeMinus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[1+idim]=0;
			_AroeMinus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroeMinus[(1+idim)*_nVar]=Kj*_vec_normal[idim] ;
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroeMinus[(1+idim)*_nVar + jdim + 1] = -kj*_vec_normal[idim]*_Vj[1+jdim];
				//matrice identite
				_AroeMinus[(1+idim)*_nVar + idim + 1] += 0;
				//derniere colonne
				_AroeMinus[(1+idim)*_nVar + _nVar-1]=kj*_vec_normal[idim];
			}

			//derniere ligne (energie)
			_AroeMinus[_nVar*(_nVar-1)] = Kj *ui_n;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[_nVar*(_nVar-1)+idim+1]= - kj*uj_n*_Vi[idim+1];
			_AroeMinus[_nVar*_nVar -1] = kj*ui_n;
		}
		else if(un<-_precision)
		{
			/***********Calcul des valeurs propres ********/
			vector<complex<double> > vp_dist(3,0);
			vp_dist[0]=uj_n-ci;vp_dist[1]=uj_n;vp_dist[2]=uj_n+ci;

			_maxvploc=fabs(uj_n)+ci;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"Valeurs propres "<<uj_n-ci<<" , "<<uj_n<<" , "<<uj_n+ci<<endl;

			/******** Construction de la matrice A^+ *********/
			//premiere ligne (masse)
			_AroePlus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[1+idim]=0;
			_AroePlus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroePlus[(1+idim)*_nVar]=Ki*_vec_normal[idim] ;
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroePlus[(1+idim)*_nVar + jdim + 1] = -ki*_vec_normal[idim]*_Vi[1+jdim];
				//matrice identite
				_AroePlus[(1+idim)*_nVar + idim + 1] += 0;
				//derniere colonne
				_AroePlus[(1+idim)*_nVar + _nVar-1]=ki*_vec_normal[idim];
			}

			//derniere ligne (energie)
			_AroePlus[_nVar*(_nVar-1)] = Ki *uj_n;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[_nVar*(_nVar-1)+idim+1]= - ki*ui_n*_Vj[idim+1];
			_AroePlus[_nVar*_nVar -1] =  ki*uj_n;

			/******** Construction de la matrice A^- *********/
			//premiere ligne (masse)
			_AroeMinus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[1+idim]=_vec_normal[idim];
			_AroeMinus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroeMinus[(1+idim)*_nVar]= - uj_n*_Vj[1+idim];
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroeMinus[(1+idim)*_nVar + jdim + 1] = _Vj[1+idim]*_vec_normal[jdim];
				//matrice identite
				_AroeMinus[(1+idim)*_nVar + idim + 1] += uj_n;
				//derniere colonne
				_AroeMinus[(1+idim)*_nVar + _nVar-1]=0;
			}

			//derniere ligne (energie)
			_AroeMinus[_nVar*(_nVar-1)] = - H*uj_n;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] ;
			_AroeMinus[_nVar*_nVar -1] = uj_n;
		}
		else//case nil velocity on the interface, apply centered scheme
		{
			double u_n=0, u_2=0;//vitesse normale et carré du module
			for(int i=0;i<_Ndim;i++)
			{
				u_2 += _Uroe[1+i]*_Uroe[1+i];
				u_n += _Uroe[1+i]*_vec_normal[i];
			}
			Vector vitesse(_Ndim);
			for(int idim=0;idim<_Ndim;idim++)
				vitesse[idim]=_Uroe[1+idim];

			double  c, H, K, k;
			/***********Calcul des valeurs propres ********/
			H = _Uroe[_nVar-1];
			c = _fluides[0]->vitesseSonEnthalpie(H-u_2/2);//vitesse du son a l'interface
			k = _fluides[0]->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
			K = u_2*k/2; //g-1/2 *|u|²

			_maxvploc=fabs(u_n)+c;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			/******** Construction de la matrice A^+ *********/
			//premiere ligne (masse)
			_AroePlus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[1+idim]=0;
			_AroePlus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroePlus[(1+idim)*_nVar]=- ui_n*_Vi[1+idim];
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroePlus[(1+idim)*_nVar + jdim + 1] = _Vi[1+idim]*_vec_normal[jdim]-0.5*ki*_vec_normal[idim]*_Vi[1+jdim];
				//matrice identite
				_AroePlus[(1+idim)*_nVar + idim + 1] += 0.5*ui_n;
				//derniere colonne
				_AroePlus[(1+idim)*_nVar + _nVar-1]=0.5*ki*_vec_normal[idim];
			}

			//derniere ligne (energie)
			_AroePlus[_nVar*(_nVar-1)] = 0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroePlus[_nVar*(_nVar-1)+idim+1]=0 ;
			_AroePlus[_nVar*_nVar -1] = 0;

			/******** Construction de la matrice A^- *********/
			//premiere ligne (masse)
			_AroeMinus[0]=0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[1+idim]=0;
			_AroeMinus[_nVar-1]=0;

			//lignes intermadiaires(qdm)
			for(int idim=0; idim<_Ndim;idim++)
			{
				//premiere colonne
				_AroeMinus[(1+idim)*_nVar]=Kj*_vec_normal[idim] - uj_n*_Vj[1+idim];
				//colonnes intermediaires
				for(int jdim=0; jdim<_Ndim;jdim++)
					_AroeMinus[(1+idim)*_nVar + jdim + 1] = -0.5*kj*_vec_normal[idim]*_Vj[1+jdim];
				//matrice identite
				_AroeMinus[(1+idim)*_nVar + idim + 1] += 0.5*uj_n;
				//derniere colonne
				_AroeMinus[(1+idim)*_nVar + _nVar-1]=0.5*kj*_vec_normal[idim];
			}

			//derniere ligne (energie)
			_AroeMinus[_nVar*(_nVar-1)] = 0;
			for(int idim=0; idim<_Ndim;idim++)
				_AroeMinus[_nVar*(_nVar-1)+idim+1]= 0;
			_AroeMinus[_nVar*_nVar -1] = 0;
		}
	}
	if(_timeScheme==Implicit)
		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinusImplicit[i] = _AroeMinus[i];
			_AroePlusImplicit[i]  = _AroePlus[i];
		}

	/******** Construction de la matrices Aroe *********/
	/*
	//premiere ligne (masse)
	_Aroe[0]=0;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[1+idim]=_vec_normal[idim];
	_Aroe[_nVar-1]=0;

	//lignes intermadiaires(qdm)
	for(int idim=0; idim<_Ndim;idim++)
	{
		//premiere colonne
		_Aroe[(1+idim)*_nVar]=Ki*_vec_normal[idim] - uj_n*_Uj[1+idim];
		//colonnes intermediaires
		for(int jdim=0; jdim<_Ndim;jdim++)
			_Aroe[(1+idim)*_nVar + jdim + 1] = _Uj[1+idim]*_vec_normal[jdim]-ki*_vec_normal[idim]*_Ui[1+jdim];
		//matrice identite
		_Aroe[(1+idim)*_nVar + idim + 1] += uj_n;
		//derniere colonne
		_Aroe[(1+idim)*_nVar + _nVar-1]=ki*_vec_normal[idim];
	}

	//derniere ligne (energie)
	_Aroe[_nVar*(_nVar-1)] = (Ki - H)*uj_n;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] - ki*ui_n*_Uj[idim+1];
	_Aroe[_nVar*_nVar -1] = (1 + ki)*uj_n;
	 */
}

void SinglePhase::staggeredVFFCMatricesPrimitiveVariables(double un)//vitesse normale de Roe en entrée
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhase::staggeredVFFCMatricesPrimitiveVariables()"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<< "SinglePhase::staggeredVFFCMatricesPrimitiveVariables: staggeredVFFCMatricesPrimitiveVariables method should be called only for VFFC formulation and staggered upwinding" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::staggeredVFFCMatricesPrimitiveVariables: staggeredVFFCMatricesPrimitiveVariables method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		double ui_n=0., ui_2=0., uj_n=0., uj_2=0.;//vitesse normale et carré du module
		double H;//enthalpie totale (expression particulière)
		consToPrim(_Ui,_Vi,_porosityi);
		consToPrim(_Uj,_Vj,_porosityj);

		for(int i=0;i<_Ndim;i++)
		{
			ui_2 += _Vi[1+i]*_Vi[1+i];
			ui_n += _Vi[1+i]*_vec_normal[i];
			uj_2 += _Vj[1+i]*_Vj[1+i];
			uj_n += _Vj[1+i]*_vec_normal[i];
		}

		if(_verbose && _nbTimeStep%_freqSave ==0){
			cout <<"SinglePhase::staggeredVFFCMatricesPrimitiveVariables " << endl;
			cout<<"Vecteur primitif _Vi" << endl;
			for(int i=0;i<_nVar;i++)
				cout<<_Vi[i]<<", ";
			cout<<endl;
			cout<<"Vecteur primitif _Vj" << endl;
			for(int i=0;i<_nVar;i++)
				cout<<_Vj[i]<<", ";
			cout<<endl;
		}

		double gamma=_fluides[0]->constante("gamma");
		double q=_fluides[0]->constante("q");

		if(fabs(un)>_precision)//non zero velocity on the interface
		{
			if(	!_useDellacherieEOS)
			{
				StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
				double cv=fluide0->constante("cv");

				if(un>_precision)
				{
					double rhoi,rhoj,pj, Ei, ei;
					double cj;//vitesse du son pour calcul valeurs propres
					rhoi=_Ui[0]/_porosityi;
					Ei= _Ui[_Ndim+1]/(rhoi*_porosityi);
					ei=Ei-0.5*ui_2;
					pj=_Vj[0];
					rhoj=_Uj[0]/_porosityj;
					cj = _fluides[0]->vitesseSonTemperature(_Vj[_Ndim+1],rhoj);

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					vp_dist[0]=ui_n-cj;vp_dist[1]=ui_n;vp_dist[2]=ui_n+cj;

					_maxvploc=fabs(ui_n)+cj;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"Valeurs propres "<<ui_n-cj<<" , "<<ui_n<<" , "<<ui_n+cj<<endl;

					/******** Construction de la matrice A^+ *********/
					//premiere ligne (masse)
					_AroePlusImplicit[0]=ui_n/((gamma-1)*(ei-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[1+idim]=rhoi*_vec_normal[idim];
					_AroePlusImplicit[_nVar-1]=-rhoi*ui_n*cv/(ei-q);

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroePlusImplicit[(1+idim)*_nVar]=ui_n/((gamma-1)*(ei-q))*_Vi[1+idim];
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroePlusImplicit[(1+idim)*_nVar + jdim + 1] = rhoi*_Vi[1+idim]*_vec_normal[jdim];
						//matrice identite
						_AroePlusImplicit[(1+idim)*_nVar + idim + 1] += rhoi*ui_n;
						//derniere colonne
						_AroePlusImplicit[(1+idim)*_nVar + _nVar-1]=-rhoi*ui_n*cv/(ei-q)*_Vi[1+idim];
					}

					//derniere ligne (energie)
					_AroePlusImplicit[_nVar*(_nVar-1)] = Ei*ui_n/((gamma-1)*(ei-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[_nVar*(_nVar-1)+idim+1]=(rhoi*Ei+pj)*_vec_normal[idim]+rhoi*ui_n*_Vi[1+idim];
					_AroePlusImplicit[_nVar*_nVar -1] = rhoi*ui_n*(1-Ei/(ei-q))*cv;

					/******** Construction de la matrice A^- *********/
					//premiere ligne (masse)
					_AroeMinusImplicit[0]=0;
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[1+idim]=0;
					_AroeMinusImplicit[_nVar-1]=0;

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroeMinusImplicit[(1+idim)*_nVar]=_vec_normal[idim] ;
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroeMinusImplicit[(1+idim)*_nVar + jdim + 1] = 0;
						//matrice identite
						_AroeMinusImplicit[(1+idim)*_nVar + idim + 1] += 0;
						//derniere colonne
						_AroeMinusImplicit[(1+idim)*_nVar + _nVar-1]=0;
					}

					//derniere ligne (energie)
					_AroeMinusImplicit[_nVar*(_nVar-1)] = ui_n;
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[_nVar*(_nVar-1)+idim+1]= 0;
					_AroeMinusImplicit[_nVar*_nVar -1] = 0;
				}
				else if(un<-_precision)
				{
					double pi, Ej, rhoi, rhoj, ej;
					double ci;//vitesse du son pour calcul valeurs propres
					rhoj=_Uj[0]/_porosityj;
					Ej= _Uj[_Ndim+1]/rhoj;
					ej=Ej-0.5*uj_2;
					pi=_Vi[0];
					rhoi=_Ui[0]/_porosityi;
					ci = _fluides[0]->vitesseSonTemperature(_Vi[_Ndim+1],rhoi);

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					vp_dist[0]=uj_n-ci;vp_dist[1]=uj_n;vp_dist[2]=uj_n+ci;

					_maxvploc=fabs(uj_n)+ci;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"Valeurs propres "<<uj_n-ci<<" , "<<uj_n<<" , "<<uj_n+ci<<endl;

					/******** Construction de la matrice A^+ *********/
					//premiere ligne (masse)
					_AroePlusImplicit[0]=0;
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[1+idim]=0;
					_AroePlusImplicit[_nVar-1]=0;

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroePlusImplicit[(1+idim)*_nVar]=0;
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroePlusImplicit[(1+idim)*_nVar + jdim + 1] = 0;
						//matrice identite
						_AroePlusImplicit[(1+idim)*_nVar + idim + 1] += 0;
						//derniere colonne
						_AroePlusImplicit[(1+idim)*_nVar + _nVar-1]=0;
					}

					//derniere ligne (energie)
					_AroePlusImplicit[_nVar*(_nVar-1)] = uj_n;
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[_nVar*(_nVar-1)+idim+1]= 0;
					_AroePlusImplicit[_nVar*_nVar -1] =  0;

					/******** Construction de la matrice A^- *********/
					//premiere ligne (masse)
					_AroeMinusImplicit[0]=uj_n/((gamma-1)*(ej-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[1+idim]=rhoj*_vec_normal[idim];
					_AroeMinusImplicit[_nVar-1]=-rhoj*uj_n*cv/(ej-q);

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroeMinusImplicit[(1+idim)*_nVar]= uj_n/((gamma-1)*(ej-q))*_Vj[1+idim];
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroeMinusImplicit[(1+idim)*_nVar + jdim + 1] = rhoj*_Vj[1+idim]*_vec_normal[jdim];
						//matrice identite
						_AroeMinusImplicit[(1+idim)*_nVar + idim + 1] += rhoj*uj_n;
						//derniere colonne
						_AroeMinusImplicit[(1+idim)*_nVar + _nVar-1]=-rhoj*uj_n*cv/(ej-q)*_Vj[1+idim];
					}

					//derniere ligne (energie)
					_AroeMinusImplicit[_nVar*(_nVar-1)] = Ej*uj_n/((gamma-1)*(ej-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[_nVar*(_nVar-1)+idim+1]=(rhoj*Ej+pi)*_vec_normal[idim]+rhoj*uj_n*_Vj[1+idim];
					_AroeMinusImplicit[_nVar*_nVar -1] = rhoj*uj_n*(1-Ej/(ej-q))*cv;
				}
				else
				{
					*_runLogFile<< "SinglePhase::staggeredVFFCMatricesPrimitiveVariables: velocity un should be non zero" << endl;
					_runLogFile->close();
					throw CdmathException("SinglePhase::staggeredVFFCMatricesPrimitiveVariables: velocity un should be non zero");
				}
			}
			else if(_useDellacherieEOS )
			{
				StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
				double cp=fluide0->constante("cp");

				if(un>_precision)
				{
					double rhoi,rhoj,pj, Ei, hi, Hi;
					double cj;//vitesse du son pour calcul valeurs propres
					rhoi=_Ui[0]/_porosityi;
					Ei= _Ui[_Ndim+1]/(rhoi*_porosityi);
					Hi=Ei+_Vi[0]/rhoi;
					hi=Ei-0.5*ui_2;
					pj=_Vj[0];
					rhoj=_Uj[0]/_porosityj;
					cj = _fluides[0]->vitesseSonTemperature(_Vj[_Ndim+1],rhoj);

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					vp_dist[0]=ui_n-cj;vp_dist[1]=ui_n;vp_dist[2]=ui_n+cj;

					_maxvploc=fabs(ui_n)+cj;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"Valeurs propres "<<ui_n-cj<<" , "<<ui_n<<" , "<<ui_n+cj<<endl;

					/******** Construction de la matrice A^+ *********/
					//premiere ligne (masse)
					_AroePlusImplicit[0]=ui_n*gamma/((gamma-1)*(hi-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[1+idim]=rhoi*_vec_normal[idim];
					_AroePlusImplicit[_nVar-1]=-rhoi*ui_n*cp/(hi-q);

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroePlusImplicit[(1+idim)*_nVar]=ui_n*gamma/((gamma-1)*(hi-q))*_Vi[1+idim];
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroePlusImplicit[(1+idim)*_nVar + jdim + 1] = rhoi*_Vi[1+idim]*_vec_normal[jdim];
						//matrice identite
						_AroePlusImplicit[(1+idim)*_nVar + idim + 1] += rhoi*ui_n;
						//derniere colonne
						_AroePlusImplicit[(1+idim)*_nVar + _nVar-1]=-rhoi*ui_n*cp/(hi-q)*_Vi[1+idim];
					}

					//derniere ligne (energie)
					_AroePlusImplicit[_nVar*(_nVar-1)] = ui_n*(Hi*gamma/((gamma-1)*(hi-q))-1);
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[_nVar*(_nVar-1)+idim+1]=(rhoi*Ei+pj)*_vec_normal[idim]+rhoi*ui_n*_Vi[1+idim];
					_AroePlusImplicit[_nVar*_nVar -1] = rhoi*ui_n*(1-Hi/(hi-q))*cp;

					/******** Construction de la matrice A^- *********/
					//premiere ligne (masse)
					_AroeMinusImplicit[0]=0;
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[1+idim]=0;
					_AroeMinusImplicit[_nVar-1]=0;

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroeMinusImplicit[(1+idim)*_nVar]=_vec_normal[idim] ;
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroeMinusImplicit[(1+idim)*_nVar + jdim + 1] = 0;
						//matrice identite
						_AroeMinusImplicit[(1+idim)*_nVar + idim + 1] += 0;
						//derniere colonne
						_AroeMinusImplicit[(1+idim)*_nVar + _nVar-1]=0;
					}

					//derniere ligne (energie)
					_AroeMinusImplicit[_nVar*(_nVar-1)] = ui_n;
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[_nVar*(_nVar-1)+idim+1]= 0;
					_AroeMinusImplicit[_nVar*_nVar -1] = 0;
				}
				else if(un<-_precision)
				{
					double pi, Ej, rhoi,rhoj, Hj, hj;
					double ci;//vitesse du son pour calcul valeurs propres
					rhoj=_Uj[0]/_porosityj;
					Ej= _Uj[_Ndim+1]/rhoj;
					Hj=Ej+_Vj[0]/rhoj;
					hj=Ej-0.5*uj_2;
					pi=_Vi[0];
					rhoi=_Ui[0]/_porosityi;
					ci = _fluides[0]->vitesseSonTemperature(_Vi[_Ndim+1],rhoi);

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					vp_dist[0]=uj_n-ci;vp_dist[1]=uj_n;vp_dist[2]=uj_n+ci;

					_maxvploc=fabs(uj_n)+ci;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"Valeurs propres "<<uj_n-ci<<" , "<<uj_n<<" , "<<uj_n+ci<<endl;

					/******** Construction de la matrice A^+ *********/
					//premiere ligne (masse)
					_AroePlusImplicit[0]=0;
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[1+idim]=0;
					_AroePlusImplicit[_nVar-1]=0;

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroePlusImplicit[(1+idim)*_nVar]=0;
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroePlusImplicit[(1+idim)*_nVar + jdim + 1] = 0;
						//matrice identite
						_AroePlusImplicit[(1+idim)*_nVar + idim + 1] += 0;
						//derniere colonne
						_AroePlusImplicit[(1+idim)*_nVar + _nVar-1]=0;
					}

					//derniere ligne (energie)
					_AroePlusImplicit[_nVar*(_nVar-1)] = uj_n;
					for(int idim=0; idim<_Ndim;idim++)
						_AroePlusImplicit[_nVar*(_nVar-1)+idim+1]= 0;
					_AroePlusImplicit[_nVar*_nVar -1] =  0;

					/******** Construction de la matrice A^- *********/
					//premiere ligne (masse)
					_AroeMinusImplicit[0]=uj_n*gamma/((gamma-1)*(hj-q));
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[1+idim]=rhoj*_vec_normal[idim];
					_AroeMinusImplicit[_nVar-1]=-rhoj*uj_n*cp/(hj-q);

					//lignes intermadiaires(qdm)
					for(int idim=0; idim<_Ndim;idim++)
					{
						//premiere colonne
						_AroeMinusImplicit[(1+idim)*_nVar]= uj_n*gamma/((gamma-1)*(hj-q))*_Vj[1+idim];
						//colonnes intermediaires
						for(int jdim=0; jdim<_Ndim;jdim++)
							_AroeMinusImplicit[(1+idim)*_nVar + jdim + 1] = rhoj*_Vj[1+idim]*_vec_normal[jdim];
						//matrice identite
						_AroeMinusImplicit[(1+idim)*_nVar + idim + 1] += rhoj*uj_n;
						//derniere colonne
						_AroeMinusImplicit[(1+idim)*_nVar + _nVar-1]=-rhoj*uj_n*cp/(hj-q)*_Vj[1+idim];
					}

					//derniere ligne (energie)
					_AroeMinusImplicit[_nVar*(_nVar-1)] = uj_n*(Hj*gamma/((gamma-1)*(hj-q))-1);
					for(int idim=0; idim<_Ndim;idim++)
						_AroeMinusImplicit[_nVar*(_nVar-1)+idim+1]=(rhoj*Ej+pi)*_vec_normal[idim]+rhoj*uj_n*_Vj[1+idim];
					_AroeMinusImplicit[_nVar*_nVar -1] = rhoj*uj_n*(1-Hj/(hj-q))*cp;
				}
				else
				{
					*_runLogFile<< "SinglePhase::staggeredVFFCMatricesPrimitiveVariables: velocity un should be non zero" << endl;
					_runLogFile->close();
					throw CdmathException("SinglePhase::staggeredVFFCMatricesPrimitiveVariables: velocity un should be non zero");
				}
			}
			else
			{
				*_runLogFile<< "SinglePhase::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie" << endl;
				_runLogFile->close();
				throw CdmathException("SinglePhase::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie");
			}
		}
		else//case nil velocity on the interface, apply centered scheme
		{
			Polynoms Poly;
			primToConsJacobianMatrix(_Vj);
			Poly.matrixProduct(_AroeMinus, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroeMinusImplicit);
			primToConsJacobianMatrix(_Vi);
			Poly.matrixProduct(_AroePlus,  _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroePlusImplicit);
		}
	}
}
void SinglePhase::applyVFRoeLowMachCorrections(bool isBord, string groupname)
{
	if(_nonLinearFormulation!=VFRoe)
	{
		*_runLogFile<< "SinglePhase::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation");
	}
	else//_nonLinearFormulation==VFRoe
	{
		if(_spaceScheme==lowMach){
			double u_2=0;
			for(int i=0;i<_Ndim;i++)
				u_2 += _Uroe[1+i]*_Uroe[1+i];
			double 	c = _maxvploc;//vitesse du son a l'interface
			double M=sqrt(u_2)/c;//Mach number
			_Vij[0]=M*_Vij[0]+(1-M)*(_Vi[0]+_Vj[0])/2;
		}
		else if(_spaceScheme==pressureCorrection)
		{//order 1 : no correction, oarder 2 : correction everywhere, order 3 : correction only inside, orders 4 and 5 : special correction at boundaries
			if(_pressureCorrectionOrder==2 || (!isBord && _pressureCorrectionOrder==3) || (!isBord && _pressureCorrectionOrder==4) || (!isBord && _pressureCorrectionOrder==5) )
			{
				double norm_uij=0, uij_n=0, ui_n=0, uj_n=0;
				for(int i=0;i<_Ndim;i++)
				{
					norm_uij += _Uroe[1+i]*_Uroe[1+i];
					uij_n += _Uroe[1+i]*_vec_normal[i];
					ui_n += _Vi[1+i]*_vec_normal[i];
					uj_n += _Vj[1+i]*_vec_normal[i];
				}
				norm_uij=sqrt(norm_uij);
				if(norm_uij>_precision)//avoid division by zero
					_Vij[0]=(_Vi[0]+_Vj[0])/2 + uij_n/norm_uij*(_Vj[0]-_Vi[0])/4 - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
				else
					_Vij[0]=(_Vi[0]+_Vj[0])/2                                    - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
			}
			else if(_pressureCorrectionOrder==4 && isBord)
				_Vij[0]=_Vi[0];
			else if(_pressureCorrectionOrder==5 && isBord)
			{
				double g_n=0;//scalar product of gravity and normal vector
				for(int i=0;i<_Ndim;i++)
					g_n += _GravityField3d[i]*_vec_normal[i];
				_Vij[0]=_Vi[0]- _Ui[0]*g_n/_inv_dxi/2;
			}
		}
		else if(_spaceScheme==staggered)
		{
			double uij_n=0;
			for(int i=0;i<_Ndim;i++)
				uij_n += _Uroe[1+i]*_vec_normal[i];

			if(uij_n>_precision){
				_Vij[0]=_Vj[0];
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=_Vi[1+i];
				_Vij[_nVar-1]=_Vi[_nVar-1];
			}
			else if(uij_n<-_precision){
				_Vij[0]=_Vi[0];
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=_Vj[1+i];
				_Vij[_nVar-1]=_Vj[_nVar-1];
			}
			else{
				_Vij[0]=(_Vi[0]+_Vi[0])/2;
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=(_Vj[1+i]+_Vj[1+i])/2;
				_Vij[_nVar-1]=(_Vj[_nVar-1]+_Vj[_nVar-1])/2;
			}
		}
		primToCons(_Vij,0,_Uij,0);
	}
}

void SinglePhase::testConservation()
{
	double SUM, DELTA, x;
	int I;
	for(int i=0; i<_nVar; i++)
	{
		{
			if(i == 0)
				cout << "Masse totale (kg): ";
			else
			{
				if(i == _nVar-1)
					cout << "Energie totale (J): ";
				else
					cout << "Quantite de mouvement totale (kg.m.s^-1): ";
			}
		}
		SUM = 0;
		I =  i;
		DELTA = 0;
		for(int j=0; j<_Nmailles; j++)
		{
			if(!_usePrimitiveVarsInNewton)
				VecGetValues(_conservativeVars, 1, &I, &x);//on recupere la valeur du champ
			else
				VecGetValues(_primitiveVars, 1, &I, &x);//on recupere la valeur du champ
			SUM += x*_mesh.getCell(j).getMeasure();
			VecGetValues(_newtonVariation, 1, &I, &x);//on recupere la variation du champ
			DELTA += x*_mesh.getCell(j).getMeasure();
			I += _nVar;
		}
		if(fabs(SUM)>_precision)
			cout << SUM << ", variation relative: " << fabs(DELTA /SUM)  << endl;
		else
			cout << " a une somme quasi nulle,  variation absolue: " << fabs(DELTA) << endl;
	}
}

void SinglePhase::getDensityDerivatives( double pressure, double temperature, double v2)
{
	double rho=_fluides[0]->getDensity(pressure,temperature);
	double gamma=_fluides[0]->constante("gamma");
	double q=_fluides[0]->constante("q");

	if(	!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		double e = fluide0->getInternalEnergy(temperature);
		double cv=fluide0->constante("cv");
		double E=e+0.5*v2;

		_drho_sur_dp=1/((gamma-1)*(e-q));
		_drho_sur_dT=-rho*cv/(e-q);
		_drhoE_sur_dp=E/((gamma-1)*(e-q));
		_drhoE_sur_dT=rho*cv*(1-E/(e-q));
	}
	else if(_useDellacherieEOS )
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		double h=fluide0->getEnthalpy(temperature);
		double H=h+0.5*v2;
		double cp=fluide0->constante("cp");

		_drho_sur_dp=gamma/((gamma-1)*(h-q));
		_drho_sur_dT=-rho*cp/(h-q);
		_drhoE_sur_dp=gamma*H/((gamma-1)*(h-q))-1;
		_drhoE_sur_dT=rho*cp*(1-H/(h-q));
	}
	else
	{
		*_runLogFile<< "SinglePhase::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"_drho_sur_dp= "<<_drho_sur_dp<<", _drho_sur_dT= "<<_drho_sur_dT<<endl;
		cout<<"_drhoE_sur_dp= "<<_drhoE_sur_dp<<", _drhoE_sur_dT= "<<_drhoE_sur_dT<<endl;
	}
}
void SinglePhase::save(){
	string prim(_path+"/SinglePhasePrim_");///Results
	string cons(_path+"/SinglePhaseCons_");
	string allFields(_path+"/");
	prim+=_fileName;
	cons+=_fileName;
	allFields+=_fileName;

	PetscInt Ii;
	for (long i = 0; i < _Nmailles; i++){
		// j = 0 : pressure; j = _nVar - 1: temperature; j = 1,..,_nVar-2: velocity
		for (int j = 0; j < _nVar; j++){
			Ii = i*_nVar +j;
			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
		}
	}
	if(_saveConservativeField){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : density; j = _nVar - 1 : energy j = 1,..,_nVar-2: momentum
			for (int j = 0; j < _nVar; j++){
				Ii = i*_nVar +j;
				VecGetValues(_conservativeVars,1,&Ii,&_UU(i,j));
			}
		}
		_UU.setTime(_time,_nbTimeStep);
	}
	_VV.setTime(_time,_nbTimeStep);

	// create mesh and component info
	if (_nbTimeStep ==0 || _restartWithNewFileName){		
		string prim_suppress ="rm -rf "+prim+"_*";
		string cons_suppress ="rm -rf "+cons+"_*";

		system(prim_suppress.c_str());//Nettoyage des précédents calculs identiques
		system(cons_suppress.c_str());//Nettoyage des précédents calculs identiques

		if(_saveConservativeField){
			_UU.setInfoOnComponent(0,"Density_(kg/m^3)");
			_UU.setInfoOnComponent(1,"Momentum_x");// (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(2,"Momentum_y");// (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(3,"Momentum_z");// (kg/m^2/s)

			_UU.setInfoOnComponent(_nVar-1,"Energy_(J/m^3)");

			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons);
				break;
			case MED :
				_UU.writeMED(cons);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}
		_VV.setInfoOnComponent(0,"Pressure_(Pa)");
		_VV.setInfoOnComponent(1,"Velocity_x_(m/s)");
		if (_Ndim>1)
			_VV.setInfoOnComponent(2,"Velocity_y_(m/s)");
		if (_Ndim>2)
			_VV.setInfoOnComponent(3,"Velocity_z_(m/s)");
		_VV.setInfoOnComponent(_nVar-1,"Temperature_(K)");

		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim);
			break;
		case MED :
			_VV.writeMED(prim);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}

	}
	// do not create mesh
	else{
		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim,false);
			break;
		case MED :
			_VV.writeMED(prim,false);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}
		if(_saveConservativeField){
			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons,false);
				break;
			case MED :
				_UU.writeMED(cons,false);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}
	}
	if(_saveVelocity || _saveAllFields){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : pressure; j = _nVar - 1: temperature; j = 1,..,_nVar-2: velocity
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +1+j;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse(i,j));
			}
			for (int j = _Ndim; j < 3; j++)//On met à zero les composantes de vitesse si la dimension est <3
				_Vitesse(i,j)=0;
		}
		_Vitesse.setTime(_time,_nbTimeStep);
		if (_nbTimeStep ==0 || _restartWithNewFileName){		
			_Vitesse.setInfoOnComponent(0,"Velocity_x_(m/s)");
			_Vitesse.setInfoOnComponent(1,"Velocity_y_(m/s)");
			_Vitesse.setInfoOnComponent(2,"Velocity_z_(m/s)");

			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity",false);
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity",false);
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
	}

	if(_saveAllFields)
	{
		double p,T,rho, h, vx,vy,vz;
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&rho);
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&T);
			Ii = i*_nVar + 1;
			VecGetValues(_primitiveVars,1,&Ii,&vx);
			if(_Ndim>1)
			{
				Ii = i*_nVar + 2;
				VecGetValues(_primitiveVars,1,&Ii,&vy);
				if(_Ndim>2){
					Ii = i*_nVar + 3;
					VecGetValues(_primitiveVars,1,&Ii,&vz);
				}
			}

			h   = _fluides[0]->getEnthalpy(T,rho);

			_Enthalpy(i)=h;
			_Density(i)=rho;
			_Pressure(i)=p;
			_Temperature(i)=T;
			_VitesseX(i)=vx;
			if(_Ndim>1)
			{
				_VitesseY(i)=vy;
				if(_Ndim>2)
					_VitesseZ(i)=vz;
			}
		}
		_Enthalpy.setTime(_time,_nbTimeStep);
		_Density.setTime(_time,_nbTimeStep);
		_Pressure.setTime(_time,_nbTimeStep);
		_Temperature.setTime(_time,_nbTimeStep);
		_VitesseX.setTime(_time,_nbTimeStep);
		if(_Ndim>1)
		{
			_VitesseY.setTime(_time,_nbTimeStep);
			if(_Ndim>2)
				_VitesseZ.setTime(_time,_nbTimeStep);
		}
		if (_nbTimeStep ==0 || _restartWithNewFileName){		
			switch(_saveFormat)
			{
			case VTK :
				_Enthalpy.writeVTK(allFields+"_Enthalpy");
				_Density.writeVTK(allFields+"_Density");
				_Pressure.writeVTK(allFields+"_Pressure");
				_Temperature.writeVTK(allFields+"_Temperature");
				_VitesseX.writeVTK(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ");
				}
				break;
			case MED :
				_Enthalpy.writeMED(allFields+"_Enthalpy");
				_Density.writeMED(allFields+"_Density");
				_Pressure.writeMED(allFields+"_Pressure");
				_Temperature.writeMED(allFields+"_Temperature");
				_VitesseX.writeMED(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ");
				}
				break;
			case CSV :
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
				_VitesseX.writeCSV(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeCSV(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeCSV(allFields+"_VelocityZ");
				}
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Enthalpy.writeVTK(allFields+"_Enthalpy",false);
				_Density.writeVTK(allFields+"_Density",false);
				_Pressure.writeVTK(allFields+"_Pressure",false);
				_Temperature.writeVTK(allFields+"_Temperature",false);
				_VitesseX.writeVTK(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ",false);
				}
				break;
			case MED :
				_Enthalpy.writeMED(allFields+"_Enthalpy",false);
				_Density.writeMED(allFields+"_Density",false);
				_Pressure.writeMED(allFields+"_Pressure",false);
				_Temperature.writeMED(allFields+"_Temperature",false);
				_VitesseX.writeMED(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ",false);
				}
				break;
			case CSV :
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
				_VitesseX.writeCSV(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeCSV(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeCSV(allFields+"_VelocityZ");
				}
				break;
			}
		}
	}

	if(_isStationary)
	{
		prim+="_Stat";
		cons+="_Stat";

		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim);
			break;
		case MED :
			_VV.writeMED(prim);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}

		if(_saveConservativeField){
			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons);
				break;
			case MED :
				_UU.writeMED(cons);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}

		if(_saveVelocity || _saveAllFields){
			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
	}

	if (_restartWithNewFileName)
		_restartWithNewFileName=false;
}

Field& SinglePhase::getPressureField()
{
	if(!_saveAllFields)
	{
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&_Pressure(i));
		}
		_Pressure.setTime(_time,_nbTimeStep);
	}
	return _Pressure;
}

Field& SinglePhase::getTemperatureField()
{
	if(!_saveAllFields)
	{
		_Temperature=Field("Temperature",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&_Temperature(i));
		}
		_Temperature.setTime(_time,_nbTimeStep);
	}
	return _Temperature;
}

Field& SinglePhase::getVelocityField()
{
	if(!_saveAllFields )
	{
		_Vitesse=Field("Vitesse",CELLS,_mesh,3);
		int Ii;
		for (long i = 0; i < _Nmailles; i++)
		{
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +1+j;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse(i,j));
			}
			for (int j = _Ndim; j < 3; j++)//On met à zero les composantes de vitesse si la dimension est <3
				_Vitesse(i,j)=0;
		}
		_Vitesse.setTime(_time,_nbTimeStep);
		_Vitesse.setInfoOnComponent(0,"Velocity_x_(m/s)");
		_Vitesse.setInfoOnComponent(1,"Velocity_y_(m/s)");
		_Vitesse.setInfoOnComponent(2,"Velocity_z_(m/s)");
	}
	
	return _Vitesse;
}

Field& SinglePhase::getVelocityXField()
{
	if(!_saveAllFields )
	{
		_VitesseX=Field("Velocity X",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++)
		{
			int Ii = i*_nVar +1;
			VecGetValues(_primitiveVars,1,&Ii,&_VitesseX(i));
		}
		_VitesseX.setTime(_time,_nbTimeStep);
		_VitesseX.setInfoOnComponent(0,"Velocity_x_(m/s)");
	}
	
	return _VitesseX;
}

Field& SinglePhase::getDensityField()
{
	if(!_saveAllFields )
	{
		_Density=Field("Density",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&_Density(i));
		}
		_Density.setTime(_time,_nbTimeStep);
	}
	return _Density;
}

Field& SinglePhase::getMomentumField()//not yet managed by parameter _saveAllFields
{
	_Momentum=Field("Momentum",CELLS,_mesh,_Ndim);
	int Ii;
	for (long i = 0; i < _Nmailles; i++)
		for (int j = 0; j < _Ndim; j++)//On récupère les composantes de qdm
		{
			int Ii = i*_nVar +1+j;
			VecGetValues(_conservativeVars,1,&Ii,&_Momentum(i,j));
		}
	_Momentum.setTime(_time,_nbTimeStep);

	return _Momentum;
}

Field& SinglePhase::getTotalEnergyField()//not yet managed by parameter _saveAllFields
{
	_TotalEnergy=Field("TotalEnergy",CELLS,_mesh,1);
	int Ii;
	for (long i = 0; i < _Nmailles; i++){
		Ii = i*_nVar +_nVar-1;
		VecGetValues(_conservativeVars,1,&Ii,&_TotalEnergy(i));
	}
	_TotalEnergy.setTime(_time,_nbTimeStep);

	return _TotalEnergy;
}

Field& SinglePhase::getEnthalpyField()
{
	if(!_saveAllFields )
	{
		_Enthalpy=Field("Enthalpy",CELLS,_mesh,1);
		int Ii;
		double p,T,rho;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&T);
			
			rho=_fluides[0]->getDensity(p,T);
			_Enthalpy(i)=_fluides[0]->getEnthalpy(T,rho);
		}
		_Enthalpy.setTime(_time,_nbTimeStep);
	}

	return _Enthalpy;
}

vector<string> SinglePhase::getOutputFieldsNames()
{
	vector<string> result(8);
	
	result[0]="Pressure";
	result[1]="Velocity";
	result[2]="Temperature";
	result[3]="Density";
	result[4]="Momentum";
	result[5]="TotalEnergy";
	result[6]="Enthalpy";
	result[7]="VelocityX";
	
	return result;
}

Field& SinglePhase::getOutputField(const string& nameField )
{
	if(nameField=="pressure" || nameField=="Pressure" || nameField=="PRESSURE" || nameField=="PRESSION" || nameField=="Pression"  || nameField=="pression" )
		return getPressureField();
	else if(nameField=="velocity" || nameField=="Velocity" || nameField=="VELOCITY" || nameField=="Vitesse" || nameField=="VITESSE" || nameField=="vitesse" )
		return getVelocityField();
	else if(nameField=="velocityX" || nameField=="VelocityX" || nameField=="VELOCITYX" || nameField=="VitesseX" || nameField=="VITESSEX" || nameField=="vitesseX" )
		return getVelocityXField();
	else if(nameField=="temperature" || nameField=="Temperature" || nameField=="TEMPERATURE" || nameField=="temperature" )
		return getTemperatureField();
	else if(nameField=="density" || nameField=="Density" || nameField=="DENSITY" || nameField=="Densite" || nameField=="DENSITE" || nameField=="densite" )
		return getDensityField();
	else if(nameField=="momentum" || nameField=="Momentum" || nameField=="MOMENTUM" || nameField=="Qdm" || nameField=="QDM" || nameField=="qdm" )
		return getMomentumField();
	else if(nameField=="enthalpy" || nameField=="Enthalpy" || nameField=="ENTHALPY" || nameField=="Enthalpie" || nameField=="ENTHALPIE" || nameField=="enthalpie" )
		return getEnthalpyField();
	else if(nameField=="totalenergy" || nameField=="TotalEnergy" || nameField=="TOTALENERGY" || nameField=="ENERGIETOTALE" || nameField=="EnergieTotale" || nameField=="energietotale" )
		return getTotalEnergyField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first" << endl;
        throw CdmathException("SinglePhase::getOutputField error : Unknown Field name");
    }
}
