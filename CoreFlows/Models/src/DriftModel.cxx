/*
 * DriftModel.cxx
 *
 *  Created on: 1 janv. 2015
 *      Author: Michael Ndjinga
 */
#include "DriftModel.hxx"

using namespace std;

DriftModel::DriftModel(pressureEstimate pEstimate, int dim, bool useDellacherieEOS){
	_Ndim=dim;
	_nVar=_Ndim+3;
	_nbPhases  = 2;
	_dragCoeffs=vector<double>(2,0);
	_fluides.resize(2);
	_saveAllFields=false;

	if( pEstimate==around1bar300K){//EOS at 1 bar and 373K
		cout<<"Fluid is water-Gas mixture around saturation point 1 bar and 373 K (100°C)"<<endl;
		*_runLogFile<<"Fluid is water-Gas mixture around saturation point 1 bar and 373 K (100°C)"<<endl;
		_Tsat=373;//saturation temperature at 1 bar
		double esatv=2.505e6;//Gas internal energy at saturation at 1 bar
		double esatl=4.174e5;//water internal energy at saturation at 1 bar
		double cv_v=1555;//Gas specific heat capacity at saturation at 1 bar
		double cv_l=3770;//water specific heat capacity at saturation at 1 bar
		double gamma_v=1.337;//Gas heat capacity ratio at saturation at 1 bar
		double rho_sat_l=958;//water density at saturation at 1 bar
		double sound_speed_l=1543;//water sound speed at saturation at 1 bar
		_fluides[0] = new StiffenedGas(gamma_v,cv_v,_Tsat,esatv);  //ideal gas law for Gas at pressure 1 bar and temperature 100°C, gamma=1.34
		_fluides[1] = new StiffenedGas(rho_sat_l,1e5,_Tsat,esatl,sound_speed_l,cv_l);  //stiffened gas law for water at pressure 1 bar and temperature 100°C
		_hsatl=4.175e5;//water enthalpy at saturation at 1 bar
		_hsatv=2.675e6;//Gas enthalpy at saturation at 1 bar

		_useDellacherieEOS=false;
	}
	else{//EOS at 155 bars and 618K
		cout<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		*_runLogFile<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		_useDellacherieEOS=useDellacherieEOS;
		if(useDellacherieEOS)
		{
			_Tsat=656;//saturation temperature used in Dellacherie EOS
			_hsatl=1.633e6;//water enthalpy at saturation at 155 bars
			_hsatv=3.006e6;//Gas enthalpy at saturation at 155 bars
			_fluides[0] = new StiffenedGasDellacherie(1.43,0  ,2.030255e6  ,1040.14); //stiffened gas law for Gas from S. Dellacherie
			_fluides[1] = new StiffenedGasDellacherie(2.35,1e9,-1.167056e6,1816.2); //stiffened gas law for water from S. Dellacherie
		}
		else
		{
			double esatv=2.444e6;//Gas internal energy at saturation at 155 bar
			double esatl=1.604e6;//water internal energy at saturation at 155 bar
			double sound_speed_v=433;//Gas sound speed at saturation at 155 bar
			double sound_speed_l=621;//water sound speed at saturation at 155 bar
			double cv_v=3633;//Gas specific heat capacity at saturation at 155 bar
			double cv_l=3100;//water specific heat capacity at saturation at 155 bar
			double rho_sat_v=102;//Gas density at saturation at 155 bar
			double rho_sat_l=594;//water density at saturation at 155 bar
			_Tsat=618;//saturation temperature at 155 bars
			_hsatl=1.63e6;//water enthalpy at saturation at 155 bars
			_hsatv=2.6e6;//Gas enthalpy at saturation at 155 bars
			_fluides[0] = new StiffenedGas(rho_sat_v,1.55e7,_Tsat,esatv, sound_speed_v,cv_v); //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
			_fluides[1] = new StiffenedGas(rho_sat_l,1.55e7,_Tsat,esatl, sound_speed_l,cv_l); //stiffened gas law for water at pressure 155 bar
		}
	}
	_latentHeat=_hsatv-_hsatl;
	cout<<"Liquid saturation enthalpy "<< _hsatl<<" J/Kg"<<endl;
	*_runLogFile<<"Liquid saturation enthalpy "<< _hsatl<<" J/Kg"<<endl;
	cout<<"Vapour saturation enthalpy "<< _hsatv<<" J/Kg"<<endl;
	*_runLogFile<<"Vapour saturation enthalpy "<< _hsatv<<" J/Kg"<<endl;
	cout<<"Latent heat "<< _latentHeat<<endl;
	*_runLogFile<<"Latent heat "<< _latentHeat<<endl;
}

void DriftModel::initialize(){
	cout<<"Initialising the drift model"<<endl;
	*_runLogFile<<"Initialising the drift model"<<endl;

	_Uroe = new double[_nVar];
	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+2]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];

	if(_saveVelocity)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 (3 composantes) pour le posttraitement des lignes de courant

	if(_saveAllFields)
	{
		_VoidFraction=Field("Void fraction",CELLS,_mesh,1);
		_Enthalpy=Field("Enthalpy",CELLS,_mesh,1);
		_Concentration=Field("Concentration",CELLS,_mesh,1);
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		_mixtureDensity=Field("Mixt density",CELLS,_mesh,1);
		_Temperature=Field("Temperature",CELLS,_mesh,1);
		_DensiteLiquide=Field("Liquid density",CELLS,_mesh,1);
		_DensiteVapeur=Field("Steam density",CELLS,_mesh,1);
		_EnthalpieLiquide=Field("Liquid enthalpy",CELLS,_mesh,1);
		_EnthalpieVapeur=Field("Steam enthalpy",CELLS,_mesh,1);
		_VitesseX=Field("Velocity x",CELLS,_mesh,1);
		if(_Ndim>1)
		{
			_VitesseY=Field("Velocity y",CELLS,_mesh,1);
			if(_Ndim>2)
				_VitesseZ=Field("Velocity z",CELLS,_mesh,1);
		}
	}

	if(_entropicCorrection)
		_entropicShift=vector<double>(3);//at most 3 distinct eigenvalues

	ProblemFluid::initialize();
}

bool DriftModel::iterateTimeStep(bool &converged)
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
				//*_runLogFileogFile<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
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
void DriftModel::computeNewtonVariation()
{
	if(_timeScheme==Explicit || !_usePrimitiveVarsInNewton)
		ProblemFluid::computeNewtonVariation();
	else
	{
		if(_verbose)
		{
			cout<<"Vecteur courant Vk "<<endl;
			VecView(_primitiveVars,PETSC_VIEWER_STDOUT_SELF);
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
			if(_verbose)
			{
				cout << "Matrice du système linéaire avant contribution delta t" << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
				cout << "Second membre du système linéaire avant contribution delta t" << endl;
				VecView(_b, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
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
				cout << "Matrice du système linéaire après contribution delta t" << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
				cout << "Second membre du système linéaire après contribution delta t" << endl;
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
			if(_verbose)
			{
				cout << "solution du systeme lineaire local:" << endl;
				VecView(_newtonVariation, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
		}
	}
}

void DriftModel::convectionState( const long &i, const long &j, const bool &IsBord){
	//	sortie: etat de Roe rho, cm, v, H,T

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
		cout<<"DriftModel::convectionState Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Ui[k]<<endl;
		cout<<"DriftModel::convectionState Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Uj[k]<<endl;
	}
	if(_Ui[0]<0||_Uj[0]<0)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!densite de melange negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!densite de melange negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("densite negative, arret de calcul");
	}
	PetscScalar ri, rj, xi, xj, pi, pj;
	PetscInt Ii;
	ri = sqrt(_Ui[0]);//racine carre de phi_i rho_i
	rj = sqrt(_Uj[0]);

	_Uroe[0] = ri*rj/sqrt(_porosityi*_porosityj);	//moyenne geometrique des densites
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Densité moyenne Roe gauche " << i << ": " << ri*ri << ", droite " << j << ": " << rj*rj << "->" << _Uroe[0] << endl;
	xi = _Ui[1]/_Ui[0];//cm gauche
	xj = _Uj[1]/_Uj[0];//cm droite

	_Uroe[1] = (xi*ri + xj*rj)/(ri + rj);//moyenne de Roe des concentrations
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Concentration de Roe  gauche " << i << ": " << xi << ", droite " << j << ": " << xj << "->" << _Uroe[1] << endl;
	for(int k=0;k<_Ndim;k++){
		xi = _Ui[k+2];//phi rho u gauche
		xj = _Uj[k+2];//phi rho u droite
		_Uroe[2+k] = (xi/ri + xj/rj)/(ri + rj);
		//"moyenne" des vitesses
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout << "Vitesse de Roe composante "<< k<<"  gauche " << i << ": " << xi/(ri*ri) << ", droite " << j << ": " << xj/(rj*rj) << "->" << _Uroe[k+2] << endl;
	}
	// H = (rho E + p)/rho
	xi = _Ui[_nVar-1];//phi rho E
	xj = _Uj[_nVar-1];
	Ii = i*_nVar+1; // correct Kieu
	VecGetValues(_primitiveVars, 1, &Ii, &pi);// _primitiveVars pour p
	if(IsBord)
	{
		consToPrim(_Uj,_Vj,_porosityj);
		pj =  _Vj[1];
	}
	else
	{
		Ii = j*_nVar+1; // correct Kieu
		VecGetValues(_primitiveVars, 1, &Ii, &pj);
	}
	xi = (xi + pi)/_Ui[0];
	xj = (xj + pj)/_Uj[0];
	_Uroe[_nVar-1] = (ri*xi + rj*xj)/(ri + rj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Enthalpie totale de Roe H  gauche " << i << ": " << xi << ", droite " << j << ": " << xj << "->" << _Uroe[_nVar-1] << endl;
	// Moyenne de Roe de Tg et Td
	Ii = i*_nVar+_nVar-1;
	VecGetValues(_primitiveVars, 1, &Ii, &xi);// _primitiveVars pour T
	if(IsBord)
	{
		//consToPrim(_Uj,_Vj,_porosityj);//Fonction déjà appelée
		xj =  _Vj[_nVar-1];
	}
	else
	{
		Ii = j*_nVar+_nVar-1;
		VecGetValues(_primitiveVars, 1, &Ii, &xj);
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
}

void DriftModel::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat de diffusion (rho, rho cm, q, T)

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
		cout << "DriftModel::diffusionStateAndMatrices cellule gauche" << i << endl;
		cout << "Ui = ";
		for(int q=0; q<_nVar; q++)
			cout << _Ui[q]  << "\t";
		cout << endl;
		cout << "DriftModel::diffusionStateAndMatrices cellule droite" << j << endl;
		cout << "Uj = ";
		for(int q=0; q<_nVar; q++)
			cout << _Uj[q]  << "\t";
		cout << endl;
	}

	for(int k=0; k<_nVar; k++)
		_Udiff[k] = (_Ui[k]/_porosityi+_Uj[k]/_porosityj)/2;
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "DriftModel::diffusionStateAndMatrices conservative diffusion state" << endl;
		cout << "_Udiff = ";
		for(int q=0; q<_nVar; q++)
			cout << _Udiff[q]  << "\t";
		cout << endl;
		cout << "porosite gauche= "<<_porosityi<< ", porosite droite= "<<_porosityj<<endl;
	}

	consToPrim(_Udiff,_phi,1);
	_Udiff[_nVar-1]=_phi[_nVar-1];

	double Tm=_phi[_nVar-1];
	double RhomEm=_Udiff[_nVar-1];
	_Udiff[_nVar-1]=Tm;

	if(_timeScheme==Implicit)
	{
		double q_2=0;
		for (int l = 0; l<_Ndim;l++)
			q_2+=_Udiff[l+2]*_Udiff[l+2];
		double pression=_phi[1];
		double m_v=_Udiff[1];
		double m_l=_Udiff[0]-_Udiff[1];
		double rho_v=_fluides[0]->getDensity(pression,Tm);
		double rho_l=_fluides[1]->getDensity(pression,Tm);
		double alpha_v=m_v/rho_v,alpha_l=1-alpha_v;

		for(int l=0; l<_nVar*_nVar;l++)
			_Diffusion[l] = 0;
		double mu = alpha_v*_fluides[0]->getViscosity(Tm)+alpha_l*_fluides[1]->getViscosity(Tm);
		for(int l=2;l<(_nVar-1);l++)
		{
			_Diffusion[l*_nVar] =  mu*_Udiff[l]/(_Udiff[0]*_Udiff[0]);
			_Diffusion[l*_nVar+l] = -mu/_Udiff[0];
		}
		double lambda = alpha_v*_fluides[0]->getConductivity(Tm)+alpha_l*_fluides[1]->getConductivity(Tm);
		double C_v=  alpha_v*_fluides[0]->constante("cv");
		double C_l=	 alpha_l*_fluides[1]->constante("cv");
		double ev0=_fluides[0]->getInternalEnergy(0,rho_v);//Corriger
		double el0=_fluides[1]->getInternalEnergy(0,rho_l);//Corriger
		double Rhomem=RhomEm-0.5*q_2/(_Udiff[0]*_Udiff[0]);
		int q = (_nVar-1)*_nVar;
		//Formules a verifier
		_Diffusion[q]=lambda*((Rhomem-m_v*ev0-m_l*el0)*C_l/((m_v*C_v+m_l*C_l)*(m_v*C_v+m_l*C_l))+el0/(m_v*C_v+m_l*C_l)-q_2/(2*_Udiff[0]*_Udiff[0]*(m_v*C_v+m_l*C_l)));
		_Diffusion[q+1]=lambda*((Rhomem-m_v*ev0-m_l*el0)*(C_v-C_l)/((m_v*C_v+m_l*C_l)*(m_v*C_v+m_l*C_l))+(ev0+el0)/(m_v*C_v+m_l*C_l));
		for(int k=2;k<(_nVar-1);k++)
		{
			_Diffusion[q+k]= lambda*_Udiff[k]/(_Udiff[0]*(m_v*C_v+m_l*C_l));
		}
		_Diffusion[_nVar*_nVar-1]=-lambda/(m_v*C_v+m_l*C_l);
		/*Affichages */
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Matrice de diffusion D, pour le couple (" << i << "," << j<< "):" << endl;
			displayMatrix( _Diffusion,_nVar," Matrice de diffusion ");
		}
	}
}

void DriftModel::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	int k;
	double v2=0, q_n=0;//q_n=quantité de mouvement normale à la face interne
	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+2)]*normale[k];

	double porosityj=_porosityField(j);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup<< ", inner cell j= "<<j << " face unit normal vector "<<endl;
		for(k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall || _limitField[nameOfGroup].bcType==InnerWall){
		//Pour la convection, inversion du sens de la vitesse normale
		for(k=0; k<_Ndim; k++)
			_externalStates[(k+2)]-= 2*q_n*normale[k];

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		//Pour la diffusion, paroi à vitesse et temperature imposees
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double concentration=_Vj[0];
		double pression=_Vj[1];
		double T=_limitField[nameOfGroup].T;
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
		{
			cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
			cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::setBoundaryState: Inlet, impossible to compute mixture density, division by zero");
		}

		_externalStates[0]=porosityj*rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
		_externalStates[1]=concentration*_externalStates[0];

		_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
		v2 +=_limitField[nameOfGroup].v_x[0]*_limitField[nameOfGroup].v_x[0];
		if(_Ndim>1)
		{
			v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
			if(_Ndim==3)
			{
				_externalStates[4]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
				v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
			}
		}
		_externalStates[_nVar-1] = _externalStates[1]*_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho_v)
																																																																																																																																					 +(_externalStates[0]-_externalStates[1])*_fluides[1]->getInternalEnergy(_limitField[nameOfGroup].T,rho_l) + _externalStates[0]*v2/2;
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){

		if(q_n<=0){
			VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
			double concentration=_limitField[nameOfGroup].conc;
			double pression=_Vj[1];
			double T=_limitField[nameOfGroup].T;
			double rho_v=_fluides[0]->getDensity(pression,T);
			double rho_l=_fluides[1]->getDensity(pression,T);
			if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
			{
				cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
				cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
				*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
				_runLogFile->close();
				throw CdmathException("DriftModel::setBoundaryState: Inlet, impossible to compute mixture density, division by zero");
			}

			_externalStates[0]=porosityj*rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
			_externalStates[1]=concentration*_externalStates[0];
			_externalStates[2]=_externalStates[0]*(_limitField[nameOfGroup].v_x[0]);
			v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
			if(_Ndim>1)
			{
				v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
				_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
				if(_Ndim==3)
				{
					_externalStates[4]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
					v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[1]*_fluides[0]->getInternalEnergy(T,rho_v)+(_externalStates[0]-_externalStates[1])*_fluides[1]->getInternalEnergy(T,rho_l) + _externalStates[0]*v2/2;
		}
		else if(_nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho

		//Building the external state
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double concentration, Tm;
		if(q_n<=0){
			concentration=_limitField[nameOfGroup].conc;
			Tm=_limitField[nameOfGroup].T;
		}
		else{
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for concentration, velocity and temperature"<<endl;
			concentration=_Vj[0];
			Tm=_Vj[_nVar-1];
		}

		double pression=_limitField[nameOfGroup].p + hydroPress;
		double rho_v=_fluides[0]->getDensity(pression, Tm);
		double rho_l=_fluides[1]->getDensity(pression, Tm);
		if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
		{
			cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
			cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::jacobian: Inlet, impossible to compute mixture density, division by zero");
		}

		_externalStates[0]=porosityj*rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
		_externalStates[1]=concentration*_externalStates[0];
		double mv=_externalStates[1], ml=_externalStates[0]-_externalStates[1];
		for(k=0; k<_Ndim; k++)
		{
			v2+=_Vj[k+2]*_Vj[k+2];
			_externalStates[k+2]=_externalStates[0]*_Vj[(k+2)] ;
		}
		_externalStates[_nVar-1] = mv*_fluides[0]->getInternalEnergy(Tm,rho_v)+ml*_fluides[1]->getInternalEnergy(Tm,rho_l) +_externalStates[0]* v2/2;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		if(q_n<=0  &&  _nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition for concentration, velocity and temperature"<<endl;

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=(Cj.x()-_gravityReferencePoint[0])*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=(Cj.y()-_gravityReferencePoint[1])*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=(Cj.z()-_gravityReferencePoint[2])*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho

		//Building the external state
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

		double concentration=_Vj[0];
		double pression=_limitField[nameOfGroup].p+hydroPress;
		double Tm=_Vj[_nVar-1];
		double rho_v=_fluides[0]->getDensity(pression, Tm);
		double rho_l=_fluides[1]->getDensity(pression, Tm);
		if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
		{
			cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
			cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::jacobian: Inlet, impossible to compute mixture density, division by zero");
		}

		_externalStates[0]=porosityj*rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
		_externalStates[1]=concentration*_externalStates[0];
		double mv=_externalStates[1], ml=_externalStates[0]-_externalStates[1];
		for(k=0; k<_Ndim; k++)
		{
			v2+=_Vj[k+2]*_Vj[k+2];
			_externalStates[k+2]=_externalStates[0]*_Vj[(k+2)] ;
		}
		_externalStates[_nVar-1] = mv*_fluides[0]->getInternalEnergy(Tm,rho_v)+ml*_fluides[1]->getInternalEnergy(Tm,rho_l) +_externalStates[0]* v2/2;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else {
		cout<<"!!!!!!!!!!!!!!! Error DriftModel::setBoundaryState !!!!!!!!!!"<<endl;
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Boundary condition not set for boundary named "<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, InnerWall, Inlet, and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named "<<nameOfGroup<<"Accepted boundary condition are Neumann,Wall, InnerWall, Inlet, and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}
}

void DriftModel::convectionMatrices()
{
	//entree: URoe = rho, cm, u, H
	//sortie: matrices Roe+  et Roe- +Roe si scheme centre

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"DriftModel::convectionMatrices()"<<endl;

	double umn=0, u_2=0; //valeur de u.normale et |u|²
	for(int i=0;i<_Ndim;i++)
	{
		u_2 += _Uroe[2+i]*_Uroe[2+i];
		umn += _Uroe[2+i]*_vec_normal[i];//vitesse normale
	}

	vector<complex<double> > vp_dist(3);

	if(_spaceScheme==staggered && _nonLinearFormulation==VFFC)//special case where we need no Roe matrix
	{
		staggeredVFFCMatricesConservativeVariables(umn);//Computation of classical upwinding matrices
		if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)//For use in implicit matrix
			staggeredVFFCMatricesPrimitiveVariables(umn);
	}
	else//In this case we build a Roe matrix
	{
		double rhom=_Uroe[0];
		double cm=_Uroe[1];
		double Hm=_Uroe[_nVar-1];
		double hm=Hm-0.5*u_2;
		double m_v=cm*rhom, m_l=rhom-m_v;
		double Tm;
		Vector vitesse(_Ndim);
		for(int idim=0;idim<_Ndim;idim++)
			vitesse[idim]=_Uroe[2+idim];

		if(cm<_precision)//pure liquid
			Tm=_fluides[1]->getTemperatureFromEnthalpy(hm,rhom);
		else if(cm>1-_precision)
			Tm=_fluides[0]->getTemperatureFromEnthalpy(hm,rhom);
		else//Hypothèse de saturation
			Tm=_Tsat;

		double pression= getMixturePressure(cm, rhom, Tm);

		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<"Roe state rhom="<<rhom<<" cm= "<<cm<<" Hm= "<<Hm<<" Tm= "<<Tm<<" pression "<<pression<<endl;

		getMixturePressureDerivatives( m_v, m_l, pression, Tm);//EOS is involved to express pressure jump and sound speed
		if(_kappa*hm+_khi+cm*_ksi<0){
			*_runLogFile<<"Calcul matrice de Roe: vitesse du son complexe"<<endl;
			_runLogFile->close();
			throw CdmathException("Calcul matrice de Roe: vitesse du son complexe");
		}
		double am=sqrt(_kappa*hm+_khi+cm*_ksi);//vitesse du son du melange
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<"DriftModel::convectionMatrices, sound speed am= "<<am<<endl;

		//On remplit la matrice de Roe
		double ecin=0.5*u_2;
		RoeMatrixConservativeVariables( cm, umn, ecin, Hm,vitesse);

		//On remplit les valeurs propres
		vp_dist[0]=umn+am;
		vp_dist[1]=umn-am;
		vp_dist[2]=umn;

		_maxvploc=fabs(umn)+am;
		if(_maxvploc>_maxvp)
			_maxvp=_maxvploc;

		/* Construction des matrices de decentrement */
		if(_spaceScheme== centered){
			if(_entropicCorrection)
			{
				*_runLogFile<<"DriftModel::convectionMatrices: entropy schemes not yet available for drift model"<<endl;
				_runLogFile->close();
				throw CdmathException("DriftModel::convectionMatrices: entropy schemes not yet available for drift model");
			}

			for(int i=0; i<_nVar*_nVar;i++)
				_absAroe[i] = 0;
			if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)
				for(int i=0; i<_nVar*_nVar;i++)
					_absAroeImplicit[i] = 0;
		}
		else if( _spaceScheme ==staggered){
			//Calcul de décentrement de type décalé pour formulation Roe
			staggeredRoeUpwindingMatrixConservativeVariables( cm, umn, ecin, Hm, vitesse);
			//staggeredRoeUpwindingMatrixPrimitiveVariables( cm, umn, ecin, Hm, vitesse);
		}
		else if(_spaceScheme == upwind || _spaceScheme ==pressureCorrection || _spaceScheme ==lowMach ){
			if(_entropicCorrection)
				entropicShift(_vec_normal);
			else
				_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

			vector< complex< double > > y (3,0);
			Polynoms Poly;
			for( int i=0 ; i<3 ; i++)
				y[i] = Poly.abs_generalise(vp_dist[i])+1*_entropicShift[i];
			Poly.abs_par_interp_directe(3,vp_dist, _Aroe, _nVar,_precision, _absAroe,y);

			if( _spaceScheme ==pressureCorrection)
			{
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(2+i)*_nVar+2+j]-=(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			}
			else if( _spaceScheme ==lowMach){
				double M=sqrt(u_2)/am;
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(2+i)*_nVar+2+j]-=(1-M)*(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			}
		}
		else
		{
			*_runLogFile<<"DriftModel::convectionMatrices: scheme not treated"<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::convectionMatrices: scheme not treated");
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
				_Vij[0]=cm;
				_Vij[1]=pression;//pressure
				_Vij[_nVar-1]=Tm;//Temperature
				for(int idim=0;idim<_Ndim; idim++)
					_Vij[2+idim]=_Uroe[2+idim];
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
		displayMatrix(_AroePlusImplicit, _nVar,"Matrice _AroePlusImplicit");
	}

	/*******Calcul de la matrice signe pour VFFC, VFRoe et décentrement des termes source***/
	if(_entropicCorrection)
	{
		InvMatriceRoe( vp_dist);
		Polynoms Poly;
		Poly.matrixProduct(_absAroe, _nVar, _nVar, _invAroe, _nVar, _nVar, _signAroe);
	}
	else if (_spaceScheme==upwind  || (_spaceScheme ==lowMach) || (_spaceScheme ==pressureCorrection))//upwind sans entropic
		SigneMatriceRoe( vp_dist);
	else if(_spaceScheme== centered)//centre  sans entropic
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
	else if(_spaceScheme ==staggered )
	{
		double signu;
		if(umn>0)
			signu=1;
		else if (umn<0)
			signu=-1;
		else
			signu=0;
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
		_signAroe[0] = signu;
		_signAroe[1+_nVar] = signu;
		for(int i=2; i<_nVar-1;i++)
			_signAroe[i*_nVar+i] = -signu;
		_signAroe[_nVar*(_nVar-1)+_nVar-1] = signu;
	}
	else
	{
		*_runLogFile<<"DriftModel::convectionMatrices: well balanced option not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::convectionMatrices: well balanced option not treated");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0 && _timeScheme==Implicit)
		displayMatrix(_signAroe, _nVar,"Signe de la matrice de Roe");
}

void DriftModel::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)
{
	double Tm=_Udiff[_nVar-1];
	double lambdal=_fluides[1]->getConductivity(Tm);
	double lambdav = _fluides[0]->getConductivity(Tm);
	double mu_l = _fluides[1]->getViscosity(Tm);
	double mu_v = _fluides[0]->getViscosity(Tm);

	if(mu_v==0 && mu_l ==0 && lambdav==0 && lambdal==0 && _heatTransfertCoeff==0)
		return;

	//extraction des valeurs
	for(int k=0; k<_nVar; k++)
		_idm[k] = _nVar*i + k;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion: variables primitives maille " << i<<endl;
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
		lambdal=max(lambdal,_heatTransfertCoeff);//wall nucleate boing -> larger heat transfer
		for(int k=0; k<_nVar; k++)
			_idn[k] = k;

		VecGetValues(_Uextdiff, _nVar, _idn, _Uj);
		consToPrim(_Uj,_Vj,1);
	}
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion: variables primitives maille " <<j <<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q] << endl;
		cout << endl;
	}
	double pression=(_Vi[1]+_Vj[1])/2;//ameliorer car traitement different pour pression et temperature
	double m_v=_Udiff[1];
	double rho_v=_fluides[0]->getDensity(pression,Tm);
	double rho_l=_fluides[1]->getDensity(pression,Tm);
	double alpha_v=m_v/rho_v,alpha_l=1-alpha_v;
	double mu = alpha_v*mu_v+alpha_l*mu_l;
	double lambda = alpha_v*lambdav+alpha_l*lambdal;

	//pas de diffusion sur la masse totale
	_phi[0]=0;
	//on n'a pas encore mis la contribution sur la masse
	_phi[1]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=2; k<_nVar-1; k++)
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
}


void DriftModel::sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i)
{
	double phirho=Ui[0],phim1=Ui[1],phim2=phirho-phim1,phirhoE=Ui[_nVar-1], cv=Vi[0],  P=Vi[1], T=Vi[_nVar-1];
	double norm_u=0;
	for(int k=0; k<_Ndim; k++)
		norm_u+=Vi[2+k]*Vi[2+k];
	norm_u=sqrt(norm_u);
	double h=(phirhoE-0.5*phirho*norm_u*norm_u+P*_porosityField(i))/phirho;//e+p/rho

	Si[0]=0;
	//if(T>_Tsat && cv<1-_precision)
	if(_hsatv>h  && h>_hsatl && cv<1-_precision)//heated boiling//
		Si[1]=_heatPowerField(i)/_latentHeat*_porosityField(i);//phi*gamma
	else if (P<_Psat && cv<1-_precision)//flash boiling
		Si[1]=-_dHsatl_over_dp*_dp_over_dt(i)/_latentHeat;
	else
		Si[1]=0;
	for(int k=2; k<_nVar-1; k++)
		Si[k]  =_gravite[k]*phirho-(phim1*_dragCoeffs[0]+phim2*_dragCoeffs[1])*norm_u*Vi[k];

	Si[_nVar-1]=_heatPowerField(i);

	for(int k=0; k<_Ndim; k++)
		Si[_nVar-1] +=(_GravityField3d[k]*phirho-(phim1*_dragCoeffs[0]+phim2*_dragCoeffs[1])*norm_u*Vi[2+k])*Vi[2+k];

	if(_timeScheme==Implicit)
	{
		for(int k=0; k<_nVar*_nVar;k++)
			_GravityImplicitationMatrix[k] = 0;
		if(!_usePrimitiveVarsInNewton)
			for(int k=0; k<_nVar;k++)
				_GravityImplicitationMatrix[k*_nVar]=-_gravite[k];
		else
		{
			getDensityDerivatives( cv, P, T ,norm_u*norm_u);
			for(int k=0; k<_nVar;k++)
			{
				_GravityImplicitationMatrix[k*_nVar+0]      =-_gravite[k]*_drho_sur_dcv;
				_GravityImplicitationMatrix[k*_nVar+1]      =-_gravite[k]*_drho_sur_dp;
				_GravityImplicitationMatrix[k*_nVar+_nVar-1]=-_gravite[k]*_drho_sur_dT;
			}
		}
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"DriftModel::sourceVector"<<endl;
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

void DriftModel::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{
	double norm_u=0, u_n=0, phirho;
	for(int i=0;i<_Ndim;i++)
		u_n += _Uroe[2+i]*_vec_normal[i];

	pressureLoss[0]=0;
	pressureLoss[1]=0;
	if(u_n>0){
		for(int i=0;i<_Ndim;i++)
			norm_u += Vi[2+i]*Vi[2+i];
		norm_u=sqrt(norm_u);
		phirho=Ui[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[2+i]=-K*phirho*norm_u*Vi[2+i];
	}
	else{
		for(int i=0;i<_Ndim;i++)
			norm_u += Vj[2+i]*Vj[2+i];
		norm_u=sqrt(norm_u);
		phirho=Uj[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[2+i]=-K*phirho*norm_u*Vj[2+i];
	}
	pressureLoss[_nVar-1]=-K*phirho*norm_u*norm_u*norm_u;


	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"DriftModel::pressureLossVector K= "<<K<<endl;
		cout<<"pressure loss vector phirho= "<< phirho<<" norm_u= "<<norm_u<<endl;
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
void DriftModel::porosityGradientSourceVector()
{
	double u_ni=0, u_nj=0, rhoi,rhoj, pi=_Vi[1], pj=_Vj[1],pij;
	for(int i=0;i<_Ndim;i++) {
		u_ni += _Vi[2+i]*_vec_normal[i];
		u_nj += _Vj[2+i]*_vec_normal[i];
	}
	_porosityGradientSourceVector[0]=0;
	_porosityGradientSourceVector[1]=0;
	rhoj=_Uj[0]/_porosityj;
	rhoi=_Ui[0]/_porosityi;
	pij=(pi+pj)/2+rhoi*rhoj/2/(rhoi+rhoj)*(u_ni-u_nj)*(u_ni-u_nj);
	for(int i=0;i<_Ndim;i++)
		_porosityGradientSourceVector[2+i]=pij*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
	_porosityGradientSourceVector[_nVar-1]=0;
}

void DriftModel::jacobian(const int &j, string nameOfGroup,double * normale)
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
	if (_limitField[nameOfGroup].bcType==Wall || _limitField[nameOfGroup].bcType==InnerWall)
	{
		for(k=0; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
		for(k=2; k<_nVar-1;k++)
			for(int l=2; l<_nVar-1;l++)
				_Jcb[k*_nVar + l] -= 2*normale[k-2]*normale[l-2];
	}
	else if (_limitField[nameOfGroup].bcType==Inlet)
	{
		return;//For the moment no implicitation, should debug inlet
		if(q_n<0){
			//Boundary state quantities
			double v[_Ndim], ve[_Ndim], v2, ve2;
			VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
			double concentration=_limitField[nameOfGroup].conc;
			double pression=_Vj[1];
			double T=_limitField[nameOfGroup].T;
			double rho_v=_fluides[0]->getDensity(pression,T);
			double rho_l=_fluides[1]->getDensity(pression,T);
			if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
			{
				cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
				cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
				*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
				_runLogFile->close();
				throw CdmathException("DriftModel::jacobian: Inlet, impossible to compute mixture density, division by zero");
			}

			double rho_e=rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
			double e_v=_fluides[0]->getInternalEnergy(T,rho_v);
			double e_l=_fluides[1]->getInternalEnergy(T,rho_l);
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
			double E_e = concentration*e_v+(1-concentration)*e_l + ve2/2;

			//Pressure differential
			double gamma_v =_fluides[0]->constante("gamma");
			double gamma_l =_fluides[1]->constante("gamma");
			double omega=concentration/((gamma_v-1)*rho_v*rho_v*e_v)+(1-concentration)/((gamma_l-1)*rho_l*rho_l*e_l);
			double rhom=_externalStates[0];
			double m_v=concentration*rhom, m_l=rhom-m_v;
			getMixturePressureDerivatives( m_v, m_l, pression, T);

			double CoeffCol1 = rho_e*rho_e*omega*(_khi+_kappa*v2/2);
			double CoeffCol2 = rho_e*rho_e*omega*_ksi;
			double CoeffCol3et4 = rho_e*rho_e*omega*_kappa;

			//Mass line
			_Jcb[0]=CoeffCol1;
			_Jcb[1]=CoeffCol2;
			for(k=0;k<_Ndim;k++)
				_Jcb[2+k]=-CoeffCol3et4*v[k];
			_Jcb[_nVar-1]=CoeffCol3et4;
			//vapour mass line
			for(k=0;k<_nVar;k++)
				_Jcb[_nVar+k]=_Jcb[k]*concentration;
			//Momentum lines
			for(int l=0;l<_Ndim;l++)
				for(k=0;k<_nVar;k++)
					_Jcb[(2+l)*_nVar+k]=_Jcb[k]*ve[l];
			//Energy line
			for(k=0;k<_nVar;k++)
				_Jcb[(_nVar-1)*_nVar+k]=_Jcb[k]*E_e;
		}
		else
			for(k=0;k<_nVar;k++)
				_Jcb[k*_nVar+k]=1;
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure && q_n<0)
	{
		return;//For the moment no implicitation, should debug inletPressure
		//Boundary state quantities
		double v[_Ndim], v2=0;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double concentration=_limitField[nameOfGroup].conc;
		double pression=_limitField[nameOfGroup].p;
		double T=_limitField[nameOfGroup].T;
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
		{
			cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
			cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::jacobian: InletPressure, impossible to compute mixture density, division by zero");
		}

		double rho_ext=rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
		double rho_int=_externalStates[0];
		double ratio_densite=rho_ext/rho_int;
		for(k=0;k<_Ndim;k++){
			v[k]=_Vj[2+k];
			v2+=v[k]*v[k];
		}
		//Momentum lines
		for(int l=0;l<_Ndim;l++){
			_Jcb[(2+l)*_nVar]=-ratio_densite*v[l];
			_Jcb[(2+l)*_nVar+2+l]=ratio_densite;
		}
		//Energy line
		_Jcb[(_nVar-1)*_nVar]=-ratio_densite*v2;
		for(int l=0;l<_Ndim;l++)
			_Jcb[(_nVar-1)*_nVar+2+l]=ratio_densite*v[l];

	}
	else if (_limitField[nameOfGroup].bcType==Outlet || (_limitField[nameOfGroup].bcType==InletPressure && q_n>=0)){
		return;//For the moment no implicitation, should debug inletPressure
		//Boundary state quantities
		double v[_Ndim], v2;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double concentration=_Vj[0];
		double pression=_limitField[nameOfGroup].p;
		double T=_Vj[_nVar-1];
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
		{
			cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
			cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::jacobian: Outlet, impossible to compute mixture density, division by zero");
		}

		double rho_ext=rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v);
		double rho_int=_externalStates[0];
		double e_v=_fluides[0]->getInternalEnergy(T,rho_v);
		double e_l=_fluides[1]->getInternalEnergy(T,rho_l);
		double ratio_densite=rho_ext/rho_int;
		for(k=0;k<_Ndim;k++){
			v[k]=_Vj[2+k];
			v2+=v[k]*v[k];
		}
		double E_m = concentration*e_v+(1-concentration)*e_l + v2/2;//total energy

		//Temperature differential
		double C_v=  _fluides[0]->constante("cv");
		double C_l=	 _fluides[1]->constante("cv");
		double ev0=_fluides[0]->getInternalEnergy(0,rho_v);//Corriger
		double el0=_fluides[1]->getInternalEnergy(0,rho_l);//Corriger

		double omega=concentration*C_v/(rho_v*e_v)+(1-concentration)*C_l/(rho_l*e_l);
		double rhomem=_externalStates[0]*(concentration*e_v+(1-concentration)*e_l);
		double m_v=concentration*rho_ext, m_l=rho_ext-m_v;
		double alpha_v=m_v/rho_v,alpha_l=1-alpha_v;
		double denom=m_v *C_v+m_l* C_l;

		double khi=(m_v*(ev0*C_l-el0*C_v)-rhomem*C_l)/(denom*denom);
		double ksi=(rho_ext*(el0*C_v-ev0*C_l)+rhomem*(C_l-C_v))/(denom*denom);
		double kappa=1/denom;

		double CoeffCol1 = rho_int*rho_int*omega*(khi+kappa*v2/2)+ratio_densite;//The +ratio_densite helps filling the lines other than total mass
		double CoeffCol2 = rho_int*rho_int*omega*ksi;
		double CoeffCol3et4 = rho_int*rho_int*omega*kappa;

		//Mass line
		_Jcb[0]=-CoeffCol1;
		_Jcb[1]=-CoeffCol2;
		for(k=0;k<_Ndim;k++)
			_Jcb[2+k]=CoeffCol3et4*v[k];
		_Jcb[_nVar-1]=-CoeffCol3et4;
		//vapour mass line
		for(k=0;k<_nVar;k++)
			_Jcb[_nVar+k]=_Jcb[k]*concentration;
		//Momentum lines
		for(int l=0;l<_Ndim;l++)
			for(k=0;k<_nVar;k++)
				_Jcb[(2+l)*_nVar+k]=_Jcb[k]*v[l];
		//Energy line
		for(k=0;k<_nVar;k++)
			_Jcb[(_nVar-1)*_nVar+k]=_Jcb[k]*E_m;

		//adding the remaining diagonal term
		for(k=0;k<_nVar;k++)
			_Jcb[k*_nVar+k]+=ratio_densite;

	}
	else  if (_limitField[nameOfGroup].bcType!=Neumann)
	{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::jacobian: The boundary condition is not recognised: neither inlet, inltPressure, outlet, wall, InnerWall, nor neumann");
	}
}

//calcule la jacobienne pour une CL de diffusion
void  DriftModel::jacobianDiff(const int &j, string nameOfGroup)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite diffusion bord "<< nameOfGroup<<endl;


	double v2=0,cd = 1,cn=0,p0, gamma;
	int idim,jdim,k;

	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;

	if (_limitField[nameOfGroup].bcType==Wall || _limitField[nameOfGroup].bcType==InnerWall){
	}
	else if (_limitField[nameOfGroup].bcType==Outlet || _limitField[nameOfGroup].bcType==Neumann
			||_limitField[nameOfGroup].bcType==Inlet || _limitField[nameOfGroup].bcType==InletPressure)
	{
		for(k=0;k<_nVar;k++)
			_JcbDiff[k*_nVar+k]=1;
	}
	else{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::jacobianDiff: This boundary condition is not recognised");
	}
}


void DriftModel::computeScaling(double maxvp)
{
	//	if(_spaceScheme!=staggered)
	{
		_blockDiag[0]=1;
		_invBlockDiag[0]=1;
		_blockDiag[1]=1;
		_invBlockDiag[1]=1;
		for(int q=2; q<_nVar-1; q++)
		{
			_blockDiag[q]=1./maxvp;//
			_invBlockDiag[q]= maxvp;//1.;//
		}
		_blockDiag[_nVar - 1]=1/(maxvp*maxvp);//1
		_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
	}
	/*
	else{//_spaceScheme==staggered
		_blockDiag[0]=maxvp*maxvp;
		_invBlockDiag[0]=1/_blockDiag[0];
		_blockDiag[1]=maxvp*maxvp;
		_invBlockDiag[1]=1/_blockDiag[1];
		for(int q=2; q<_nVar-1; q++)
		{
			_blockDiag[q]=1;//
			_invBlockDiag[q]= 1;//1.;//
		}
		_blockDiag[_nVar - 1]=1;//1
		_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
	}
	 */
}
Vector DriftModel::computeExtendedPrimState(double *V)
{
	Vector Vext(7+2*_Ndim);

	double C=V[0], P=V[1], T=V[_nVar-1];
	double rho_v=_fluides[0]->getDensity(P,T);
	double rho_l=_fluides[1]->getDensity(P,T);
	double e_v=_fluides[0]->getInternalEnergy(T,rho_v);
	double e_l=_fluides[1]->getInternalEnergy(T,rho_l);
	if(fabs(rho_l*C+rho_v*(1-C))<_precision)
	{
		cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<", concentration= "<<C<<endl;
		*_runLogFile<<"DriftModel::computeExtendedPrimState: impossible to compute void fraction, division by zero"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::computeExtendedPrimState: impossible to compute void fraction, division by zero");
	}

	_externalStates[0]=rho_v*rho_l/(C*rho_l+(1-C)*rho_v);
	double alpha_v=rho_l*C/(rho_l*C+rho_v*(1-C)), alpha_l=1-alpha_v;
	double rho_m=alpha_v*rho_v+alpha_l*rho_l;
	double h_m=(alpha_v*rho_v*e_v+alpha_l*rho_l*e_l+P)/rho_m;
	Vector Vit(_Ndim);
	for(int i=0;i<_Ndim;i++)
		Vit(i)=V[2+i];
	Vector u_r=relative_velocity(C, Vit,rho_m);

	Vext(0)=C;
	Vext(1)=P;
	for(int i=0;i<_Ndim;i++)
		Vext(2+i)=Vit(i);
	Vext(2+_Ndim)=T;
	Vext(3+_Ndim)=alpha_v;
	for(int i=0;i<_Ndim;i++)
		Vext(4+_Ndim+i)=u_r(i);
	Vext((4+2*_Ndim))=rho_v;
	Vext((5+2*_Ndim))=rho_l;
	Vext(6+2*_Ndim)=h_m;

	return Vext;
}


void DriftModel::primToCons(const double *P, const int &i, double *W, const int &j){
	double concentration=P[i*_nVar];
	double pression=P[i*_nVar+1];
	double temperature=P[i*_nVar+_nVar-1];
	double rho_v=_fluides[0]->getDensity(pression,temperature);
	double rho_l=_fluides[1]->getDensity(pression,temperature);
	double e_v = _fluides[0]->getInternalEnergy(temperature,rho_v);
	double e_l = _fluides[1]->getInternalEnergy(temperature,rho_l);
	if(fabs(concentration*rho_l+(1-concentration)*rho_v)<_precision)
	{
		cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<endl;
		cout<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
		*_runLogFile<<"concentration*rho_l+(1-concentration)*rho_v= "<<concentration*rho_l+(1-concentration)*rho_v<<endl;
		*_runLogFile<<"DriftModel::primToCons: impossible to compute mixture density, division by zero"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::primToCons: impossible to compute mixture density, division by zero");
	}
	W[j*(_nVar)]  =(rho_v*rho_l/(concentration*rho_l+(1-concentration)*rho_v))*_porosityField(j);//phi*rho
	W[j*(_nVar)+1]  =W[j*(_nVar)] *concentration;//phi *rho*c_v
	for(int k=0; k<_Ndim; k++)
		W[j*_nVar+(k+2)] = W[j*(_nVar)] *P[i*_nVar+(k+2)];//phi*rho*u

	W[j*_nVar+_nVar-1] = W[j*(_nVar)+1]* e_v + (W[j*(_nVar)]-W[j*(_nVar)+1])* e_l;//phi rhom e_m
	for(int k=0; k<_Ndim; k++)
		W[j*_nVar+_nVar-1] += W[j*_nVar]*P[i*_nVar+(k+2)]*P[i*_nVar+(k+2)]*0.5;//phi rhom e_m + phi rho u*u

	if(_verbose && _nbTimeStep%_freqSave ==0){
		cout<<"DriftModel::primToCons: T="<<temperature<<", pression= "<<pression<<endl;
		cout<<"rhov= "<<rho_v<<", rhol= "<<rho_l<<", rhom= "<<W[j*(_nVar)]<<" e_v="<<e_v<<" e_l="<<e_l<<endl;
	}
}
void DriftModel::primToConsJacobianMatrix(double *V)
{
	double concentration=V[0];
	double pression=V[1];
	double temperature=V[_nVar-1];
	double vitesse[_Ndim];
	for(int idim=0;idim<_Ndim;idim++)
		vitesse[idim]=V[2+idim];
	double v2=0;
	for(int idim=0;idim<_Ndim;idim++)
		v2+=vitesse[idim]*vitesse[idim];

	double rho_v=_fluides[0]->getDensity(pression,temperature);
	double rho_l=_fluides[1]->getDensity(pression,temperature);
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");

	double rho=concentration*rho_v+(1-concentration)*rho_l;;

	for(int k=0;k<_nVar*_nVar; k++)
		_primToConsJacoMat[k]=0;

	if(		!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);
		double e_v = fluide0->getInternalEnergy(temperature);
		double e_l = fluide0->getInternalEnergy(temperature);
		double cv_v=fluide0->constante("cv");
		double cv_l=fluide1->constante("cv");
		double e=concentration*e_v+(1-concentration)*e_l;
		double E=e+0.5*v2;

		/******* Masse totale **********/
		_primToConsJacoMat[0]=-rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[1]=rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		_primToConsJacoMat[_nVar-1]=-rho*rho*(
				cv_v*   concentration /(rho_v*(e_v-q_v))
				+cv_l*(1-concentration)/(rho_l*(e_l-q_l))
		);

		/******* Masse partielle **********/
		_primToConsJacoMat[_nVar]=rho-concentration*rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[_nVar+1]=concentration*rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		_primToConsJacoMat[_nVar+_nVar-1]=-concentration*rho*rho*(
				cv_v*   concentration /(rho_v*(e_v-q_v))
				+cv_l*(1-concentration)/(rho_l*(e_l-q_l))
		);
		/******* Quantité de mouvement **********/
		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[2*_nVar+idim*_nVar]=-vitesse[idim]*rho*rho*(1/rho_v-1/rho_l);
			_primToConsJacoMat[2*_nVar+idim*_nVar+1]=vitesse[idim]*rho*rho*(
					concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
					+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
			);
			_primToConsJacoMat[2*_nVar+idim*_nVar+2+idim]=rho;
			_primToConsJacoMat[2*_nVar+idim*_nVar+_nVar-1]=-vitesse[idim]*rho*rho*(
					cv_v*   concentration /(rho_v*(e_v-q_v))
					+cv_l*(1-concentration)/(rho_l*(e_l-q_l))
			);
		}
		/******* Energie totale **********/
		_primToConsJacoMat[(_nVar-1)*_nVar]=rho*(e_v-e_l)-E*rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[(_nVar-1)*_nVar+1]=E*rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+2+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*(cv_v*concentration + cv_l*(1-concentration))
																																																																																											-rho*rho*E*( cv_v*   concentration /(rho_v*(e_v-q_v))
																																																																																													+cv_l*(1-concentration)/(rho_l*(e_l-q_l)));
	}
	else if(_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);
		double h_v=fluide0->getEnthalpy(temperature);
		double h_l=fluide1->getEnthalpy(temperature);
		double h=concentration*h_v+(1-concentration)*h_l;
		double H=h+0.5*v2;
		double cp_v=fluide0->constante("cp");
		double cp_l=fluide1->constante("cp");

		/******* Masse totale **********/
		_primToConsJacoMat[0]=-rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[1]=rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		);
		_primToConsJacoMat[_nVar-1]=-rho*rho*(
				cp_v*   concentration /(rho_v*(h_v-q_v))
				+cp_l*(1-concentration)/(rho_l*(h_l-q_l))
		);

		/******* Masse partielle **********/
		_primToConsJacoMat[_nVar]=rho-concentration*rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[_nVar+1]=concentration*rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		);
		_primToConsJacoMat[_nVar+_nVar-1]=-concentration*rho*rho*(
				cp_v*   concentration /(rho_v*(h_v-q_v))
				+cp_l*(1-concentration)/(rho_l*(h_l-q_l))
		);
		/******* Quantité de mouvement **********/
		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[2*_nVar+idim*_nVar]=-vitesse[idim]*rho*rho*(1/rho_v-1/rho_l);
			_primToConsJacoMat[2*_nVar+idim*_nVar+1]=vitesse[idim]*rho*rho*(
					gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
					+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
			);
			_primToConsJacoMat[2*_nVar+idim*_nVar+2+idim]=rho;
			_primToConsJacoMat[2*_nVar+idim*_nVar+_nVar-1]=-vitesse[idim]*rho*rho*(
					cp_v*   concentration /(rho_v*(h_v-q_v))
					+cp_l*(1-concentration)/(rho_l*(h_l-q_l))
			);
		}
		/******* Energie totale **********/
		_primToConsJacoMat[(_nVar-1)*_nVar]=rho*(h_v-h_l)-H*rho*rho*(1/rho_v-1/rho_l);
		_primToConsJacoMat[(_nVar-1)*_nVar+1]=H*rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		)-1;
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+2+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*(cp_v*concentration + cp_l*(1-concentration))
																																																																																											-rho*rho*H*(cp_v*   concentration /(rho_v*(h_v-q_v))
																																																																																													+cp_l*(1-concentration)/(rho_l*(h_l-q_l)));
	}
	else
	{
		*_runLogFile<<"DriftModel::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" DriftModel::primToConsJacobianMatrix" << endl;
		displayVector(_Vi,_nVar," _Vi " );
		cout<<"rho_v= "<<rho_v<<", rho_l= "<<rho_l<<endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobienne primToCons: ");
	}
}

void DriftModel::consToPrim(const double *Wcons, double* Wprim,double porosity)
{
	double m_v=Wcons[1]/porosity;
	double m_l=(Wcons[0]-Wcons[1])/porosity;
	_minm1=min(m_v,_minm1);
	_minm2=min(m_l,_minm2);
	if(m_v<-_precision || m_l<-_precision)
		_nbMaillesNeg+=1;
	else if(m_v< 0 && m_v>-_precision )
		m_v=0;
	else if(m_l< 0 && m_l>-_precision )
		m_l=0;
	double concentration=m_v/(m_v+m_l);
	double q_2 = 0;
	for(int k=0;k<_Ndim;k++)
		q_2 += Wcons[k+2]*Wcons[k+2];
	q_2 /= Wcons[0];	//phi rho u²

	double rho_m_e_m=(Wcons[_nVar-1] -0.5*q_2)/porosity;
	double pression,Temperature;

	if(concentration<_precision)
	{
		pression=_fluides[1]->getPressure(rho_m_e_m,m_l);
		Temperature=_fluides[1]->getTemperatureFromPressure(pression,m_l);
	}
	else if(concentration>1-_precision)
	{
		pression=_fluides[0]->getPressure(rho_m_e_m,m_v);
		Temperature=_fluides[0]->getTemperatureFromPressure(pression,m_v);
	}
	else//Hypothèses de saturation
	{
		//Première approche
		getMixturePressureAndTemperature(concentration, m_v+m_l, rho_m_e_m, pression, Temperature);
		//cout<<"Temperature= "<<Temperature<<", pression= "<<pression<<endl;
		//throw CdmathException("DriftModel::consToPrim: Apparition de vapeur");

		//Seconde approche : on impose la température directement
		//Temperature=_Tsat;
		//pression=getMixturePressure(concentration, m_v+m_l, Temperature);

		//Troisieme approche : on impose les enthalpies de saturation
		//double rho_m_h_m= m_v*_hsatv + m_l*_hsatl;
		//pression = rho_m_h_m - rho_m_e_m;
		//Temperature = getMixtureTemperature(concentration, m_v+m_l, pression);
	}

	if (pression<0){
		cout << "pressure = "<< pression << " < 0 " << endl;
		*_runLogFile << "pressure = "<< pression << " < 0 " << endl;
		cout<<"Vecteur conservatif"<<endl;
		*_runLogFile<<"Vecteur conservatif"<<endl;
		for(int k=0;k<_nVar;k++){
			cout<<Wcons[k]<<endl;
			*_runLogFile<<Wcons[k]<<endl;
		}
		_runLogFile->close();
		throw CdmathException("DriftModel::consToPrim: negative pressure");
	}

	Wprim[0]=concentration;
	Wprim[1] =  pression;
	for(int k=0;k<_Ndim;k++)
		Wprim[k+2] = Wcons[k+2]/Wcons[0];
	Wprim[_nVar-1] =  Temperature;
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

double DriftModel::getMixturePressure(double c_v, double rhom, double temperature)
{
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double Pinf_v=_fluides[0]->constante("p0");
	double Pinf_l=_fluides[1]->constante("p0");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");
	double c_l=1-c_v;
	double a=1., b, c;

	if(	!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);
		double e_v=fluide0->getInternalEnergy(temperature);
		double e_l=fluide1->getInternalEnergy(temperature);
		b=gamma_v*Pinf_v+gamma_l*Pinf_l -rhom*(c_l*(gamma_l-1)*(e_l-q_l)+c_v*(gamma_v-1)*(e_v-q_v));
		c=	gamma_v*Pinf_v*gamma_l*Pinf_l
				-rhom*(c_l*(gamma_l-1)*(e_l-q_l)*gamma_v*Pinf_v + c_v*(gamma_v-1)*(e_v-q_v)*gamma_l*Pinf_l);
	}
	else if(_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);
		double h_v=fluide0->getEnthalpy(temperature);
		double h_l=fluide1->getEnthalpy(temperature);
		b=Pinf_v+Pinf_l -rhom*(c_l*(gamma_l-1)*(h_l-q_l)/gamma_l+c_v*(gamma_v-1)*(h_v-q_v)/gamma_v);
		c=Pinf_v*Pinf_l-rhom*(c_l*(gamma_l-1)*(h_l-q_l)*Pinf_v + c_v*(gamma_v-1)*(h_v-q_v)*Pinf_l);
	}
	else
	{
		*_runLogFile<<"DriftModel::getMixturePressure: phases must have the same eos"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::getMixturePressure: phases must have the same eos");
	}

	double delta= b*b-4*a*c;

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"DriftModel::getMixturePressure: a= "<<a<<", b= "<<b<<", c= "<<c<<", delta= "<<delta<<endl;

	if(delta<0){
		*_runLogFile<<"DriftModel::getMixturePressure: cannot compute pressure, delta<0"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::getMixturePressure: cannot compute pressure, delta<0");
	}
	else
		return (-b+sqrt(delta))/(2*a);
}

void DriftModel::getMixturePressureAndTemperature(double c_v, double rhom, double rhom_em, double& pression, double& temperature)
{
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double Pinf_v=_fluides[0]->constante("p0");
	double Pinf_l=_fluides[1]->constante("p0");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");
	double c_l=1-c_v, m_v=c_v*rhom, m_l=rhom-m_v;
	double a, b, c, delta;

	if(	!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);

		temperature= (rhom_em-m_v*fluide0->getInternalEnergy(0)-m_l*fluide1->getInternalEnergy(0))
																																																																																																																																			/(m_v*fluide0->constante("cv")+m_l*fluide1->constante("cv"));

		double e_v=fluide0->getInternalEnergy(temperature);
		double e_l=fluide1->getInternalEnergy(temperature);
		a=1.;
		b=gamma_v*Pinf_v+gamma_l*Pinf_l -rhom*(c_l*(gamma_l-1)*(e_l-q_l)+c_v*(gamma_v-1)*(e_v-q_v));
		c=	gamma_v*Pinf_v*gamma_l*Pinf_l
				-rhom*(c_l*(gamma_l-1)*(e_l-q_l)*gamma_v*Pinf_v + c_v*(gamma_v-1)*(e_v-q_v)*gamma_l*Pinf_l);

		delta= b*b-4*a*c;

		if(delta<0){
			*_runLogFile<<"DriftModel::getMixturePressureAndTemperature: cannot compute pressure, delta<0"<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::getMixturePressureAndTemperature: cannot compute pressure, delta<0");
		}
		else
			pression= (-b+sqrt(delta))/(2*a);

	}
	else if(_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);

		double h0_v=fluide0->getEnthalpy(0);
		double h0_l=fluide1->getEnthalpy(0);
		double cp_v = _fluides[0]->constante("cp");
		double cp_l = _fluides[1]->constante("cp");
		double denom=m_v*cp_v+m_l*cp_l;
		double num_v=cp_v*rhom_em+m_l*(cp_l* h0_v-cp_v* h0_l);
		double num_l=cp_l*rhom_em+m_v*(cp_v* h0_l-cp_l* h0_v);

		a=gamma_v*gamma_l-(m_v*(gamma_v-1)*gamma_l*cp_v+m_l*(gamma_l-1)*gamma_v*cp_l)/denom;
		b=gamma_v*gamma_l*(Pinf_v+Pinf_l)-(m_v*(gamma_v-1)*gamma_l*cp_v*Pinf_l+m_l*(gamma_l-1)*gamma_v*cp_l*Pinf_v)/denom
				-m_v*(gamma_v-1)*gamma_l*(num_v/denom -q_v)-m_l*(gamma_l-1)*gamma_v*(num_l/denom -q_l);
		c=gamma_v*gamma_l*Pinf_v*Pinf_l
				-m_v*(gamma_v-1)*gamma_l*(num_v/denom -q_v)*Pinf_l-m_l*(gamma_l-1)*gamma_v*(num_l/denom -q_l)*Pinf_v;

		delta= b*b-4*a*c;

		if(delta<0)
		{
			*_runLogFile<<"DriftModel::getMixturePressureAndTemperature: cannot compute pressure, delta<0"<<endl;
			_runLogFile->close();
			throw CdmathException("DriftModel::getMixturePressureAndTemperature: cannot compute pressure, delta<0");
		}
		else
			pression= (-b+sqrt(delta))/(2*a);

		temperature=(rhom_em+pression-m_v*h0_v-m_l*h0_l)/denom;
	}
	else
	{
		_runLogFile->close();
		throw CdmathException("DriftModel::getMixturePressureAndTemperature: phases must have the same eos");
	}


	if(_verbose && _nbTimeStep%_freqSave ==0){
		cout<<"DriftModel::getMixturePressureAndTemperature: a= "<<a<<", b= "<<b<<", c= "<<c<<", delta= "<<delta<<endl;
		cout<<"pressure= "<<pression<<", temperature= "<<temperature<<endl;
	}

}
double DriftModel::getMixtureTemperature(double c_v, double rhom, double pression)
{
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double Pinf_v = _fluides[0]->constante("p0");
	double Pinf_l = _fluides[1]->constante("p0");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");
	double c_l=1-c_v;

	if(	!_useDellacherieEOS)
	{
		double cv_v = _fluides[0]->constante("cv");
		double cv_l = _fluides[1]->constante("cv");
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);
		double e0_v=fluide0->getInternalEnergy(0);
		double e0_l=fluide1->getInternalEnergy(0);

		double numerator=(pression+gamma_v*Pinf_v)*(pression+gamma_l*Pinf_l)/rhom
				- c_l*(pression+gamma_v*Pinf_v)*(gamma_l-1)*(e0_l-q_v)
				- c_v*(pression+gamma_l*Pinf_l)*(gamma_v-1)*(e0_v-q_l);
		double denominator=  c_l*(pression+gamma_v*Pinf_v)*(gamma_l-1)*cv_l
				+c_v*(pression+gamma_l*Pinf_l)*(gamma_v-1)*cv_v;
		return numerator/denominator;
	}
	else if(_useDellacherieEOS)
	{
		double cp_v = _fluides[0]->constante("cp");
		double cp_l = _fluides[1]->constante("cp");
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);
		double h0_v=fluide0->getEnthalpy(0);
		double h0_l=fluide1->getEnthalpy(0);

		double numerator= gamma_v*(pression+Pinf_v)* gamma_l*(pression+Pinf_l)/rhom
				- c_l*gamma_v*(pression+Pinf_v)*(gamma_l-1)*(h0_l-q_l)
				- c_v*gamma_l*(pression+Pinf_l)*(gamma_v-1)*(h0_v-q_v);
		double denominator=  c_l*gamma_v*(pression+Pinf_v)*(gamma_l-1)*cp_l
				+c_v*gamma_l*(pression+Pinf_l)*(gamma_v-1)*cp_v;
		return numerator/denominator;
	}
	else
	{
		*_runLogFile<<"DriftModel::getMixtureTemperature: phases must have the same eos"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::getMixtureTemperature: phases must have the same eos");
	}

}

void DriftModel::getMixturePressureDerivatives(double m_v, double m_l, double pression, double Tm)
{
	double rho_v=_fluides[0]->getDensity(pression,Tm);
	double rho_l=_fluides[1]->getDensity(pression,Tm);
	double alpha_v=m_v/rho_v,alpha_l=1-alpha_v;
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");
	double temp1, temp2, denom;

	if(	   !_useDellacherieEOS)
	{//Classical stiffened gas with linear law e(T)
		double cv_v = _fluides[0]->constante("cv");
		double cv_l = _fluides[1]->constante("cv");

		double e_v=_fluides[0]->getInternalEnergy(Tm, rho_v);//Actually does not depends on rho_v
		double e_l=_fluides[1]->getInternalEnergy(Tm, rho_l);//Actually does not depends on rho_l

		//On estime les dérivées discrètes de la pression (cf doc)
		temp1= m_v* cv_v + m_l* cv_l;
		denom=(alpha_v/((gamma_v-1)*rho_v*(e_v-q_v))+alpha_l/((gamma_l-1)*rho_l*(e_l-q_l)))*temp1;
		temp2=alpha_v*cv_v/ (e_v-q_v)+alpha_l*cv_l/ (e_l-q_l);

		_khi=(temp1/rho_l-e_l*temp2)/denom;
		_ksi=(temp1*(rho_l-rho_v)/(rho_v*rho_l)+(e_l-e_v)*temp2)/denom;
		_kappa=temp2/denom;
	}

	else if( _useDellacherieEOS)
	{//S. Dellacherie stiffened gas with linear law h(T)
		double cp_v = _fluides[0]->constante("cp");
		double cp_l = _fluides[1]->constante("cp");

		double h_v=_fluides[0]->getEnthalpy(Tm, rho_v);//Actually does not depends on rho_v
		double h_l=_fluides[1]->getEnthalpy(Tm, rho_l);//Actually does not depends on rho_l

		//On estime les dérivées discrètes de la pression (cf doc)
		temp1= m_v* cp_v + m_l* cp_l;
		temp2= alpha_v*cp_v/(h_v-q_v)+alpha_l*cp_l/(h_l-q_l);
		//denom=temp1*(alpha_v/(pression+Pinf_v)+alpha_l/(pression+Pinf_l))-temp2;
		denom=temp1*(alpha_v*gamma_v/((gamma_v-1)*rho_v*(h_v-q_v))+alpha_l*gamma_l/((gamma_l-1)*rho_l*(h_l-q_l)))-temp2;

		//On estime les dérivées discrètes de la pression (cf doc)
		_khi=(temp1/rho_l - h_l*temp2)/denom;
		_ksi=(temp1*(1/rho_v-1/rho_l)+(h_l-h_v)*temp2)/denom;
		_kappa=temp2/denom;
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"rho_l= "<<rho_l<<", temp1= "<<temp1<<", temp2= "<<temp2<<", denom= "<<denom<<endl;
		cout<<"khi= "<<_khi<<", ksi= "<<_ksi<<", kappa= "<<_kappa<<endl;
	}
}

void DriftModel::entropicShift(double* n)//TO do: make sure _Vi and _Vj are well set
{
	//_Vi is (cm, p, u, T)
	//_Ui is (rhom, rhom cm, rhom um, rhom Em)
	consToPrim(_Ui,_Vi,_porosityi);
	double umnl=0, ul_2=0, umnr=0, ur_2=0; //valeur de u.normale et |u|²
	for(int i=0;i<_Ndim;i++)
	{
		ul_2 += _Vi[2+i]*_Vi[2+i];
		umnl += _Vi[2+i]*n[i];//vitesse normale
		ur_2 += _Vj[2+i]*_Vj[2+i];
		umnr += _Vj[2+i]*n[i];//vitesse normale
	}

	//Left sound speed
	double rhom=_Ui[0];
	double cm=_Vi[0];
	double Hm=(_Ui[_nVar-1]+_Vi[1])/rhom;
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Entropic shift left state rhom="<<rhom<<" cm= "<<cm<<"Hm= "<<Hm<<endl;
	double Tm=_Vi[_nVar-1];
	double hm=Hm-0.5*ul_2;
	double m_v=cm*rhom, m_l=rhom-m_v;
	double pression=getMixturePressure( cm, rhom, Tm);

	getMixturePressureDerivatives( m_v, m_l, pression, Tm);

	if(_kappa*hm+_khi+cm*_ksi<0)
	{
		*_runLogFile<<"DriftModel::entropicShift : vitesse du son cellule gauche complexe"<<endl;
		_runLogFile->close();
		throw CdmathException("Valeurs propres jacobienne: vitesse du son complexe");
	}
	double aml=sqrt(_kappa*hm+_khi+cm*_ksi);//vitesse du son du melange
	//cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", am= "<<am<<endl;

	//Right sound speed
	consToPrim(_Uj,_Vj,_porosityj);
	rhom=_Uj[0];
	cm=_Vj[0];
	Hm=(_Uj[_nVar-1]+_Vj[1])/rhom;
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Entropic shift right state rhom="<<rhom<<" cm= "<<cm<<"Hm= "<<Hm<<endl;
	Tm=_Vj[_nVar-1];
	hm=Hm-0.5*ul_2;
	m_v=cm*rhom;
	m_l=rhom-m_v;
	pression=getMixturePressure( cm, rhom, Tm);

	getMixturePressureDerivatives( m_v, m_l, pression, Tm);

	if(_kappa*hm+_khi+cm*_ksi<0){
		*_runLogFile<<"DriftModel::entropicShift: vitesse du son cellule droite complexe"<<endl;
		_runLogFile->close();
		throw CdmathException("Valeurs propres jacobienne: vitesse du son complexe");
	}

	double amr=sqrt(_kappa*hm+_khi+cm*_ksi);//vitesse du son du melange
	//cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", am= "<<am<<endl;

	_entropicShift[0]=abs(umnl-aml - (umnr-amr));
	_entropicShift[1]=abs(umnl - umnr);
	_entropicShift[2]=abs(umnl+aml - (umnr+amr));

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Entropic shift left eigenvalues: "<<endl;
		cout<<"("<< umnl-aml<< ", " <<umnl<<", "<<umnl+aml << ")";
		cout<<endl;
		cout<<"Entropic shift right eigenvalues: "<<endl;
		cout<<"("<< umnr-amr<< ", " <<umnr<<", "<<umnr+amr << ")";
		cout<< endl;
		cout<<"eigenvalue jumps "<<endl;
		cout<< _entropicShift[0] << ", " << _entropicShift[1] << ", "<< _entropicShift[2] <<endl;
	}
}

Vector DriftModel::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"DriftModel::convectionFlux start"<<endl;
		cout<<"Ucons"<<endl;
		cout<<U<<endl;
		cout<<"Vprim"<<endl;
		cout<<V<<endl;
	}

	double rho_m=U(0);//phi rhom
	double m_v=U(1), m_l=rho_m-m_v;//phi rhom cv, phi rhom cl
	Vector q_m(_Ndim);
	for(int i=0;i<_Ndim;i++)
		q_m(i)=U(2+i);

	double concentration_v=V(0);
	double concentration_l= 1 - concentration_v;
	double pression=V(1);
	Vector vitesse_melange(_Ndim);
	for(int i=0;i<_Ndim;i++)
		vitesse_melange(i)=V(2+i);
	double Temperature= V(2+_Ndim);//(rho_m_e_m-m_v*internal_energy(0, e0v,c0v,T0)-m_l*internal_energy(0, e0l,c0l,T0))/(m_v*c0v+m_l*c0l);

	Vector vitesse_relative=relative_velocity(concentration_v,vitesse_melange,rho_m);
	Vector vitesse_v=vitesse_melange+concentration_l*vitesse_relative;
	Vector vitesse_l=vitesse_melange-concentration_v*vitesse_relative;
	double vitesse_vn=vitesse_v*normale;
	double vitesse_ln=vitesse_l*normale;
	//double rho_m_e_m=rho_m_E_m -0.5*(q_m*vitesse_melange + rho_m*concentration_v*concentration_l*vitesse_relative*vitesse_relative)
	double rho_v=_fluides[0]->getDensity(pression,Temperature);
	double rho_l=_fluides[1]->getDensity(pression,Temperature);
	double e_v=_fluides[0]->getInternalEnergy(Temperature,rho_v);
	double e_l=_fluides[1]->getInternalEnergy(Temperature,rho_l);

	Vector F(_nVar);
	F(0)=m_v*vitesse_vn+m_l*vitesse_ln;
	F(1)=m_v*vitesse_vn;
	for(int i=0;i<_Ndim;i++)
		F(2+i)=m_v*vitesse_vn*vitesse_v(i)+m_l*vitesse_ln*vitesse_l(i)+pression*normale(i)*porosity;
	F(2+_Ndim)=m_v*(e_v+0.5*vitesse_v*vitesse_v+pression/rho_v)*vitesse_vn+m_l*(e_l+0.5*vitesse_l*vitesse_l+pression/rho_l)*vitesse_ln;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"DriftModel::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

Vector DriftModel::staggeredVFFCFlux()
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"DriftModel::staggeredVFFCFlux()start"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<<"DriftModel::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered
	{
		Vector Fij(_nVar);

		double uijn=0, phiqn=0, uin=0, ujn=0;
		for(int idim=0;idim<_Ndim;idim++)
		{
			uijn+=_vec_normal[idim]*_Uroe[2+idim];//URoe = rho, cm, u, H
			uin+=_vec_normal[idim]*_Ui[2+idim];
			ujn+=_vec_normal[idim]*_Uj[2+idim];
		}
		if( (uin>0 && ujn >0) || (uin>=0 && ujn <=0 && uijn>0) ) // formerly (uijn>_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Ui[2+idim];//phi rho u n
			Fij(0)=phiqn;
			Fij(1)=_Vi[0]*phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+idim)=phiqn*_Vi[2+idim]+_Vj[1]*_vec_normal[idim]*_porosityj;
			Fij(_nVar-1)=phiqn/_Ui[0]*(_Ui[_nVar-1]+_Vj[1]*sqrt(_porosityj/_porosityi));
		}
		else if( (uin<0 && ujn <0) || (uin>=0 && ujn <=0 && uijn<0) ) // formerly (uijn<-_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Uj[2+idim];//phi rho u n
			Fij(0)=phiqn;
			Fij(1)=_Vj[0]*phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+idim)=phiqn*_Vj[2+idim]+_Vi[1]*_vec_normal[idim]*_porosityi;
			Fij(_nVar-1)=phiqn/_Uj[0]*(_Uj[_nVar-1]+_Vi[1]*sqrt(_porosityi/_porosityj));
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
			Fij=(Fi+Fj)/2+_maxvploc*(Ui-Uj)/2;

			//case nil velocity on the interface, apply smoothed scheme
			/*
			double theta=uijn/(.01);
			Vector F1(_nVar),F2(_nVar);
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Ui[2+idim];//phi rho u n
			F1(0)=phiqn;
			F1(1)=_Vi[0]*phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				F1(2+idim)=phiqn*_Vi[2+idim]+_Vj[1]*_vec_normal[idim]*_porosityj;
			F1(_nVar-1)=phiqn/_Ui[0]*(_Ui[_nVar-1]+_Vj[1]*sqrt(_porosityj/_porosityi));

			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Uj[2+idim];//phi rho u n
			F2(0)=phiqn;
			F2(1)=_Vj[0]*phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				F2(2+idim)=phiqn*_Vj[2+idim]+_Vi[1]*_vec_normal[idim]*_porosityi;
			F2(_nVar-1)=phiqn/_Uj[0]*(_Uj[_nVar-1]+_Vi[1]*sqrt(_porosityi/_porosityj));

			Fij=(1+theta)/2*F1+(1-theta)/2*F2;
			 */
		}
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"DriftModel::staggeredVFFCFlux() end uijn="<<uijn<<endl;
			cout<<Fij<<endl;
		}
		return Fij;
	}
}

void DriftModel::staggeredVFFCMatricesConservativeVariables(double u_mn)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"DriftModel::staggeredVFFCMatricesConservativeVariables()"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<<"DriftModel::staggeredVFFCMatricesConservativeVariables: staggeredVFFCMatrices method should be called only for VFFC formulation and staggered upwinding"<<endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::staggeredVFFCMatricesConservativeVariables: staggeredVFFCMatrices method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		double ui_n=0, ui_2=0, uj_n=0, uj_2=0;//vitesse normale et carré du module
		double H;//enthalpie totale (expression particulière)
		consToPrim(_Ui,_Vi,_porosityi);
		consToPrim(_Uj,_Vj,_porosityj);

		for(int i=0;i<_Ndim;i++)
		{
			ui_2 += _Vi[2+i]*_Vi[2+i];
			ui_n += _Vi[2+i]*_vec_normal[i];
			uj_2 += _Vj[2+i]*_Vj[2+i];
			uj_n += _Vj[2+i]*_vec_normal[i];
		}

		double rhomi=_Ui[0]/_porosityi;
		double mi_v=_Ui[1]/_porosityi;
		double mi_l=rhomi-mi_v;
		double cmi=_Vi[0];
		double pi=_Vi[1];
		double Tmi=_Vi[_Ndim+2];
		double Emi=_Ui[_Ndim+2]/(rhomi*_porosityi);
		double ecini=0.5*ui_2;
		double emi=Emi-0.5*ui_2;
		double hmi=emi+pi/rhomi;

		double pj=_Vj[1];
		double rhomj=_Uj[0]/_porosityj;
		double mj_v =_Uj[1]/_porosityj;
		double mj_l=rhomj-mj_v;
		double cmj=_Vj[0];
		double Tmj=_Vj[_Ndim+2];
		double Emj=_Uj[_Ndim+2]/(rhomj*_porosityj);
		double ecinj=0.5*uj_2;
		double emj=Emj-0.5*uj_2;
		double hmj=emj+pj/rhomj;

		double Hm;//value depends on the sign of interfacial velocity

		if(u_mn>_precision)
		{
			Hm=Emi+pj/rhomi;
			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"VFFC Staggered state rhomi="<<rhomi<<" cmi= "<<cmi<<" Hm= "<<Hm<<" Tmi= "<<Tmi<<" pj= "<<pj<<endl;

			/***********Calcul des valeurs propres ********/
			vector<complex<double> > vp_dist(3,0);
			getMixturePressureDerivatives( mj_v, mj_l, pj, Tmj);
			if(_kappa*hmj+_khi+cmj*_ksi<0){
				*_runLogFile<<"staggeredVFFCMatricesConservativeVariables: vitesse du son complexe"<<endl;
				_runLogFile->close();
				throw CdmathException("staggeredVFFCMatricesConservativeVariables: vitesse du son complexe");
			}
			double amj=sqrt(_kappa*hmj+_khi+cmj*_ksi);//vitesse du son du melange

			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", amj= "<<amj<<endl;

			//On remplit les valeurs propres
			vp_dist[0]=ui_n+amj;
			vp_dist[1]=ui_n-amj;
			vp_dist[2]=ui_n;

			_maxvploc=fabs(ui_n)+amj;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			/******** Construction de la matrice A^+ *********/
			_AroePlus[0*_nVar+0]=0;
			_AroePlus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[0*_nVar+2+i]=_vec_normal[i];
			_AroePlus[0*_nVar+2+_Ndim]=0;
			_AroePlus[1*_nVar+0]=-ui_n*cmi;
			_AroePlus[1*_nVar+1]=ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[1*_nVar+2+i]=cmi*_vec_normal[i];
			_AroePlus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroePlus[(2+i)*_nVar+0]=-ui_n*_Vi[2+i];
				_AroePlus[(2+i)*_nVar+1]=0;
				for(int j=0;j<_Ndim;j++)
					_AroePlus[(2+i)*_nVar+2+j]=_Vi[2+i]*_vec_normal[j];
				_AroePlus[(2+i)*_nVar+2+i]+=ui_n;
				_AroePlus[(2+i)*_nVar+2+_Ndim]=0;
			}
			_AroePlus[(2+_Ndim)*_nVar+0]=-Hm*ui_n;
			_AroePlus[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i];
			_AroePlus[(2+_Ndim)*_nVar+2+_Ndim]=ui_n;

			/******** Construction de la matrice A^- *********/
			_AroeMinus[0*_nVar+0]=0;
			_AroeMinus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[0*_nVar+2+i]=0;
			_AroeMinus[0*_nVar+2+_Ndim]=0;
			_AroeMinus[1*_nVar+0]=0;
			_AroeMinus[1*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[1*_nVar+2+i]=0;
			_AroeMinus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroeMinus[(2+i)*_nVar+0]=(_khi+_kappa*ecinj)*_vec_normal[i];
				_AroeMinus[(2+i)*_nVar+1]=_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroeMinus[(2+i)*_nVar+2+j]=-_kappa*_vec_normal[i]*_Vj[2+j];
				_AroeMinus[(2+i)*_nVar+2+i]+=0;
				_AroeMinus[(2+i)*_nVar+2+_Ndim]=_kappa*_vec_normal[i];
			}
			_AroeMinus[(2+_Ndim)*_nVar+0]=(_khi+_kappa*ecinj)*ui_n;
			_AroeMinus[(2+_Ndim)*_nVar+1]=_ksi*ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[(2+_Ndim)*_nVar+2+i]=-_kappa*uj_n*_Vi[2+i];
			_AroeMinus[(2+_Ndim)*_nVar+2+_Ndim]=_kappa*ui_n;
		}
		else if(u_mn<-_precision)
		{
			Hm=Emj+pi/rhomj;
			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"VFFC Staggered state rhomj="<<rhomj<<" cmj= "<<cmj<<" Hm= "<<Hm<<" Tmj= "<<Tmj<<" pi= "<<pi<<endl;

			/***********Calcul des valeurs propres ********/
			vector<complex<double> > vp_dist(3,0);
			getMixturePressureDerivatives( mi_v, mi_l, pi, Tmi);
			if(_kappa*hmi+_khi+cmi*_ksi<0)
			{
				*_runLogFile<<"staggeredVFFCMatricesConservativeVariables: vitesse du son complexe"<<endl;
				_runLogFile->close();
				throw CdmathException("staggeredVFFCMatricesConservativeVariables: vitesse du son complexe");
			}
			double ami=sqrt(_kappa*hmi+_khi+cmi*_ksi);//vitesse du son du melange
			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", ami= "<<ami<<endl;

			//On remplit les valeurs propres
			vp_dist[0]=uj_n+ami;
			vp_dist[1]=uj_n-ami;
			vp_dist[2]=uj_n;

			_maxvploc=fabs(uj_n)+ami;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			/******** Construction de la matrice A^+ *********/
			_AroePlus[0*_nVar+0]=0;
			_AroePlus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[0*_nVar+2+i]=0;
			_AroePlus[0*_nVar+2+_Ndim]=0;
			_AroePlus[1*_nVar+0]=0;
			_AroePlus[1*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[1*_nVar+2+i]=0;
			_AroePlus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroePlus[(2+i)*_nVar+0]=(_khi+_kappa*ecini)*_vec_normal[i];
				_AroePlus[(2+i)*_nVar+1]=_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroePlus[(2+i)*_nVar+2+j]=-_kappa*_vec_normal[i]*_Vi[2+j];
				_AroePlus[(2+i)*_nVar+2+i]+=0;
				_AroePlus[(2+i)*_nVar+2+_Ndim]=_kappa*_vec_normal[i];
			}
			_AroePlus[(2+_Ndim)*_nVar+0]=(_khi+_kappa*ecini)*ui_n;
			_AroePlus[(2+_Ndim)*_nVar+1]=_ksi*ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[(2+_Ndim)*_nVar+2+i]=-_kappa*uj_n*_Vj[2+i];
			_AroePlus[(2+_Ndim)*_nVar+2+_Ndim]=_kappa*ui_n;

			/******** Construction de la matrice A^- *********/
			_AroeMinus[0*_nVar+0]=0;
			_AroeMinus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[0*_nVar+2+i]=_vec_normal[i];
			_AroeMinus[0*_nVar+2+_Ndim]=0;
			_AroeMinus[1*_nVar+0]=-uj_n*cmj;
			_AroeMinus[1*_nVar+1]=uj_n;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[1*_nVar+2+i]=cmj*_vec_normal[i];
			_AroeMinus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroeMinus[(2+i)*_nVar+0]=-uj_n*_Vj[2+i];
				_AroeMinus[(2+i)*_nVar+1]=0;
				for(int j=0;j<_Ndim;j++)
					_AroeMinus[(2+i)*_nVar+2+j]=_Vj[2+i]*_vec_normal[j];
				_AroeMinus[(2+i)*_nVar+2+i]+=uj_n;
				_AroeMinus[(2+i)*_nVar+2+_Ndim]=0;
			}
			_AroeMinus[(2+_Ndim)*_nVar+0]=-Hm*uj_n;
			_AroeMinus[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i];
			_AroeMinus[(2+_Ndim)*_nVar+2+_Ndim]=uj_n;
		}
		else//case nil velocity on the interface
		{
			double rhom=_Uroe[0];
			double cm=_Uroe[1];
			double Hm=_Uroe[_nVar-1];
			double m_v=cm*rhom, m_l=rhom-m_v;
			double Tm;
			double umn=0,u_2=0;
			Vector vitesse(_Ndim);
			for(int idim=0;idim<_Ndim;idim++)
			{
				vitesse[idim]=_Uroe[2+idim];
				umn += _Uroe[2+idim]*_vec_normal[idim];//vitesse normale
				u_2 += _Uroe[2+idim]*_Uroe[2+idim];
			}
			double hm=Hm-0.5*u_2;

			if(cm<_precision)//pure liquid
				Tm=_fluides[1]->getTemperatureFromEnthalpy(hm,rhom);
			else if(cm>1-_precision)
				Tm=_fluides[0]->getTemperatureFromEnthalpy(hm,rhom);
			else//Hypothèse de saturation
				Tm=_Tsat;

			double pression= getMixturePressure(cm, rhom, Tm);

			getMixturePressureDerivatives( m_v, m_l, pression, Tm);//EOS is involved to express pressure jump and sound speed
			if(_kappa*hm+_khi+cm*_ksi<0){
				*_runLogFile<<"Calcul matrice de Roe: vitesse du son complexe"<<endl;
				_runLogFile->close();
				throw CdmathException("Calcul matrice de Roe: vitesse du son complexe");
			}
			double am=sqrt(_kappa*hm+_khi+cm*_ksi);//vitesse du son du melange
			_maxvploc=fabs(umn)+am;
			if(_maxvploc>_maxvp)
				_maxvp=_maxvploc;

			//Rusanov scheme
			/******* Construction de la matrice A^+ ********/
			//A^+=0.5*jacobienne Flux (U_i)+0.5 _maxvploc Id
			Hm=Emi+pi/rhomi;
			getMixturePressureDerivatives( mi_v, mi_l, pi, Tmi);
			_AroePlus[0*_nVar+0]=0;
			_AroePlus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[0*_nVar+2+i]=0.5*_vec_normal[i];
			_AroePlus[0*_nVar+2+_Ndim]=0;
			_AroePlus[1*_nVar+0]=-0.5*ui_n*cmi;
			_AroePlus[1*_nVar+1]=0.5*ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[1*_nVar+2+i]=0.5*cmi*_vec_normal[i];
			_AroePlus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroePlus[(2+i)*_nVar+0]=0.5*(_khi+_kappa*ecini)*_vec_normal[i]-0.5*ui_n*_Vi[2+i];
				_AroePlus[(2+i)*_nVar+1]=0.5*_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroePlus[(2+i)*_nVar+2+j]=0.5*_Vi[2+i]*_vec_normal[j]-0.5*_kappa*_vec_normal[i]*_Vi[2+j];
				_AroePlus[(2+i)*_nVar+2+i]+=0.5*ui_n;
				_AroePlus[(2+i)*_nVar+2+_Ndim]=0.5*_kappa*_vec_normal[i];
			}
			_AroePlus[(2+_Ndim)*_nVar+0]=-0.5*Hm*ui_n;
			_AroePlus[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlus[(2+_Ndim)*_nVar+2+i]=0.5*Hm*_vec_normal[i];
			_AroePlus[(2+_Ndim)*_nVar+2+_Ndim]=0.5*ui_n;

			for(int i=0;i<_nVar;i++)
				_AroePlus[i*_nVar+i]+=_maxvploc/2;
			/****** Construction de la matrice A^- *******/
			//A^-=0.5*jacobienne Flux (U_j)-0.5 _maxvploc Id
			Hm=Emj+pj/rhomj;
			getMixturePressureDerivatives( mj_v, mj_l, pj, Tmj);
			_AroeMinus[0*_nVar+0]=0;
			_AroeMinus[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[0*_nVar+2+i]=0.5*_vec_normal[i];
			_AroeMinus[0*_nVar+2+_Ndim]=0;
			_AroeMinus[1*_nVar+0]=-0.5*uj_n*cmj;
			_AroeMinus[1*_nVar+1]=0.5*uj_n;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[1*_nVar+2+i]=0.5*cmj*_vec_normal[i];
			_AroeMinus[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroeMinus[(2+i)*_nVar+0]=0.5*(_khi+_kappa*ecinj)*_vec_normal[i]-0.5*uj_n*_Vj[2+i];
				_AroeMinus[(2+i)*_nVar+1]=0.5*_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroeMinus[(2+i)*_nVar+2+j]=0.5*_Vj[2+i]*_vec_normal[j]-0.5*_kappa*_vec_normal[i]*_Vj[2+j];
				_AroeMinus[(2+i)*_nVar+2+i]+=0.5*uj_n;
				_AroeMinus[(2+i)*_nVar+2+_Ndim]=0.5*_kappa*_vec_normal[i];
			}
			_AroeMinus[(2+_Ndim)*_nVar+0]=-0.5*Hm*uj_n;
			_AroeMinus[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinus[(2+_Ndim)*_nVar+2+i]=0.5*Hm*_vec_normal[i];
			_AroeMinus[(2+_Ndim)*_nVar+2+_Ndim]=0.5*uj_n;

			for(int i=0;i<_nVar;i++)
				_AroeMinus[i*_nVar+i]-=_maxvploc/2;

			/*
			double _AroePlusi[_nVar*_nVar], _AroePlusj[_nVar*_nVar];
			double _AroeMinusi[_nVar*_nVar], _AroeMinusj[_nVar*_nVar];

			//Matrices côté gauche
			Hm=Emi+pj/rhomi;
			 ****** Construction de la matrice A^+ i (gauche) *******
			_AroePlusi[0*_nVar+0]=0;
			_AroePlusi[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlusi[0*_nVar+2+i]=_vec_normal[i];
			_AroePlusi[0*_nVar+2+_Ndim]=0;
			_AroePlusi[1*_nVar+0]=-ui_n*cmi;
			_AroePlusi[1*_nVar+1]=ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroePlusi[1*_nVar+2+i]=cmi*_vec_normal[i];
			_AroePlusi[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroePlusi[(2+i)*_nVar+0]=-ui_n*_Vi[2+i];
				_AroePlusi[(2+i)*_nVar+1]=0;
				for(int j=0;j<_Ndim;j++)
					_AroePlusi[(2+i)*_nVar+2+j]=_Vi[2+i]*_vec_normal[j];
				_AroePlusi[(2+i)*_nVar+2+i]+=ui_n;
				_AroePlusi[(2+i)*_nVar+2+_Ndim]=0;
			}
			_AroePlusi[(2+_Ndim)*_nVar+0]=-Hm*ui_n;
			_AroePlusi[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlusi[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i];
			_AroePlusi[(2+_Ndim)*_nVar+2+_Ndim]=ui_n;

			 ****** Construction de la matrice A^- i (coté gauche)*******
			_AroeMinusi[0*_nVar+0]=0;
			_AroeMinusi[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusi[0*_nVar+2+i]=0;
			_AroeMinusi[0*_nVar+2+_Ndim]=0;
			_AroeMinusi[1*_nVar+0]=0;
			_AroeMinusi[1*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusi[1*_nVar+2+i]=0;
			_AroeMinusi[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroeMinusi[(2+i)*_nVar+0]=(_khi+_kappa*ecinj)*_vec_normal[i];
				_AroeMinusi[(2+i)*_nVar+1]=_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroeMinusi[(2+i)*_nVar+2+j]=-_kappa*_vec_normal[i]*_Vj[2+j];
				_AroeMinusi[(2+i)*_nVar+2+i]+=0;
				_AroeMinusi[(2+i)*_nVar+2+_Ndim]=_kappa*_vec_normal[i];
			}
			_AroeMinusi[(2+_Ndim)*_nVar+0]=(_khi+_kappa*ecinj)*ui_n;
			_AroeMinusi[(2+_Ndim)*_nVar+1]=_ksi*ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusi[(2+_Ndim)*_nVar+2+i]=-_kappa*uj_n*_Vi[2+i];
			_AroeMinusi[(2+_Ndim)*_nVar+2+_Ndim]=_kappa*ui_n;

			//Matrices côté droit
			Hm=Emj+pi/rhomj;

			 ****** Construction de la matrice A^+ j (coté droit)*******
			_AroePlusj[0*_nVar+0]=0;
			_AroePlusj[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlusj[0*_nVar+2+i]=0;
			_AroePlusj[0*_nVar+2+_Ndim]=0;
			_AroePlusj[1*_nVar+0]=0;
			_AroePlusj[1*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroePlusj[1*_nVar+2+i]=0;
			_AroePlusj[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroePlusj[(2+i)*_nVar+0]=(_khi+_kappa*ecini)*_vec_normal[i];
				_AroePlusj[(2+i)*_nVar+1]=_ksi*_vec_normal[i];
				for(int j=0;j<_Ndim;j++)
					_AroePlusj[(2+i)*_nVar+2+j]=-_kappa*_vec_normal[i]*_Vi[2+j];
				_AroePlusj[(2+i)*_nVar+2+i]+=0;
				_AroePlusj[(2+i)*_nVar+2+_Ndim]=_kappa*_vec_normal[i];
			}
			_AroePlusj[(2+_Ndim)*_nVar+0]=(_khi+_kappa*ecini)*ui_n;
			_AroePlusj[(2+_Ndim)*_nVar+1]=_ksi*ui_n;
			for(int i=0;i<_Ndim;i++)
				_AroePlusj[(2+_Ndim)*_nVar+2+i]=-_kappa*uj_n*_Vj[2+i];
			_AroePlusj[(2+_Ndim)*_nVar+2+_Ndim]=_kappa*ui_n;

			 ****** Construction de la matrice A^- j (coté droit)*******
			_AroeMinusj[0*_nVar+0]=0;
			_AroeMinusj[0*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusj[0*_nVar+2+i]=_vec_normal[i];
			_AroeMinusj[0*_nVar+2+_Ndim]=0;
			_AroeMinusj[1*_nVar+0]=-uj_n*cmj;
			_AroeMinusj[1*_nVar+1]=uj_n;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusj[1*_nVar+2+i]=cmj*_vec_normal[i];
			_AroeMinusj[1*_nVar+2+_Ndim]=0;
			for(int i=0;i<_Ndim;i++)
			{
				_AroeMinusj[(2+i)*_nVar+0]=-uj_n*_Vj[2+i];
				_AroeMinusj[(2+i)*_nVar+1]=0;
				for(int j=0;j<_Ndim;j++)
					_AroeMinusj[(2+i)*_nVar+2+j]=_Vj[2+i]*_vec_normal[j];
				_AroeMinusj[(2+i)*_nVar+2+i]+=uj_n;
				_AroeMinusj[(2+i)*_nVar+2+_Ndim]=0;
			}
			_AroeMinusj[(2+_Ndim)*_nVar+0]=-Hm*uj_n;
			_AroeMinusj[(2+_Ndim)*_nVar+1]=0;
			for(int i=0;i<_Ndim;i++)
				_AroeMinusj[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i];
			_AroeMinusj[(2+_Ndim)*_nVar+2+_Ndim]=uj_n;

			double theta=u_mn/(.01);
			for(int i=0; i<_nVar*_nVar;i++)
			{
				_AroePlus[i] =(1+theta)/2*_AroePlusi[i] +(1-theta)/2*_AroePlusj[i];
				_AroeMinus[i]=(1+theta)/2*_AroeMinusi[i]+(1-theta)/2*_AroeMinusj[i];
			}
			 */
		}
	}
	for(int i=0; i<_nVar*_nVar;i++)
	{
		_Aroe[i]=_AroePlus[i]+_AroeMinus[i];
		_absAroe[i]=_AroePlus[i]-_AroeMinus[i];
	}

	if(_timeScheme==Implicit)
		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinusImplicit[i] = _AroeMinus[i];
			_AroePlusImplicit[i]  = _AroePlus[i];
		}
}

void DriftModel::staggeredVFFCMatricesPrimitiveVariables(double u_mn)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"DriftModel::staggeredVFFCMatricesPrimitiveVariables()"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<< "DriftModel::staggeredVFFCMatricesPrimitiveVariables: staggeredVFFCMatricesPrimitiveVariables method should be called only for VFFC formulation and staggered upwinding" << endl;
		_runLogFile->close();
		throw CdmathException("DriftModel::staggeredVFFCMatricesPrimitiveVariables: staggeredVFFCMatricesPrimitiveVariables method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		//Calls to getDensityDerivatives needed
		double ui_n=0, ui_2=0, uj_n=0, uj_2=0;//vitesse normale et carré du module
		double H;//enthalpie totale (expression particulière)
		consToPrim(_Ui,_Vi,_porosityi);
		consToPrim(_Uj,_Vj,_porosityj);

		for(int i=0;i<_Ndim;i++)
		{
			ui_2 += _Vi[2+i]*_Vi[2+i];
			ui_n += _Vi[2+i]*_vec_normal[i];
			uj_2 += _Vj[2+i]*_Vj[2+i];
			uj_n += _Vj[2+i]*_vec_normal[i];
		}

		double rhomi=_Ui[0]/_porosityi;
		double mi_v=_Ui[1]/_porosityi;
		double mi_l=rhomi-mi_v;
		double cmi=_Vi[0];
		double pi=_Vi[1];
		double Tmi=_Vi[_Ndim+2];
		double Emi=_Ui[_Ndim+2]/(rhomi*_porosityi);
		double ecini=0.5*ui_2;
		double emi=Emi-0.5*ui_2;
		double hmi=emi+pi/rhomi;
		double rho_vi=_fluides[0]->getDensity(pi,Tmi);
		double rho_li=_fluides[1]->getDensity(pi,Tmi);

		double pj=_Vj[1];
		double rhomj=_Uj[0]/_porosityj;
		double mj_v =_Uj[1]/_porosityj;
		double mj_l=rhomj-mj_v;
		double cmj=_Vj[0];
		double Tmj=_Vj[_Ndim+2];
		double Emj=_Uj[_Ndim+2]/(rhomj*_porosityj);
		double ecinj=0.5*uj_2;
		double emj=Emj-0.5*uj_2;
		double hmj=emj+pj/rhomj;
		double rho_vj=_fluides[0]->getDensity(pj,Tmj);
		double rho_lj=_fluides[1]->getDensity(pj,Tmj);

		double gamma_v=_fluides[0]->constante("gamma");
		double gamma_l=_fluides[1]->constante("gamma");
		double q_v=_fluides[0]->constante("q");
		double q_l=_fluides[1]->constante("q");

		if(fabs(u_mn)>_precision)//non zero velocity on the interface
		{
			if(		!_useDellacherieEOS)
			{
				StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
				StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);
				double cv_v=fluide0->constante("cv");
				double cv_l=fluide1->constante("cv");

				double e_vi=_fluides[0]->getInternalEnergy(Tmi,rho_vi);
				double e_li=_fluides[1]->getInternalEnergy(Tmi,rho_li);
				double e_vj=_fluides[0]->getInternalEnergy(Tmj,rho_vj);
				double e_lj=_fluides[1]->getInternalEnergy(Tmj,rho_lj);

				if(u_mn>_precision)
				{
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"VFFC Staggered state rhomi="<<rhomi<<" cmi= "<<cmi<<" Emi= "<<Emi<<" Tmi= "<<Tmi<<" pj= "<<pj<<endl;

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					getMixturePressureDerivatives( mj_v, mj_l, pj, Tmj);
					if(_kappa*hmj+_khi+cmj*_ksi<0)
					{
						*_runLogFile<<"staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe"<<endl;
						_runLogFile->close();
						throw CdmathException("staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe");
					}
					double amj=sqrt(_kappa*hmj+_khi+cmj*_ksi);//vitesse du son du melange

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", amj= "<<amj<<endl;

					//On remplit les valeurs propres
					vp_dist[0]=ui_n+amj;
					vp_dist[1]=ui_n-amj;
					vp_dist[2]=ui_n;

					_maxvploc=fabs(ui_n)+amj;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					/******** Construction de la matrice A^+ *********/
					_AroePlusImplicit[0*_nVar+0]=-rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li);
					_AroePlusImplicit[0*_nVar+1]= rhomi*rhomi*ui_n*(cmi/(rho_vi*rho_vi*(gamma_v-1)*(e_vi-q_v))+(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(e_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[0*_nVar+2+i]=rhomi*_vec_normal[i];
					_AroePlusImplicit[0*_nVar+2+_Ndim]=-rhomi*rhomi*ui_n*(cv_v*cmi/(rho_vi*(e_vi-q_v))+cv_l*(1-cmi)/(rho_li*(e_li-q_l)));
					_AroePlusImplicit[1*_nVar+0]=rhomi*ui_n-cmi*rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li);
					_AroePlusImplicit[1*_nVar+1]=-cmi*rhomi*rhomi*ui_n*(cmi/(rho_vi*rho_vi*(gamma_v-1)*(e_vi-q_v))+(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(e_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[1*_nVar+2+i]=cmi*rhomi*_vec_normal[i];
					_AroePlusImplicit[1*_nVar+2+_Ndim]=-cmi*rhomi*rhomi*ui_n*(cv_v*cmi/(rho_vi*(e_vi-q_v))+cv_l*(1-cmi)/(rho_li*(e_li-q_l)));
					for(int i=0;i<_Ndim;i++)
					{
						_AroePlusImplicit[(2+i)*_nVar+0]=-rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li)*_Vi[2+i];
						_AroePlusImplicit[(2+i)*_nVar+1]=rhomi*rhomi*ui_n*(cmi/(rho_vi*rho_vi*(gamma_v-1)*(e_vi-q_v))+(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(e_li-q_l)))*_Vi[2+i];
						for(int j=0;j<_Ndim;j++)
							_AroePlusImplicit[(2+i)*_nVar+2+j]=rhomi*_Vi[2+i]*_vec_normal[j];
						_AroePlusImplicit[(2+i)*_nVar+2+i]+=rhomi*ui_n;
						_AroePlusImplicit[(2+i)*_nVar+2+_Ndim]=-rhomi*rhomi*ui_n*(cv_v*cmi/(rho_vi*(e_vi-q_v))+(1-cmi)/(rho_li*(e_li-q_l)))*_Vi[2+i];
					}
					_AroePlusImplicit[(2+_Ndim)*_nVar+0]=ui_n*(rhomi*(e_vi-e_li)-Emi*rhomi*rhomi*(1/rho_vi-1/rho_li));
					_AroePlusImplicit[(2+_Ndim)*_nVar+1]=rhomi*rhomi*Emi*(cmi/(rho_vi*rho_vi*(gamma_v-1)*(e_vi-q_v))+(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(e_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[(2+_Ndim)*_nVar+2+i]=(rhomi*Emi+pj)*_vec_normal[i] + ui_n*rhomi*_Vi[2+i];
					_AroePlusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=ui_n*(rhomi*(cv_v*cmi+cv_l*(1-cmi))-Emi*rhomi*rhomi*(cv_v*cmi/(rho_vi*(e_vi-q_v))+cv_l*(1-cmi)/(rho_li*(e_li-q_l))));


					/******** Construction de la matrice A^- *********/
					_AroeMinusImplicit[0*_nVar+0]=0;
					_AroeMinusImplicit[0*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[0*_nVar+2+i]=0;
					_AroeMinusImplicit[0*_nVar+2+_Ndim]=0;
					_AroeMinusImplicit[1*_nVar+0]=0;
					_AroeMinusImplicit[1*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[1*_nVar+2+i]=0;
					_AroeMinusImplicit[1*_nVar+2+_Ndim]=0;
					for(int i=0;i<_Ndim;i++)
					{
						_AroeMinusImplicit[(2+i)*_nVar+0]=0;
						_AroeMinusImplicit[(2+i)*_nVar+1]=_vec_normal[i];
						for(int j=0;j<_Ndim;j++)
							_AroeMinusImplicit[(2+i)*_nVar+2+j]=0;
						_AroeMinusImplicit[(2+i)*_nVar+2+i]+=0;
						_AroeMinusImplicit[(2+i)*_nVar+2+_Ndim]=0;
					}
					_AroeMinusImplicit[(2+_Ndim)*_nVar+0]=0;
					_AroeMinusImplicit[(2+_Ndim)*_nVar+1]=ui_n;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[(2+_Ndim)*_nVar+2+i]=0;
					_AroeMinusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=0;
				}
				else if(u_mn<-_precision)
				{
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"VFFC Staggered state rhomj="<<rhomj<<" cmj= "<<cmj<<" Emj= "<<Emj<<" Tmj= "<<Tmj<<" pi= "<<pi<<endl;

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					getMixturePressureDerivatives( mi_v, mi_l, pi, Tmi);
					if(_kappa*hmi+_khi+cmi*_ksi<0)
					{
						*_runLogFile<<"staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe"<<endl;
						throw CdmathException("staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe");
					}
					double ami=sqrt(_kappa*hmi+_khi+cmi*_ksi);//vitesse du son du melange
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", ami= "<<ami<<endl;

					//On remplit les valeurs propres
					vp_dist[0]=uj_n+ami;
					vp_dist[1]=uj_n-ami;
					vp_dist[2]=uj_n;

					_maxvploc=fabs(uj_n)+ami;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					/******** Construction de la matrice A^+ *********/
					_AroePlusImplicit[0*_nVar+0]=0;
					_AroePlusImplicit[0*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[0*_nVar+2+i]=0;
					_AroePlusImplicit[0*_nVar+2+_Ndim]=0;
					_AroePlusImplicit[1*_nVar+0]=0;
					_AroePlusImplicit[1*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[1*_nVar+2+i]=0;
					_AroePlusImplicit[1*_nVar+2+_Ndim]=0;
					for(int i=0;i<_Ndim;i++)
					{
						_AroePlusImplicit[(2+i)*_nVar+0]=0;
						_AroePlusImplicit[(2+i)*_nVar+1]=_vec_normal[i];
						for(int j=0;j<_Ndim;j++)
							_AroePlusImplicit[(2+i)*_nVar+2+j]=0;
						_AroePlusImplicit[(2+i)*_nVar+2+i]+=0;
						_AroePlusImplicit[(2+i)*_nVar+2+_Ndim]=0;
					}
					_AroePlusImplicit[(2+_Ndim)*_nVar+0]=0;
					_AroePlusImplicit[(2+_Ndim)*_nVar+1]=uj_n;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[(2+_Ndim)*_nVar+2+i]=0;
					_AroePlusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=0;

					/******** Construction de la matrice A^- *********/
					_AroeMinusImplicit[0*_nVar+0]=-rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj);
					_AroeMinusImplicit[0*_nVar+1]= rhomj*rhomj*uj_n*(cmj/(rho_vj*rho_vj*(gamma_v-1)*(e_vj-q_v))+(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(e_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[0*_nVar+2+i]=rhomj*_vec_normal[i];
					_AroeMinusImplicit[0*_nVar+2+_Ndim]=-rhomj*rhomj*uj_n*(cv_v*cmj/(rho_vj*(e_vj-q_v))+cv_l*(1-cmj)/(rho_lj*(e_lj-q_l)));
					_AroeMinusImplicit[1*_nVar+0]=rhomj*uj_n-cmj*rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj);
					_AroeMinusImplicit[1*_nVar+1]=-cmj*rhomj*rhomj*uj_n*(cmj/(rho_vj*rho_vj*(gamma_v-1)*(e_vj-q_v))+(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(e_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[1*_nVar+2+i]=cmj*rhomj*_vec_normal[i];
					_AroeMinusImplicit[1*_nVar+2+_Ndim]=-cmj*rhomj*rhomj*uj_n*(cv_v*cmj/(rho_vj*(e_vj-q_v))+cv_l*(1-cmj)/(rho_lj*(e_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
					{
						_AroeMinusImplicit[(2+i)*_nVar+0]=-rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj)*_Vj[2+i];
						_AroeMinusImplicit[(2+i)*_nVar+1]=rhomj*rhomj*uj_n*(cmj/(rho_vj*rho_vj*(gamma_v-1)*(e_vj-q_v))+(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(e_lj-q_l)))*_Vj[2+i];
						for(int j=0;j<_Ndim;j++)
							_AroeMinusImplicit[(2+i)*_nVar+2+j]=rhomj*_Vj[2+i]*_vec_normal[j];
						_AroeMinusImplicit[(2+i)*_nVar+2+i]+=rhomj*uj_n;
						_AroeMinusImplicit[(2+i)*_nVar+2+_Ndim]=-rhomj*rhomj*uj_n*(cv_v*cmj/(rho_vj*(e_vj-q_v))+cv_l*(1-cmj)/(rho_lj*(e_lj-q_l)))*_Vj[2+i];
					}
					_AroeMinusImplicit[(2+_Ndim)*_nVar+0]=uj_n*(rhomj*(e_vj-e_lj)-Emj*rhomj*rhomj*(1/rho_vj-1/rho_lj));
					_AroeMinusImplicit[(2+_Ndim)*_nVar+1]=rhomj*rhomj*Emj*(cmj/(rho_vj*rho_vj*(gamma_v-1)*(e_vj-q_v))+(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(e_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[(2+_Ndim)*_nVar+2+i]=(rhomj*Emj+pi)*_vec_normal[i] + uj_n*rhomj*_Vj[2+i];
					_AroeMinusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=uj_n*(rhomj*(cv_v*cmj+cv_l*(1-cmj))-Emj*rhomj*rhomj*(cv_v*cmj/(rho_vj*(e_vj-q_v))+cv_l*(1-cmj)/(rho_lj*(e_lj-q_l))));
				}
				else
				{
					*_runLogFile<< "DriftModel::staggeredVFFCMatricesPrimitiveVariables: velocity umn should be non zero" << endl;
					_runLogFile->close();
					throw CdmathException("DriftModel::staggeredVFFCMatricesPrimitiveVariables: velocity umn should be non zero");
				}
			}
			else if(_useDellacherieEOS)
			{
				StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
				StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);
				double cp_v=fluide0->constante("cp");
				double cp_l=fluide1->constante("cp");

				double h_vi=_fluides[0]->getEnthalpy(Tmi,rho_vi);
				double h_li=_fluides[1]->getEnthalpy(Tmi,rho_li);
				double h_vj=_fluides[0]->getEnthalpy(Tmj,rho_vj);
				double h_lj=_fluides[1]->getEnthalpy(Tmj,rho_lj);
				double Hmi=Emi+pi/rho_vi;
				double Hmj=Emj+pj/rho_vj;

				if(u_mn>_precision)
				{
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"VFFC Staggered state rhomi="<<rhomi<<" cmi= "<<cmi<<" Hmi= "<<Hmi<<" Tmi= "<<Tmi<<" pj= "<<pj<<endl;

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					getMixturePressureDerivatives( mj_v, mj_l, pj, Tmj);
					if(_kappa*hmj+_khi+cmj*_ksi<0)
					{
						*_runLogFile<<"staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe"<<endl;
						throw CdmathException("staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe");
					}
					double amj=sqrt(_kappa*hmj+_khi+cmj*_ksi);//vitesse du son du melange

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", amj= "<<amj<<endl;

					//On remplit les valeurs propres
					vp_dist[0]=ui_n+amj;
					vp_dist[1]=ui_n-amj;
					vp_dist[2]=ui_n;

					_maxvploc=fabs(ui_n)+amj;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					/******** Construction de la matrice A^+ *********/
					_AroePlusImplicit[0*_nVar+0]=-rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li);
					_AroePlusImplicit[0*_nVar+1]= rhomi*rhomi*ui_n*(gamma_v*cmi/(rho_vi*rho_vi*(gamma_v-1)*(h_vi-q_v))+gamma_l*(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(h_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[0*_nVar+2+i]=rhomi*_vec_normal[i];
					_AroePlusImplicit[0*_nVar+2+_Ndim]=-rhomi*rhomi*ui_n*(cp_v*cmi/(rho_vi*(h_vi-q_v))+cp_l*(1-cmi)/(rho_li*(h_li-q_l)));
					_AroePlusImplicit[1*_nVar+0]=rhomi*ui_n-cmi*rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li);
					_AroePlusImplicit[1*_nVar+1]=-cmi*rhomi*rhomi*ui_n*(gamma_v*cmi/(rho_vi*rho_vi*(gamma_v-1)*(h_vi-q_v))+gamma_l*(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(h_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[1*_nVar+2+i]=cmi*rhomi*_vec_normal[i];
					_AroePlusImplicit[1*_nVar+2+_Ndim]=-cmi*rhomi*rhomi*ui_n*(cp_v*cmi/(rho_vi*(h_vi-q_v))+cp_l*(1-cmi)/(rho_li*(h_li-q_l)));
					for(int i=0;i<_Ndim;i++)
					{
						_AroePlusImplicit[(2+i)*_nVar+0]=-rhomi*rhomi*ui_n*(1/rho_vi-1/rho_li)*_Vi[2+i];
						_AroePlusImplicit[(2+i)*_nVar+1]=rhomi*rhomi*ui_n*(gamma_v*cmi/(rho_vi*rho_vi*(gamma_v-1)*(h_vi-q_v))+gamma_l*(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(h_li-q_l)))*_Vi[2+i];
						for(int j=0;j<_Ndim;j++)
							_AroePlusImplicit[(2+i)*_nVar+2+j]=rhomi*_Vi[2+i]*_vec_normal[j];
						_AroePlusImplicit[(2+i)*_nVar+2+i]+=rhomi*ui_n;
						_AroePlusImplicit[(2+i)*_nVar+2+_Ndim]=-rhomi*rhomi*ui_n*(cp_v*cmi/(rho_vi*(h_vi-q_v))+cp_l*(1-cmi)/(rho_li*(h_li-q_l)))*_Vi[2+i];
					}
					_AroePlusImplicit[(2+_Ndim)*_nVar+0]=ui_n*(rhomi*(h_vi-h_li)-Hmi*rhomi*rhomi*(1/rho_vi-1/rho_li));
					_AroePlusImplicit[(2+_Ndim)*_nVar+1]=rhomi*rhomi*Hmi*(gamma_v*cmi/(rho_vi*rho_vi*(gamma_v-1)*(h_vi-q_v))+gamma_l*(1-cmi)/(rho_li*rho_li*(gamma_l-1)*(h_li-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[(2+_Ndim)*_nVar+2+i]=(rhomi*Emi+pj)*_vec_normal[i] + ui_n*rhomi*_Vi[2+i];
					_AroePlusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=ui_n*(rhomi*(cp_v*cmi+cp_l*(1-cmi))-Hmi*rhomi*rhomi*(cp_v*cmi/(rho_vi*(h_vi-q_v))+cp_l*(1-cmi)/(rho_li*(h_li-q_l))));


					/******** Construction de la matrice A^- *********/
					_AroeMinusImplicit[0*_nVar+0]=0;
					_AroeMinusImplicit[0*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[0*_nVar+2+i]=0;
					_AroeMinusImplicit[0*_nVar+2+_Ndim]=0;
					_AroeMinusImplicit[1*_nVar+0]=0;
					_AroeMinusImplicit[1*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[1*_nVar+2+i]=0;
					_AroeMinusImplicit[1*_nVar+2+_Ndim]=0;
					for(int i=0;i<_Ndim;i++)
					{
						_AroeMinusImplicit[(2+i)*_nVar+0]=0;
						_AroeMinusImplicit[(2+i)*_nVar+1]=_vec_normal[i];
						for(int j=0;j<_Ndim;j++)
							_AroeMinusImplicit[(2+i)*_nVar+2+j]=0;
						_AroeMinusImplicit[(2+i)*_nVar+2+i]+=0;
						_AroeMinusImplicit[(2+i)*_nVar+2+_Ndim]=0;
					}
					_AroeMinusImplicit[(2+_Ndim)*_nVar+0]=0;
					_AroeMinusImplicit[(2+_Ndim)*_nVar+1]=ui_n;
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[(2+_Ndim)*_nVar+2+i]=0;
					_AroeMinusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=0;
				}
				else if(u_mn<-_precision)
				{
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"VFFC Staggered state rhomj="<<rhomj<<" cmj= "<<cmj<<" Hmj= "<<Hmj<<" Tmj= "<<Tmj<<" pi= "<<pi<<endl;

					/***********Calcul des valeurs propres ********/
					vector<complex<double> > vp_dist(3,0);
					getMixturePressureDerivatives( mi_v, mi_l, pi, Tmi);
					if(_kappa*hmi+_khi+cmi*_ksi<0)
					{
						*_runLogFile<<"staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe"<<endl;
						throw CdmathException("staggeredVFFCMatricesPrimitiveVariables: vitesse du son complexe");
					}
					double ami=sqrt(_kappa*hmi+_khi+cmi*_ksi);//vitesse du son du melange
					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"_khi= "<<_khi<<", _kappa= "<< _kappa << ", _ksi= "<<_ksi <<", ami= "<<ami<<endl;

					//On remplit les valeurs propres
					vp_dist[0]=uj_n+ami;
					vp_dist[1]=uj_n-ami;
					vp_dist[2]=uj_n;

					_maxvploc=fabs(uj_n)+ami;
					if(_maxvploc>_maxvp)
						_maxvp=_maxvploc;

					/******** Construction de la matrice A^+ *********/
					_AroePlusImplicit[0*_nVar+0]=0;
					_AroePlusImplicit[0*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[0*_nVar+2+i]=0;
					_AroePlusImplicit[0*_nVar+2+_Ndim]=0;
					_AroePlusImplicit[1*_nVar+0]=0;
					_AroePlusImplicit[1*_nVar+1]=0;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[1*_nVar+2+i]=0;
					_AroePlusImplicit[1*_nVar+2+_Ndim]=0;
					for(int i=0;i<_Ndim;i++)
					{
						_AroePlusImplicit[(2+i)*_nVar+0]=0;
						_AroePlusImplicit[(2+i)*_nVar+1]=_vec_normal[i];
						for(int j=0;j<_Ndim;j++)
							_AroePlusImplicit[(2+i)*_nVar+2+j]=0;
						_AroePlusImplicit[(2+i)*_nVar+2+i]+=0;
						_AroePlusImplicit[(2+i)*_nVar+2+_Ndim]=0;
					}
					_AroePlusImplicit[(2+_Ndim)*_nVar+0]=0;
					_AroePlusImplicit[(2+_Ndim)*_nVar+1]=uj_n;
					for(int i=0;i<_Ndim;i++)
						_AroePlusImplicit[(2+_Ndim)*_nVar+2+i]=0;
					_AroePlusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=0;

					/******** Construction de la matrice A^- *********/
					_AroeMinusImplicit[0*_nVar+0]=-rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj);
					_AroeMinusImplicit[0*_nVar+1]= rhomj*rhomj*uj_n*(gamma_v*cmj/(rho_vj*rho_vj*(gamma_v-1)*(h_vj-q_v))+gamma_l*(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(h_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[0*_nVar+2+i]=rhomj*_vec_normal[i];
					_AroeMinusImplicit[0*_nVar+2+_Ndim]=-rhomj*rhomj*uj_n*(cp_v*cmj/(rho_vj*(h_vj-q_v))+cp_l*(1-cmj)/(rho_lj*(h_lj-q_l)));
					_AroeMinusImplicit[1*_nVar+0]=rhomj*uj_n-cmj*rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj);
					_AroeMinusImplicit[1*_nVar+1]=-cmj*rhomj*rhomj*uj_n*(gamma_v*cmj/(rho_vj*rho_vj*(gamma_v-1)*(h_vj-q_v))+gamma_l*(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(h_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[1*_nVar+2+i]=cmj*rhomj*_vec_normal[i];
					_AroeMinusImplicit[1*_nVar+2+_Ndim]=-cmj*rhomj*rhomj*uj_n*(cp_v*cmj/(rho_vj*(h_vj-q_v))+cp_l*(1-cmj)/(rho_lj*(h_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
					{
						_AroeMinusImplicit[(2+i)*_nVar+0]=-rhomj*rhomj*uj_n*(1/rho_vj-1/rho_lj)*_Vj[2+i];
						_AroeMinusImplicit[(2+i)*_nVar+1]=rhomj*rhomj*uj_n*(gamma_v*cmj/(rho_vj*rho_vj*(gamma_v-1)*(h_vj-q_v))+gamma_l*(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(h_lj-q_l)))*_Vj[2+i];
						for(int j=0;j<_Ndim;j++)
							_AroeMinusImplicit[(2+i)*_nVar+2+j]=rhomj*_Vj[2+i]*_vec_normal[j];
						_AroeMinusImplicit[(2+i)*_nVar+2+i]+=rhomj*uj_n;
						_AroeMinusImplicit[(2+i)*_nVar+2+_Ndim]=-rhomj*rhomj*uj_n*(cp_v*cmj/(rho_vj*(h_vj-q_v))+cp_l*(1-cmj)/(rho_lj*(h_lj-q_l)))*_Vj[2+i];
					}
					_AroeMinusImplicit[(2+_Ndim)*_nVar+0]=uj_n*(rhomj*(h_vj-h_lj)-Hmj*rhomj*rhomj*(1/rho_vj-1/rho_lj));
					_AroeMinusImplicit[(2+_Ndim)*_nVar+1]=rhomj*rhomj*Hmj*(gamma_v*cmj/(rho_vj*rho_vj*(gamma_v-1)*(h_vj-q_v))+gamma_l*(1-cmj)/(rho_lj*rho_lj*(gamma_l-1)*(h_lj-q_l)));
					for(int i=0;i<_Ndim;i++)
						_AroeMinusImplicit[(2+_Ndim)*_nVar+2+i]=(rhomj*Emj+pi)*_vec_normal[i] + uj_n*rhomj*_Vj[2+i];
					_AroeMinusImplicit[(2+_Ndim)*_nVar+2+_Ndim]=uj_n*(rhomj*(cp_v*cmj+cp_l*(1-cmj))-Hmj*rhomj*rhomj*(cp_v*cmj/(rho_vj*(h_vj-q_v))+cp_l*(1-cmj)/(rho_lj*(h_lj-q_l))));
				}
				else
				{
					*_runLogFile<< "DriftModel::staggeredVFFCMatricesPrimitiveVariables: velocity umn should be non zero" << endl;
					_runLogFile->close();
					throw CdmathException("DriftModel::staggeredVFFCMatricesPrimitiveVariables: velocity umn should be non zero");
				}
			}
			else
				throw CdmathException("DriftModel::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie");
		}
		else//case nil velocity on the interface, multiply by jacobian matrix
		{
			Polynoms Poly;
			primToConsJacobianMatrix(_Vj);
			Poly.matrixProduct(_AroeMinus, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroeMinusImplicit);
			primToConsJacobianMatrix(_Vi);
			Poly.matrixProduct(_AroePlus,  _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroePlusImplicit);
		}
	}
}

void DriftModel::applyVFRoeLowMachCorrections(bool isBord, string nameOfGroup)
{
	if(_nonLinearFormulation!=VFRoe)
		throw CdmathException("DriftModel::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation");
	else//_nonLinearFormulation==VFRoe
	{
		if(_spaceScheme==lowMach){
			double u_2=0;
			for(int i=0;i<_Ndim;i++)
				u_2 += _Uroe[2+i]*_Uroe[2+i];

			double 	c = _maxvploc;//mixture sound speed
			double M=sqrt(u_2)/c;//Mach number
			_Vij[1]=M*_Vij[1]+(1-M)*(_Vi[1]+_Vj[1])/2;
		}
		else if(_spaceScheme==pressureCorrection)
		{
			/* 
			 * option 1 : no pressure correction except for inner walls where p^*=p_int
			 * option 2 : Clerc pressure correction everywhere, except for inner walls where p^*=p_int
			 * option 3 : Clerc pressure correction everywhere, even for inner walls
			 * option 4 : Clerc pressure correction only inside, upwind on wall and innerwall boundaries
			 * option 5 : Clerc pressure correction inside the domain and special pressure correction involving gravity at wall and inner wall boundaries 
			 * */
			bool isWall=false, isInnerWall=false;
			if(isBord)
			{
				isWall = (_limitField[nameOfGroup].bcType==Wall);
				isInnerWall = (_limitField[nameOfGroup].bcType==InnerWall);
			}

			if(	(!isInnerWall && _pressureCorrectionOrder==2)  || _pressureCorrectionOrder==3 ||
				(!isBord &&	_pressureCorrectionOrder==4) || (!isWall && !isInnerWall && _pressureCorrectionOrder==5)   )//Clerc pressure correction
			{
				double norm_uij=0, uij_n=0, ui_n=0, uj_n=0;
				for(int i=0;i<_Ndim;i++)
				{
					norm_uij += _Uroe[2+i]*_Uroe[2+i];
					uij_n += _Uroe[2+i]*_vec_normal[i];
					ui_n += _Vi[2+i]*_vec_normal[i];
					uj_n += _Vj[2+i]*_vec_normal[i];
				}
				norm_uij=sqrt(norm_uij);
					if(norm_uij>_precision)//avoid division by zero
					_Vij[1]=(_Vi[1]+_Vj[1])/2 + uij_n/norm_uij*(_Vj[1]-_Vi[1])/4 - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
				else
					_Vij[1]=(_Vi[1]+_Vj[1])/2                                    - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
			}
			else if((isWall || isInnerWall) && _pressureCorrectionOrder==5)//correction de pression gravitaire
			{
				double g_n=0;//scalar product of gravity and normal vector
				for(int i=0;i<_Ndim;i++)
					g_n += _GravityField3d[i]*_vec_normal[i];
				_Vij[1]=_Vi[1]+ _Ui[0]*g_n/_inv_dxi/2;
			}
			else if(isInnerWall && (_pressureCorrectionOrder==1 || _pressureCorrectionOrder==2) )
					_Vij[1]=_Vi[1];
		}
		else if(_spaceScheme==staggered)
		{
			double uij_n=0;
			for(int i=0;i<_Ndim;i++)
				uij_n += _Uroe[1+i]*_vec_normal[i];

			if(uij_n>_precision){
				_Vij[0]=_Vi[0];
				_Vij[1]=_Vj[1];
				for(int i=0;i<_Ndim;i++)
					_Vij[2+i]=_Vi[2+i];
				_Vij[_nVar-1]=_Vi[_nVar-1];
			}
			else if(uij_n<-_precision){
				_Vij[0]=_Vj[0];
				_Vij[1]=_Vi[1];
				for(int i=0;i<_Ndim;i++)
					_Vij[2+i]=_Vj[2+i];
				_Vij[_nVar-1]=_Vj[_nVar-1];
			}
			else{
				_Vij[0]=(_Vi[0]+_Vj[0])/2;
				_Vij[1]=(_Vi[1]+_Vj[1])/2;
				for(int i=0;i<_Ndim;i++)
					_Vij[2+i]=(_Vi[2+i]+_Vj[2+i])/2;
				_Vij[_nVar-1]=(_Vi[_nVar-1]+_Vj[_nVar-1])/2;
			}
		}
		primToCons(_Vij,0,_Uij,0);
	}
}
void DriftModel::RoeMatrixConservativeVariables(double cm, double umn,double ecin, double Hm,Vector vitesse)
{
	//On remplit la matrice de Roe en variables conservatives : F(u_L)-F(U_R)=Aroe (U_L-U_R)
	_Aroe[0*_nVar+0]=0;
	_Aroe[0*_nVar+1]=0;
	for(int i=0;i<_Ndim;i++)
		_Aroe[0*_nVar+2+i]=_vec_normal[i];
	_Aroe[0*_nVar+2+_Ndim]=0;
	_Aroe[1*_nVar+0]=-umn*cm;
	_Aroe[1*_nVar+1]=umn;
	for(int i=0;i<_Ndim;i++)
		_Aroe[1*_nVar+2+i]=cm*_vec_normal[i];
	_Aroe[1*_nVar+2+_Ndim]=0;
	for(int i=0;i<_Ndim;i++)
	{
		_Aroe[(2+i)*_nVar+0]=(_khi+_kappa*ecin)*_vec_normal[i]-umn*vitesse[i];
		_Aroe[(2+i)*_nVar+1]=_ksi*_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_Aroe[(2+i)*_nVar+2+j]=vitesse[i]*_vec_normal[j]-_kappa*_vec_normal[i]*vitesse[j];
		_Aroe[(2+i)*_nVar+2+i]+=umn;
		_Aroe[(2+i)*_nVar+2+_Ndim]=_kappa*_vec_normal[i];
	}
	_Aroe[(2+_Ndim)*_nVar+0]=(_khi+_kappa*ecin-Hm)*umn;
	_Aroe[(2+_Ndim)*_nVar+1]=_ksi*umn;
	for(int i=0;i<_Ndim;i++)
		_Aroe[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i]-_kappa*umn*vitesse[i];
	_Aroe[(2+_Ndim)*_nVar+2+_Ndim]=(_kappa+1)*umn;
}

void DriftModel::staggeredRoeUpwindingMatrixConservativeVariables(double cm, double umn,double ecin, double Hm, Vector vitesse)
{
	//Calcul de décentrement de type décalé pour formulation Roe
	if(abs(umn)>_precision)//non zero velocity at the interface
	{
		_absAroe[0*_nVar+0]=0;
		_absAroe[0*_nVar+1]=0;
		for(int i=0;i<_Ndim;i++)
			_absAroe[0*_nVar+2+i]=_vec_normal[i];
		_absAroe[0*_nVar+2+_Ndim]=0;
		_absAroe[1*_nVar+0]=-umn*cm;
		_absAroe[1*_nVar+1]=umn;
		for(int i=0;i<_Ndim;i++)
			_absAroe[1*_nVar+2+i]=cm*_vec_normal[i];
		_absAroe[1*_nVar+2+_Ndim]=0;
		for(int i=0;i<_Ndim;i++)
		{
			_absAroe[(2+i)*_nVar+0]=-(_khi+_kappa*ecin)*_vec_normal[i]-umn*vitesse[i];
			_absAroe[(2+i)*_nVar+1]=-_ksi*_vec_normal[i];
			for(int j=0;j<_Ndim;j++)
				_absAroe[(2+i)*_nVar+2+j]=vitesse[i]*_vec_normal[j]+_kappa*_vec_normal[i]*vitesse[j];
			_absAroe[(2+i)*_nVar+2+i]+=umn;
			_absAroe[(2+i)*_nVar+2+_Ndim]=-_kappa*_vec_normal[i];
		}
		_absAroe[(2+_Ndim)*_nVar+0]=(-_khi-_kappa*ecin-Hm)*umn;
		_absAroe[(2+_Ndim)*_nVar+1]=-_ksi*umn;
		for(int i=0;i<_Ndim;i++)
			_absAroe[(2+_Ndim)*_nVar+2+i]=Hm*_vec_normal[i]+_kappa*umn*vitesse[i];
		_absAroe[(2+_Ndim)*_nVar+2+_Ndim]=(-_kappa+1)*umn;

		double signu;
		if(umn>_precision)
			signu=1;
		else // umn<-_precision
			signu=-1;

		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] *= signu;
	}
	else//umn=0>centered scheme
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] =0;
		//for(int i=0; i<_nVar;i++)
		//	_absAroe[i+_nVar*i] =_maxvploc;
	}
}

void DriftModel::staggeredRoeUpwindingMatrixPrimitiveVariables(double concentration, double rhom, double umn, double Hm, Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//Calcul de décentrement de type décalé pour formulation Roe
	_absAroeImplicit[0*_nVar+0]= _drho_sur_dcv*umn;
	_absAroeImplicit[0*_nVar+1]= _drho_sur_dp *umn;
	for(int i=0;i<_Ndim;i++)
		_absAroeImplicit[0*_nVar+2+i]=rhom*_vec_normal[i];
	_absAroeImplicit[0*_nVar+2+_Ndim]=_drho_sur_dT*umn;
	_absAroeImplicit[1*_nVar+0]= _drhocv_sur_dcv*umn;
	_absAroeImplicit[1*_nVar+1]= _drhocv_sur_dp *umn;
	for(int i=0;i<_Ndim;i++)
		_absAroeImplicit[1*_nVar+2+i]=rhom*concentration*_vec_normal[i];
	_absAroeImplicit[1*_nVar+2+_Ndim]=_drhocv_sur_dT*umn;
	for(int i=0;i<_Ndim;i++)
	{
		_absAroeImplicit[(2+i)*_nVar+0]= _drho_sur_dcv*umn*vitesse[i];
		_absAroeImplicit[(2+i)*_nVar+1]= _drho_sur_dp *umn*vitesse[i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_absAroeImplicit[(2+i)*_nVar+2+j]=rhom*vitesse[i]*_vec_normal[j];
		_absAroeImplicit[(2+i)*_nVar+2+i]+=rhom*umn;
		_absAroeImplicit[(2+i)*_nVar+2+_Ndim]=_drho_sur_dT*umn*vitesse[i];
	}
	_absAroeImplicit[(2+_Ndim)*_nVar+0]=  _drhoE_sur_dcv  *umn;
	_absAroeImplicit[(2+_Ndim)*_nVar+1]=( _drhoE_sur_dp+1)*umn;
	for(int i=0;i<_Ndim;i++)
		_absAroeImplicit[(2+_Ndim)*_nVar+2+i]=rhom*(Hm*_vec_normal[i]+umn*vitesse[i]);
	_absAroeImplicit[(2+_Ndim)*_nVar+2+_Ndim]=_drhoE_sur_dT*umn;
}

void DriftModel::convectionMatrixPrimitiveVariables(double concentration, double rhom, double umn, double Hm, Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//On remplit la matrice de Roe en variables primitives : F(V_L)-F(V_R)=Aroe (V_L-V_R)
	//EOS is more involved with primitive variables
	//first call to getDensityDerivatives(double concentration, double pression, double temperature,double v2) needed
	_AroeImplicit[0*_nVar+0]=_drho_sur_dcv*umn;
	_AroeImplicit[0*_nVar+1]=_drho_sur_dp*umn;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+2+i]=rhom*_vec_normal[i];
	_AroeImplicit[0*_nVar+2+_Ndim]=_drho_sur_dT*umn;
	_AroeImplicit[1*_nVar+0]=_drhocv_sur_dcv*umn;
	_AroeImplicit[1*_nVar+1]=_drhocv_sur_dp*umn;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[1*_nVar+2+i]=rhom*concentration*_vec_normal[i];
	_AroeImplicit[1*_nVar+2+_Ndim]=_drhocv_sur_dT*umn;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(2+i)*_nVar+0]=_drho_sur_dcv*umn*vitesse[i];
		_AroeImplicit[(2+i)*_nVar+1]=_drho_sur_dp *umn*vitesse[i]+_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(2+i)*_nVar+2+j]=rhom*vitesse[i]*_vec_normal[j];
		_AroeImplicit[(2+i)*_nVar+2+i]+=rhom*umn;
		_AroeImplicit[(2+i)*_nVar+2+_Ndim]=_drho_sur_dT*umn*vitesse[i];
	}
	_AroeImplicit[(2+_Ndim)*_nVar+0]= _drhoE_sur_dcv  *umn;
	_AroeImplicit[(2+_Ndim)*_nVar+1]=(_drhoE_sur_dp+1)*umn;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(2+_Ndim)*_nVar+2+i]=rhom*(Hm*_vec_normal[i]+umn*vitesse[i]);
	_AroeImplicit[(2+_Ndim)*_nVar+2+_Ndim]=_drhoE_sur_dT*umn;
}

void DriftModel::getDensityDerivatives(double concentration, double pression, double temperature,double v2)
{
	//EOS is more involved with primitive variables

	double rho_v=_fluides[0]->getDensity(pression,temperature);
	double rho_l=_fluides[1]->getDensity(pression,temperature);
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double q_v=_fluides[0]->constante("q");
	double q_l=_fluides[1]->constante("q");

	double rho=concentration*rho_v+(1-concentration)*rho_l;;

	if(	!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
		StiffenedGas* fluide1=dynamic_cast<StiffenedGas*>(_fluides[1]);
		double e_v = fluide0->getInternalEnergy(temperature);
		double e_l = fluide1->getInternalEnergy(temperature);
		double cv_v=fluide0->constante("cv");
		double cv_l=fluide1->constante("cv");
		double e=concentration*e_v+(1-concentration)*e_l;
		double E=e+0.5*v2;

		_drho_sur_dcv=-rho*rho*(1/rho_v-1/rho_l);
		_drho_sur_dp=rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		_drho_sur_dT=-rho*rho*(
				cv_v*   concentration /(rho_v*(e_v-q_v))
				+cv_l*(1-concentration)/(rho_l*(e_l-q_l))
		);

		_drhocv_sur_dcv=rho-concentration*rho*rho*(1/rho_v-1/rho_l);
		_drhocv_sur_dp=concentration*rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		_drhocv_sur_dT=-concentration*rho*rho*(
				cv_v*   concentration /(rho_v*(e_v-q_v))
				+cv_l*(1-concentration)/(rho_l*(e_l-q_l))
		);
		_drhoE_sur_dcv=rho*(e_v-e_l)-E*rho*rho*(1/rho_v-1/rho_l);
		_drhoE_sur_dp=E*rho*rho*(
				concentration /(rho_v*rho_v*(gamma_v-1)*(e_v-q_v))
				+(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(e_l-q_l))
		);
		_drhoE_sur_dT=rho*(cv_v*concentration + cv_l*(1-concentration))
																																																																																			-rho*rho*E*( cv_v*   concentration /(rho_v*(e_v-q_v))
																																																																																					+cv_l*(1-concentration)/(rho_l*(e_l-q_l)));
	}
	else if(_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_fluides[0]);
		StiffenedGasDellacherie* fluide1=dynamic_cast<StiffenedGasDellacherie*>(_fluides[1]);
		double h_v=fluide0->getEnthalpy(temperature);
		double h_l=fluide1->getEnthalpy(temperature);
		double h=concentration*h_v+(1-concentration)*h_l;
		double H=h+0.5*v2;
		double cp_v=fluide0->constante("cp");
		double cp_l=fluide1->constante("cp");

		_drho_sur_dcv=-rho*rho*(1/rho_v-1/rho_l);
		_drho_sur_dp =rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		);
		_drho_sur_dT=-rho*rho*(
				cp_v*   concentration /(rho_v*(h_v-q_v))
				+cp_l*(1-concentration)/(rho_l*(h_l-q_l))
		);

		_drhocv_sur_dcv=rho-concentration*rho*rho*(1/rho_v-1/rho_l);
		_drhocv_sur_dp=    concentration*rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		);
		_drhocv_sur_dT=-concentration*rho*rho*(
				cp_v*   concentration /(rho_v*(h_v-q_v))
				+cp_l*(1-concentration)/(rho_l*(h_l-q_l))
		);
		_drhoE_sur_dcv=rho*(h_v-h_l)-H*rho*rho*(1/rho_v-1/rho_l);
		_drhoE_sur_dp=H*rho*rho*(
				gamma_v*   concentration /(rho_v*rho_v*(gamma_v-1)*(h_v-q_v))
				+gamma_l*(1-concentration)/(rho_l*rho_l*(gamma_l-1)*(h_l-q_l))
		)-1;
		_drhoE_sur_dT=rho*(cp_v*concentration + cp_l*(1-concentration))
		           	    																																																																   -rho*rho*H*( cp_v*   concentration /(rho_v*(h_v-q_v))
		           	    																																																																		   +cp_l*(1-concentration)/(rho_l*(h_l-q_l)));
	}
	else
		throw CdmathException("DriftModel::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie");

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"_drho_sur_dcv= "<<_drho_sur_dcv<<", _drho_sur_dp= "<<_drho_sur_dp<<", _drho_sur_dT= "<<_drho_sur_dT<<endl;
		cout<<"_drhocv_sur_dcv= "<<_drhocv_sur_dcv<<", _drhocv_sur_dp= "<<_drhocv_sur_dp<<", _drhocv_sur_dT= "<<_drhocv_sur_dT<<endl;
		cout<<"_drhoE_sur_dcv= "<<_drhoE_sur_dcv<<", _drhoE_sur_dp= "<<_drhoE_sur_dp<<", _drhoE_sur_dT= "<<_drhoE_sur_dT<<endl;
	}
}

void DriftModel::save(){
	string prim(_path+"/DriftModelPrim_");
	string cons(_path+"/DriftModelCons_");
	string allFields(_path+"/");
	prim+=_fileName;
	cons+=_fileName;
	allFields+=_fileName;

	PetscInt Ii;
	for (long i = 0; i < _Nmailles; i++){
		// j = 0 : concentration, j=1 : pressure; j = _nVar - 1: temperature; j = 2,..,_nVar-2: velocity
		for (int j = 0; j < _nVar; j++){
			Ii = i*_nVar +j;
			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
		}
	}
	if(_saveConservativeField){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : total density; j = 1 : vapour density; j = _nVar - 1 : energy j = 2,..,_nVar-2: momentum
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
		if (_restartWithNewFileName)
			_restartWithNewFileName=false;
		string suppress_previous_runs ="rm -rf *"+_fileName+"_*";
		system(suppress_previous_runs.c_str());//Nettoyage des précédents calculs identiques

		if(_saveConservativeField){
			_UU.setInfoOnComponent(0,"Total_Density");// (kg/m^3)

			_UU.setInfoOnComponent(1,"Partial_Density");// (kg/m^3)
			_UU.setInfoOnComponent(2,"Momentum_x");// (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(3,"Momentum_y");// (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(4,"Momentum_z");// (kg/m^2/s)

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

		_VV.setInfoOnComponent(0,"Concentration");
		_VV.setInfoOnComponent(1,"Pressure_(Pa)");
		_VV.setInfoOnComponent(2,"Velocity_x_(m/s)");
		if (_Ndim>1)
			_VV.setInfoOnComponent(3,"Velocity_y_(m/s)");
		if (_Ndim>2)
			_VV.setInfoOnComponent(4,"Velocity_z_(m/s)");
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
	if(_saveVelocity){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : concentration, j=1 : pressure; j = _nVar - 1: temperature; j = 2,..,_nVar-2: velocity
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +2+j;
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
		double p,Tm,cv,alpha_v,rhom,rho_v,rho_l, m_v, m_l, h_v, h_l, vx,vy,vz;
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&rhom);
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&cv);
			Ii = i*_nVar +1;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&Tm);
			Ii = i*_nVar +2;
			VecGetValues(_primitiveVars,1,&Ii,&vx);
			if(_Ndim>1)
			{
				Ii = i*_nVar +3;
				VecGetValues(_primitiveVars,1,&Ii,&vy);
				if(_Ndim>2){
					Ii = i*_nVar +4;
					VecGetValues(_primitiveVars,1,&Ii,&vz);
				}
			}

			rho_v=_fluides[0]->getDensity(p,Tm);
			rho_l=_fluides[1]->getDensity(p,Tm);
			alpha_v=cv*rhom/rho_v;
			m_v=cv*rhom;
			m_l=(1-cv)*rhom;
			h_v=_fluides[0]->getEnthalpy(Tm,rho_v);
			h_l=_fluides[1]->getEnthalpy(Tm,rho_l);

			_VoidFraction(i)=alpha_v;
			_Enthalpy(i)=(m_v*h_v+m_l*h_l)/rhom;
			_Concentration(i)=cv;
			_mixtureDensity(i)=rhom;
			_Pressure(i)=p;
			_Temperature(i)=Tm;
			_DensiteLiquide(i)=rho_l;
			_DensiteVapeur(i)=rho_v;
			_EnthalpieLiquide(i)=h_l;
			_EnthalpieVapeur(i)=h_v;
			_VitesseX(i)=vx;
			if(_Ndim>1)
			{
				_VitesseY(i)=vy;
				if(_Ndim>2)
					_VitesseZ(i)=vz;
			}
		}
		_VoidFraction.setTime(_time,_nbTimeStep);
		_Enthalpy.setTime(_time,_nbTimeStep);
		_Concentration.setTime(_time,_nbTimeStep);
		_mixtureDensity.setTime(_time,_nbTimeStep);
		_Pressure.setTime(_time,_nbTimeStep);
		_Temperature.setTime(_time,_nbTimeStep);
		_DensiteLiquide.setTime(_time,_nbTimeStep);
		_DensiteVapeur.setTime(_time,_nbTimeStep);
		_EnthalpieLiquide.setTime(_time,_nbTimeStep);
		_EnthalpieVapeur.setTime(_time,_nbTimeStep);
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
				_VoidFraction.writeVTK(allFields+"_VoidFraction");
				_Enthalpy.writeVTK(allFields+"_Enthalpy");
				_Concentration.writeVTK(allFields+"_Concentration");
				_mixtureDensity.writeVTK(allFields+"_Density");
				_Pressure.writeVTK(allFields+"_Pressure");
				_Temperature.writeVTK(allFields+"_Temperature");
				_DensiteLiquide.writeVTK(allFields+"_LiquidDensity");
				_DensiteVapeur.writeVTK(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeVTK(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeVTK(allFields+"_SteamEnthalpy");
				_VitesseX.writeVTK(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ");
				}
				break;
			case MED :
				_VoidFraction.writeMED(allFields+"_VoidFraction");
				_Enthalpy.writeMED(allFields+"_Enthalpy");
				_Concentration.writeMED(allFields+"_Concentration");
				_mixtureDensity.writeMED(allFields+"_Density");
				_Pressure.writeMED(allFields+"_Pressure");
				_Temperature.writeMED(allFields+"_Temperature");
				_DensiteLiquide.writeMED(allFields+"_LiquidDensity");
				_DensiteVapeur.writeMED(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeMED(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeMED(allFields+"_SteamEnthalpy");
				_VitesseX.writeMED(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ");
				}
				break;
			case CSV :
				_VoidFraction.writeCSV(allFields+"_VoidFraction");
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Concentration.writeCSV(allFields+"_Concentration");
				_mixtureDensity.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
				_DensiteLiquide.writeCSV(allFields+"_LiquidDensity");
				_DensiteVapeur.writeCSV(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeCSV(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeCSV(allFields+"_SteamEnthalpy");
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
				_VoidFraction.writeVTK(allFields+"_VoidFraction",false);
				_Enthalpy.writeVTK(allFields+"_Enthalpy",false);
				_Concentration.writeVTK(allFields+"_Concentration",false);
				_mixtureDensity.writeVTK(allFields+"_Density",false);
				_Pressure.writeVTK(allFields+"_Pressure",false);
				_Temperature.writeVTK(allFields+"_Temperature",false);
				_DensiteLiquide.writeVTK(allFields+"_LiquidDensity",false);
				_DensiteVapeur.writeVTK(allFields+"_SteamDensity",false);
				_EnthalpieLiquide.writeVTK(allFields+"_LiquidEnthalpy",false);
				_EnthalpieVapeur.writeVTK(allFields+"_SteamEnthalpy",false);
				_VitesseX.writeVTK(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ",false);
				}
				break;
			case MED :
				_VoidFraction.writeMED(allFields+"_VoidFraction",false);
				_Enthalpy.writeMED(allFields+"_Enthalpy",false);
				_Concentration.writeMED(allFields+"_Concentration",false);
				_mixtureDensity.writeMED(allFields+"_Density",false);
				_Pressure.writeMED(allFields+"_Pressure",false);
				_Temperature.writeMED(allFields+"_Temperature",false);
				_DensiteLiquide.writeMED(allFields+"_LiquidDensity",false);
				_DensiteVapeur.writeMED(allFields+"_SteamDensity",false);
				_EnthalpieLiquide.writeMED(allFields+"_LiquidEnthalpy",false);
				_EnthalpieVapeur.writeMED(allFields+"_SteamEnthalpy",false);
				_VitesseX.writeMED(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ",false);
				}
				break;
			case CSV :
				_VoidFraction.writeCSV(allFields+"_VoidFraction");
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Concentration.writeCSV(allFields+"_Concentration");
				_mixtureDensity.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
				_DensiteLiquide.writeCSV(allFields+"_LiquidDensity");
				_DensiteVapeur.writeCSV(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeCSV(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeCSV(allFields+"_SteamEnthalpy");
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

		if(_saveVelocity){
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

		if(_saveAllFields)
		{
			allFields+="_Stat";
			switch(_saveFormat)
			{
			case VTK :
				_VoidFraction.writeVTK(allFields+"_VoidFraction");
				_Enthalpy.writeVTK(allFields+"_Enthalpy");
				_Concentration.writeVTK(allFields+"_Concentration");
				_mixtureDensity.writeVTK(allFields+"_Density");
				_Pressure.writeVTK(allFields+"_Pressure");
				_Temperature.writeVTK(allFields+"_Temperature");
				_DensiteLiquide.writeVTK(allFields+"_LiquidDensity");
				_DensiteVapeur.writeVTK(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeVTK(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeVTK(allFields+"_SteamEnthalpy");
				_VitesseX.writeVTK(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ");
				}
				break;
			case MED :
				_VoidFraction.writeMED(allFields+"_VoidFraction");
				_Enthalpy.writeMED(allFields+"_Enthalpy");
				_Concentration.writeMED(allFields+"_Concentration");
				_mixtureDensity.writeMED(allFields+"_Density");
				_Pressure.writeMED(allFields+"_Pressure");
				_Temperature.writeMED(allFields+"_Temperature");
				_DensiteLiquide.writeMED(allFields+"_LiquidDensity");
				_DensiteVapeur.writeMED(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeMED(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeMED(allFields+"_SteamEnthalpy");
				_VitesseX.writeMED(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ");
				}
				break;
			case CSV :
				_VoidFraction.writeCSV(allFields+"_VoidFraction");
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Concentration.writeCSV(allFields+"_Concentration");
				_mixtureDensity.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
				_DensiteLiquide.writeCSV(allFields+"_LiquidDensity");
				_DensiteVapeur.writeCSV(allFields+"_SteamDensity");
				_EnthalpieLiquide.writeCSV(allFields+"_LiquidEnthalpy");
				_EnthalpieVapeur.writeCSV(allFields+"_SteamEnthalpy");
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

	if (_restartWithNewFileName)
		_restartWithNewFileName=false;
}

void DriftModel::testConservation()
{
	double SUM, DELTA, x;
	int I;
	for(int i=0; i<_nVar; i++)
	{
		{
			if(i == 0)
				cout << "Masse totale (kg): ";
			else if(i == 1)
				cout << "Masse partielle 1 (kg): ";
			else
			{
				if(i == _nVar-1)
					cout << "Energie totale (J): ";
				else
					cout << "Quantite de mouvement direction "<<i-1<<" (kg.m.s^-1): ";
			}
		}
		SUM = 0;
		DELTA = 0;
		I=i;
		for(int j=0; j<_Nmailles; j++)
		{
			if(!_usePrimitiveVarsInNewton)
				VecGetValues(_old, 1, &I, &x);//on recupere la valeur du champ
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
			cout << " a une somme nulle,  variation absolue: " << fabs(DELTA) << endl;
	}
}

