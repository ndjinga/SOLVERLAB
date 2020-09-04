/*
 * SinglePhaseStaggered.cxx
 */

#include "SinglePhaseStaggered.hxx"

using namespace std;

void
computeVelocityMCells(const Field& velocity,
                      Field& velocityMCells)
{
    Mesh myMesh=velocity.getMesh();
    int nbCells=myMesh.getNumberOfCells();

    for(int i=0;i<nbCells;i++)
    {
        std::vector< int > facesId=myMesh.getCell(i).getFacesId();
        velocityMCells(i)=(velocity(facesId[0])+velocity(facesId[1]))/2.;
    }
}

SinglePhaseStaggered::SinglePhaseStaggered(phaseType fluid, pressureEstimate pEstimate, int dim, bool useDellacherieEOS){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1;
	_dragCoeffs=vector<double>(1,0);
	_fluides.resize(1);
	_useDellacherieEOS=useDellacherieEOS;

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		if(fluid==Gas){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(1.4,743,300,2.22e5);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, e=2.22e5, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(996,1e5,300,1.12e5,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C, e=1.12e5, c_v=4130
		}
	}
	else{//EOS at 155 bars and 618K 
		if(fluid==Gas){
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_fluides[0] = new StiffenedGas(102,1.55e7,618,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
		}
		else{//To do : change to normal regime: 155 bars and 573K
			cout<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			if(_useDellacherieEOS)
				_fluides[0]= new StiffenedGasDellacherie(2.35,1e9,-1.167e6,1816); //stiffened gas law for water from S. Dellacherie
			else
				_fluides[0]= new StiffenedGas(594.,1.55e7,618.,1.6e6, 621.,3100.);  //stiffened gas law for water at pressure 155 bar, and temperature 345°C
		}
	}
}
void SinglePhaseStaggered::initialize(){
	cout<<"Initialising the Navier-Stokes model"<<endl;
	*_runLogFile<<"Initialising the Navier-Stokes model"<<endl;

	_Uroe = new double[_nVar];
	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+1]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];
	if(_saveVelocity)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant

	if(_entropicCorrection)
		_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

	ProblemFluid::initialize();
}

bool SinglePhaseStaggered::iterateTimeStep(bool &converged)
{
	if(_timeScheme == Explicit || !_usePrimitiveVarsInNewton)
		ProblemFluid::iterateTimeStep(converged);
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
void SinglePhaseStaggered::computeNewtonVariation()
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
void SinglePhaseStaggered::convectionState( const long &i, const long &j, const bool &IsBord){

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

void SinglePhaseStaggered::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	int k;
	double v2=0, q_n=0;//q_n=quantité de mouvement normale à la face frontière;
	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+1)]*normale[k];

	double porosityj=_porosityField(j);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< " face unit normal vector "<<endl;
		for(k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall){
		//Pour la convection, inversion du sens de la vitesse normale
		for(k=0; k<_Ndim; k++)
			_externalStates[(k+1)]-= 2*q_n*normale[k];

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
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double pression=_externalStates[0];
		double T=_limitField[nameOfGroup].T;
		double rho=_fluides[0]->getDensity(pression,T);

		_externalStates[0]=porosityj*rho;
		_externalStates[1]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
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
			VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_fluides[0]->getDensity(pression,T);

			_externalStates[0]=porosityj*rho;
			_externalStates[1]=_externalStates[0]*(_limitField[nameOfGroup].v_x[0]);
			v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
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
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		//Building the external state
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		if(q_n<=0){
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress,_limitField[nameOfGroup].T);
		}
		else{
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature"<<endl;
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		}

		for(k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);


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
		if(q_n<=0 &&  _nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature"<<endl;

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
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		for(k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}
}

void SinglePhaseStaggered::convectionMatrices()
{
	//entree: URoe = rho, u, H
	//sortie: matrices Roe+  et Roe-

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhaseStaggered::convectionMatrices()"<<endl;

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
			cout<<"SinglePhaseStaggered::convectionMatrices Eigenvalues "<<u_n-c<<" , "<<u_n<<" , "<<u_n+c<<endl;

		RoeMatrixConservativeVariables( u_n, H,vitesse,k,K);

		/******** Construction des matrices de decentrement ********/
		if( _spaceScheme ==centered){
			if(_entropicCorrection)
			{
				*_runLogFile<<"SinglePhaseStaggered::convectionMatrices: entropy scheme not available for centered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhaseStaggered::convectionMatrices: entropy scheme not available for centered scheme");
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
				*_runLogFile<<"SinglePhaseStaggered::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhaseStaggered::convectionMatrices: entropy scheme not available for staggered scheme");
			}

			staggeredRoeUpwindingMatrixConservativeVariables( u_n, H, vitesse, k, K);
		}
		else
		{
			*_runLogFile<<"SinglePhaseStaggered::convectionMatrices: scheme not treated"<<endl;
			_runLogFile->close();
			throw CdmathException("SinglePhaseStaggered::convectionMatrices: scheme not treated");
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
		*_runLogFile<<"SinglePhaseStaggered::convectionMatrices: well balanced option not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("SinglePhaseStaggered::convectionMatrices: well balanced option not treated");
	}
}
void SinglePhaseStaggered::computeScaling(double maxvp)
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


void SinglePhaseStaggered::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
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
		cout<<"SinglePhaseStaggered::sourceVector"<<endl;
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

void SinglePhaseStaggered::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
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
		cout<<"SinglePhaseStaggered::pressureLossVector K= "<<K<<endl;
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

void SinglePhaseStaggered::porosityGradientSourceVector()
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

void SinglePhaseStaggered::jacobian(const int &j, string nameOfGroup,double * normale)
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
		throw CdmathException("SinglePhaseStaggered::jacobian: This boundary condition is not treated");
	}
}


Vector SinglePhaseStaggered::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhaseStaggered::convectionFlux start"<<endl;
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
		cout<<"SinglePhaseStaggered::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void SinglePhaseStaggered::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
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

void SinglePhaseStaggered::getDensityDerivatives( double pressure, double temperature, double v2)
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
		*_runLogFile<< "SinglePhaseStaggered::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhaseStaggered::staggeredVFFCMatricesPrimitiveVariables: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"_drho_sur_dp= "<<_drho_sur_dp<<", _drho_sur_dT= "<<_drho_sur_dT<<endl;
		cout<<"_drhoE_sur_dp= "<<_drhoE_sur_dp<<", _drhoE_sur_dT= "<<_drhoE_sur_dT<<endl;
	}
}
void SinglePhaseStaggered::save(){
	string prim(_path+"/SinglePhaseStaggeredPrim_");///Results
	string cons(_path+"/SinglePhaseStaggeredCons_");
	prim+=_fileName;
	cons+=_fileName;

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
	if (_nbTimeStep ==0){
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
	if(_saveVelocity){
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
		if (_nbTimeStep ==0){
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
	}
}
