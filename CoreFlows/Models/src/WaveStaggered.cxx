/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"

using namespace std;


WaveStaggered::WaveStaggered(phaseType fluid, int dim){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1
	_fluides.resize(1);

	// To do : delete EOS ?
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


void WaveStaggered::initialize(){
	cout<<"\n Initialising the Wave System model\n"<<endl;
	*_runLogFile<<"\n Initialising the Wave Sytem model\n"<<endl;

	_globalNbUnknowns = (_nVar-1)*_Nmailles + _Nfaces;//Staggered discretisation : velocity is on faces

	if(!_initialDataSet)
	{
		*_runLogFile<<"!!!!!!!!WaveStaggered::initialize() set initial data first"<<endl;
		_runLogFile->close();
		throw CdmathException("!!!!!!!!WaveStaggered::initialize() set initial data first");
	}
	cout << "Number of Phases = " << _nbPhases << " mesh dimension = "<<_Ndim<<" number of variables = "<<_nVar<<endl;
	*_runLogFile << "Number of Phases = " << _nbPhases << " spaceDim= "<<_Ndim<<" number of variables= "<<_nVar<<endl;

	_vec_normal = new double[_Ndim];
	double kappa = 3;
	double rho0 = 2; //  TODO : à définir dans le cas général
	double c = 4;
	double d = 1/2.


	//primitive field used only for saving results
	_VV=Field ("Primitive vec", FACES, _mesh, 1); //TODO comment fonctionnent Field ?

	//Construction des champs primitifs initiaux comme avant dans ParaFlow
	double * initialFieldPrim = new double[_globalNbUnknowns];
	for(int i =0; i<_Nmailles;i++)
		for(int j =0; j<_nVar;j++)
			initialFieldPrim[i*_nVar+j]=_VV(i,j); //  TODO : à corriger _VV doit être un field de size glbal unkonw et non nVar*nmailles


	/**********Petsc structures:  ****************/

	// creation des matrices en explicite et en implicite : 
	Mat B, Btopo; 
	// matrice Q tq U^n+1 = (Id + dt V^-1 Q)U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _Q); 
	MatSetSize(_Q, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatZeroEntries(_Q);
	// matrice des Inverses Volumes
	MatCreate(PETSC_COMM_SELF, & _InvVol); 
	MatSetSize(_InvVol, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatZeroEntries(_InvVol);
	// matrice |K|div(u)
	MatCreate(PETSC_COMM_SELF, & B); 
	MatSetSize(B, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatZeroEntries(B);
	
	if(_timeScheme == Implicit)
		MatCreateSeqBAIJ(PETSC_COMM_SELF, _nVar, _nVar*_Nmailles, _nVar*_Nmailles, (1+_neibMaxNbCells), PETSC_NULL, &_A);
	if(_timeScheme == Explicit){
		Mat _InvSurface, Laplacian, DivTilde, GraDivTilde;

		// matrice inverse des surfaces du bord des cellules (pour opérateur) div tilde
		MatCreate(PETSC_COMM_SELF, & _InvSurface); 
		MatSetSize(_InvSurface, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
		// matrice |K|div(u) sans |sigma| dans div(u) -> pour définir facilement le laplacien à partir de l'opérateur divergence
		MatCreate(PETSC_COMM_SELF, & Btopo); 
		MatSetSize(Btopo, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces ); 
		//  matrice gradient de pression, qui contiendra les conditions aux limites de presion
		MatCreate(PETSC_COMM_SELF, & Btpressure); 
		MatSetSize(Btpressure, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles ); 
		//  matrice gradient de div tilde, qui contiendra les conditions aux limites div tilde
		MatCreate(PETSC_COMM_SELF, & Btdiv); 
		MatSetSize(Btdiv, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles ); 
		// Matrice laplacien de pression :
		MatCreate(PETSC_COMM_SELF, & Laplacian); 
		MatSetSize(Laplacian, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles );
		// Matrice DivTilde :=  |K|/|dK|div()_h
		MatCreate(PETSC_COMM_SELF, & DivTilde); 
		MatSetSize(DivTilde, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
		// Matrice GradDivTil := [| divTild |]_\sigma
		MatCreate(PETSC_COMM_SELF, & DivTilde); 
		MatSetSize(DivTilde, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );

		// TODO : penser à détruire ces matrices 
		MatZeroEntries(_InvSurface);
		MatZeroEntries(Btopo);

		// Remplissage des matrices car indépendantes du temps
		for (int j=0; j<nbFaces;j++){
			Fj = _mesh.getFace(j);
			_isBoundary=Fj.isBorder();
			idCells = Fj.getCellsId();
			double D_sigma;
			double Perimeter1 = 0;
			double Perimeter2 = 0;
			double FaceArea;
		
			if (Fj.getNumberOfCells()==2 ){	// Fj is inside the domain and has two neighours (no junction)
				// compute the normal vector corresponding to face j : from Ctemp1 to Ctemp2
				Ctemp1 = _mesh.getCell(idCells[0]);//origin of the normal vector
				Ctemp2 = _mesh.getCell(idCells[1]);
				if (_Ndim = 1){
					if(!_sectionFieldSet)
					{
						if (Fj.x()<Ctemp1.x())
							_vec_normal[0] = -1;
						else
							_vec_normal[0] = 1;
					}
					else
					{
						if(idCells[0]>idCells[1])
							_vec_normal[0] = -1;
						else
							_vec_normal[0] = 1;
					}
					PetscReal det = Ctemp2.x() - Ctemp1.x();
					D_sigma = PetscAbsReal(det);
					Perimeter1 = 1;
					Perimeter2 = 1;
					FaceArea = 1;
				} 
				else{
					for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
						sigma1 = _mesh.getFace( Ctemp1.getFacesId()[l] );
						Perimeter1 += sigma1.getMeasure();
							if (j == Ctemp1.getFacesId()[l]){
								for (int idim = 0; idim < _Ndim; ++idim)
									_vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
								break;
							}
						}
					for(int l=0; l<Ctemp2.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
						sigma2 = _mesh.getFace( Ctemp2.getFacesId()[l] );
						Perimeter2 += sigma2.getMeasure();
						}
						
					FaceArea = Fj.getMeasure();

					if (_Ndim = 2){
						std::vector< int > nodes =  Fj.getNodesId();
						vertex = _mesh.getNode( nodes[0] );
						PetscReal det = (Ctemp1.x() - vertex.x() )* (Ctemp2.y() - vertex.y() ) - (Ctemp1.y() - vertex.y() )* (Ctemp2.x() - vertex.x() );
						// determinant of the vectors forming the diamond cell around the face sigma
					}
					//  TODO : 3D case
					D_sigma = PetscAbsReal(det);
					MatSetValues(B, 1, idCells[0], 1, j, FaceArea, ADD_VALUES ); 
					MatSetValues(B, 1, idCells[1], 1, j, -FaceArea, ADD_VALUES );  //sign minus because the exterior normal to Ctemp2 is the opposite of the Ctemp1's one

					MatSetValues(Btpressure, 1, j 1, idCells[0], FaceArea, ADD_VALUES ); 
					MatSetValues(Btpressure, 1, j 1, idCells[1], -FaceArea, ADD_VALUES ); 

					MatSetValues(Btdiv, 1, j 1, idCells[0], FaceArea, ADD_VALUES ); 
					MatSetValues(Btdiv, 1, j 1, idCells[1], -FaceArea, ADD_VALUES );

					MatSetValues(Btopo,1, idCells[0], 1, j, 1, ADD_VALUES); 
					MatSetValues(Btopo,1, idCells[0], 1, j, -1, ADD_VALUES ); 

					MatSetValues(_InvSurface,1, idCells[0],1, idCells[0], 1/(Perimeter1*Ctemp1.getNumberOfFaces()), ADD_VALUES );
					MatSetValues(_InvSurface,1, idCells[1],1, idCells[1], 1/(Perimeter2*Ctemp2.getNumberOfFaces()), ADD_VALUES );

					MatSetValues(_InvVol, 1, idCells[0],1 ,idCells[0], 1/( Ctemp1.getMeasure()* Ctemp1.getNumberOfFaces()), ADD_VALUES );
					MatSetValues(_InvVol, 1, idCells[1],1 ,idCells[1], 1/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces()), ADD_VALUES );
					MatSetValues(_InvVol, 1, _Nmailles + j, 1, _Nmailles + j,  1.0/D_sigma, ADD_VALUES); 	
				}		
			}
			else
			{
				Cint = _mesh.getCell(idCells[0]);// TODO : vérifier que l'on obtient la cellule intérieure
				if (_Ndim = 1){
					_vec_normal[0] = 1;
					PetscReal det = Fj.x() - Cint.x();
					D_sigma = PetscAbsReal(det);
					Perimeter1 = 1;
					FaceArea = 1;
				} 
				else{
					for(int l=0; l<Cint.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
						sigma1 = _mesh.getFace( Cint.getFacesId()[l] );
						Perimeter1 += sigma1.getMeasure();
							if (j == Cint.getFacesId()[l]){
								for (int idim = 0; idim < _Ndim; ++idim)
									_vec_normal[idim] = Cint.getNormalVector(l,idim);
								break;
							}
						}
					FaceArea = Fj.getMeasure();

					if (_Ndim = 2){
						std::vector< int > nodes =  Fj.getNodesId();
						vertex1 = _mesh.getNode( nodes[0] );
						vertex2 = _mesh.getNode( nodes[1] );
						PetscReal det = (Cint.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (vertex2.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
						// determinant of the vectors forming the interior half diamond cell around the face sigma
					}
					//  TODO : 3D case
				}
				D_sigma = PetscAbsReal(det);	
				MatSetValues(B, 1, idCells[0], 1, j, FaceArea, ADD_VALUES ); //TODO : vérifier l'orientation
				MatSetValues(Btopo,1, idCells[0], 1, j, 1, ADD_VALUES); 
				MatSetValues(_InvSurface,1, idCells[0],1, idCells[0], 1/(Perimeter1*Cint.getNumberOfFaces()), ADD_VALUES ),
				MatSetValues(_InvVol, 1, idCells[0],1 ,idCells[0], 1/( Cint.getMeasure()* Cint.getNumberOfFaces()), ADD_VALUES );
				MatSetValues(_InvVol, 1, _Nmailles + j, 1, _Nmailles + j,  1.0/D_sigma, ADD_VALUES); 	
	
				double pInt = primitiveVars[idCells[0]];
				double pExt = _VV(idCells[0]); // TODO : à modifirer pour pouvoir imposer conditions aux limites
				MatSetValues(Btdiv, 1, j, 1, idCells[0], 0, INSERT_VALUES );  // TODO à vérifier 
				MatSetValues(Btpressure, 1, j, 1, idCells[0], -FaceArea * Pext/Pint , ADD_VALUES );  // TODO : orientation ok ?
			 
			}	
		}
		int *indices1 = new int[_Nmailles];
		std::iota(indices1, indices1 +_Nmailles, 0);
		int *indices2 = new int[_nFaces];
		std::iota(indices2, indices2 +_nFaces, _Nmailles);

		MatMatMult(Btopo, Btpressure, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Laplacian); 
		MatMatMult(_InvSurface, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DivTilde); 
		MatMatMult(Btdiv, DivTilde, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &GraDivTilde); 

		MatSetValuesBlocked(_Q, _nMailles, indices1, _Nmailles, indices1, -d*c* Laplacian , INSERT_VALUES);  
		MatSetValuesBlocked(_Q, _nMailles, indices1, _nFaces, indices2, B/rho0 , INSERT_VALUES);
		MatSetValuesBlocked(_Q, _nFaces, indices2, _Nmailles, indices1, -kappa * Bt, INSERT_VALUES);
		MatSetValuesBlocked(_Q, _nFaces, indices2, _nFaces, indices2, -d*c * GradDivTilde , INSERT_VALUES);

		MatMatMult(_InvVol, _Q, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); // Faut-il créer la matrice Prod ?
		MatCopy(Prod,_Q, SAME_NONZERO_PATTERN); // TODO : SAME_NONZERO_PATTERN ?
		delete[] indices1, indices2;
	}

	//creation des vecteurs
	VecCreate(PETSC_COMM_SELF, & _primitiveVars);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_primitiveVars, PETSC_DECIDE, _globalNbUnknowns);
	VecSetFromOptions(_primitiveVars);
	VecDuplicate(_primitiveVars, &_old);//Old primitive variables at time step n
	VecDuplicate(_primitiveVars, &_newtonVariation);//Newton variation Uk+1-Uk 
	VecDuplicate(_primitiveVars, &_b);//Right hand side of Newton method

	// transfer information de condition initial vers primitiveVars
	int *indices = new int[_globalNbUnknowns];
	std::iota(indices, indices + _globalNbUnknowns, 0);
	VecSetValues(_primitiveVars, _globalNbUnknowns, indices, initialFieldPrim, INSERT_VALUES); 
	VecAssemblyBegin(_primitiveVars);
	VecAssemblyEnd(_primitiveVars);
	VecCopy(_primitiveVars, _old);
	VecAssemblyBegin(_old);
	VecAssemblyEnd(_old);
	if(_system)
	{
		cout << "Variables primitives initiales : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}

	delete[] initialFieldPrim, indices, indices1, indices2;

	createKSP();
	PetscPrintf(PETSC_COMM_WORLD,"SOLVERLAB Newton solver ");
	*_runLogFile << "SOLVERLAB Newton solver" << endl;
	_runLogFile->close();

	_initializedMemory=true;
	save();//save initial data
}


double WaveStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "WaveStaggered::computeTimeStep : Début calcul matrice implicite et second membre"<<endl;
		cout << endl;
	}
	if(_timeScheme == Implicit)
		MatZeroEntries(_A);
	if (_timeScheme == Explicit){
		VecAssemblyBegin(_b);
		MatMult(_Q,_old, _b); //TODO : _old = U^n ?
		VecAssemblyEnd(_b); //TODO : à quoi sert VecAssembly ?
		double maxSigma =0; 
		double minCell = 0;
		for (int j=0; j < _globalNbUnknowns; j++){
			if (j < _Nmailles){
				if (minCell > 1.0/_InvVol[j,j]){ // Primal cells
					minCell = 1.0/_InvVol[j,j];
				}
				if  (maxSigma < 1.0/_InvSurface[j,j] ){ // Perimeters
					maxSigma = 1.0/_InvSurface[j,j];
				}
			}
			else{
				if (minCell > 1.0/_InvVol[j,j]){ //Dual Cells 
					minCell = 1.0/_InvVol[j,j];
				}
			}
		}
		
	}
	
	stop=false;

	return _cfl * minCell * / (maxSigma * _c)  //TODO calculer vp de _InvVol et _InvSurf

}

bool WaveStaggered::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		_maxvp=0.;
		computeTimeStep(stop);//This compute timestep is just to update the linear system. The time step was imposed before starting the Newton iterations
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

		KSPConvergedReason reason;
		KSPGetConvergedReason(_ksp,&reason);

		if(reason<0)//solving the linear system failed
		{
			if( reason == -3)
			    cout<<"Maximum number of iterations "<<_maxPetscIts<<" reached"<<endl;
			else if( reason == -11)
			    cout<<"!!!!!!! Construction of preconditioner failed !!!!!!"<<endl;
			else if( reason == -5)
				cout<<"!!!!!!! Generic breakdown of the linear solver (Could be due to a singular matrix or preconditioner) !!!!!!"<<endl;
			else
			{
			    cout<<"PETSc divergence reason  "<< reason <<endl;
				cout<<"Nombre d'itérations effectuées "<< _PetscIts<<" nombre maximal Itérations autorisées "<<_maxPetscIts<<endl;
			}
			*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Raison PETSC numéro "<<reason<<endl;
			*_runLogFile<<"Nombre d'itérations effectuées "<< _PetscIts<<" nombre maximal Itérations autorisées "<<_maxPetscIts<<endl;
			converged=false;
			return false;
		}
		else{//solving the linear system succeeded
			//Calcul de la variation relative Uk+1-Uk ou Vkp1-Vk
			_erreur_rel = 0.;
			double x, dx;
			int I;
			for(int j=1; j<=_Nmailles; j++)
			{
				for(int k=0; k<_nVar; k++)
				{
					I = (j-1)*_nVar + k;
					VecGetValues(_newtonVariation, 1, &I, &dx);
					if( !_usePrimitiveVarsInNewton)
						VecGetValues(_conservativeVars, 1, &I, &x);
					else
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

	//Change the relaxation coefficient to ease convergence
	double relaxation=1;
	VecAXPY(_primitiveVars,     relaxation, _newtonVariation);//Vk+1=Vk+relaxation*deltaV
	if(_system)
	{
		cout<<"Vecteur Vkp1-Vk "<<endl;
		VecView(_newtonVariation,  PETSC_VIEWER_STDOUT_SELF);
		cout << "Nouvel etat courant Vk de l'iteration Newton: " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_SELF);
	}


	return true;
}

void WaveStaggered::validateTimeStep()
{
	if(_system)
	{
		cout<<" Vecteur Un"<<endl;
		VecView(_old,  PETSC_VIEWER_STDOUT_WORLD);
		cout<<" Vecteur Un+1"<<endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
	}
	VecAXPY(_old,  -1, _primitiveVars);//old contient old-courant

	//Calcul de la variation Un+1-Un
	_erreur_rel= 0;
	double x, dx;

	for(int j=1; j<=_globalNbUnknowns; j++){
		VecGetValues(_old, 1, &j, &dx);
		VecGetValues(_conservativeVars, 1, &j, &x);
		if (fabs(x)< _precision)
		{
			if(_erreur_rel < fabs(dx))
				_erreur_rel = fabs(dx);
		}
		else if(_erreur_rel < fabs(dx/x))
			_erreur_rel = fabs(dx/x);
	}

	_isStationary =_erreur_rel <_precision;
	VecCopy(_primitiveVars, _old);

	_time+=_dt;
	_nbTimeStep++;
	if (_nbTimeStep%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
		save();
}

void WaveStaggered::computeNewtonVariation()
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
		if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
		{
			cout<<"Vecteur _newtonVariation =_b*dt"<<endl;
			VecView(_newtonVariation,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
	}
}

void WaveStaggered::convectionState( const long &i, const long &j, const bool &IsBord)
{}

void WaveStaggered::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	int k;
	double v2=0, q_n=0;//q_n=quantité de mouvement normale à la face frontière;
	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+1)]*normale[k];

	double porosityj=_porosityField(j);

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
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
		else if((_nbTimeStep-1)%_freqSave ==0)
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
			if((_nbTimeStep-1)%_freqSave ==0)
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
		if(q_n<=0 &&  (_nbTimeStep-1)%_freqSave ==0)
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

		if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
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

void WaveStaggered::convectionMatrices()
{}

void WaveStaggered::computeScaling(double maxvp)
{}

void WaveStaggered::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{}

void WaveStaggered::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{}

void WaveStaggered::porosityGradientSourceVector()
{}

void WaveStaggered::jacobian(const int &j, string nameOfGroup,double * normale)
{}

Vector WaveStaggered::convectionFlux(Vector U,Vector V, Vector normale, double porosity)
{}

void WaveStaggered::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
{}

void WaveStaggered::getDensityDerivatives( double pressure, double temperature, double v2)
{}

void WaveStaggered::terminate(){ 
	VecDestroy(&_newtonVariation);
	VecDestroy(&_b);
	VecDestroy(&_primitiveVars);


	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);
	for(int i=0;i<_nbPhases;i++)
		delete _fluides[i];

	// Destruction du solveur de Newton de PETSc
	if(_timeScheme == Implicit && _nonLinearSolver != Newton_SOLVERLAB)
		SNESDestroy(&_snes);
}

void WaveStaggered::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/WaveStaggeredPrim_");///Results
	string cons(_path+"/WaveStaggeredCons_");
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
