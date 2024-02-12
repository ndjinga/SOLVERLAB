/*
 * WaveStaggered.cxx
 */

#include "../inc/WaveStaggered.hxx"
#include "StiffenedGas.hxx"

using namespace std;


WaveStaggered::WaveStaggered(phaseType fluid, int dim, double kappa, double rho){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1;
	_fluides.resize(1);
	_kappa = kappa;
	_rho = rho;
	_c = sqrt(kappa/rho);
	// TODO : appeler constructeur pb fluid ?
}

void WaveStaggered::setInitialField(const Field &field)
{
	if(_Ndim != field.getSpaceDimension()){
		*_runLogFile<<"WaveStaggered::setInitialField: mesh has incorrect space dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered::setInitialField: mesh has incorrect space dimension");
	}
	if  (field.getName() == "pressure"){ //TODO : field->getName() ? , penser à donner un nom à l'initialisation dans cas test 
		_Pressure = field;
		_Pressure.setName("Pressure results");
		_time=_Pressure.getTime();
		_mesh=_Pressure.getMesh();
	}
	else {
		_Velocity = field;
		_Velocity.setName("Velocity results");
		_time=_Velocity.getTime();
		_mesh=_Velocity.getMesh();

	}

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
	
    _neibMaxNbCells=_mesh.getMaxNbNeighbours(CELLS);
	
	/*** MPI distribution of parameters ***/
	MPI_Allreduce(&_initialDataSet, &_initialDataSet, 1, MPIU_BOOL, MPI_LOR, PETSC_COMM_WORLD);
	
	int nbVoisinsMax;
	MPI_Bcast(&_Nmailles      , 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&_neibMaxNbCells, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	nbVoisinsMax = _neibMaxNbCells;
	
    _d_nnz = (nbVoisinsMax+1)*_nVar;
    _o_nnz =  nbVoisinsMax   *_nVar;
}

void WaveStaggered::initialize(){
	cout<<"\n Initialising the Wave System model\n"<<endl;
	*_runLogFile<<"\n Initialising the Wave Sytem model\n"<<endl;

	_globalNbUnknowns = _Nmailles + _Nfaces;//Staggered discretisation : velocity is on faces

	if(!_initialDataSet)
	{
		*_runLogFile<<"!!!!!!!!WaveStaggered::initialize() set initial data first"<<endl;
		_runLogFile->close();
		throw CdmathException("!!!!!!!!WaveStaggered::initialize() set initial data first");
	}
	cout << "Number of Phases = " << _nbPhases << " mesh dimension = "<<_Ndim<<" number of variables = "<<_nVar<<endl;
	*_runLogFile << "Number of Phases = " << _nbPhases << " spaceDim= "<<_Ndim<<" number of variables= "<<_nVar<<endl;

	_vec_normal = new double[_Ndim];
	_d = 1/(2* sqrt(_neibMaxNbCells) );


	//primitive field used only for saving results
	_Velocity=Field ("Primitive vec", FACES, _mesh, 1); //TODO comment fonctionnent Field ?
	_Pressure=Field ("Primitive vec", CELLS, _mesh, 1); //TODO comment fonctionnent Field ?

	//Construction des champs primitifs initiaux comme avant dans ParaFlow
	double * initialFieldPrim = new double[_globalNbUnknowns];
	for(int i =0; i<_globalNbUnknowns; i++)
		initialFieldPrim[i]=_VV(i); //  TODO : à corriger _VV doit être un field de size glbal unkonw et non nVar*nmailles

	/**********Petsc structures:  ****************/

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

	delete[] initialFieldPrim, indices;

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
	double maxPerim = 0; 
	double minCell = 0;
	
	if (_timeScheme == Explicit && _dt == 0){ // The matrices are assembled only in the first time step since linear problem
		Mat B, Btopo, Btpressure, Bt, InvSurface, InvVol;; 
		// matrice Q tq U^n+1 = (Id + dt V^-1 Q)U^n pour schéma explicite
		MatCreate(PETSC_COMM_SELF, & _Q); 
		MatSetSizes(_Q, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
		MatZeroEntries(_Q);

		// matrice des Inverses Volumes
		MatCreate(PETSC_COMM_SELF, & InvVol); 
		MatSetSizes(InvVol, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
		MatZeroEntries(InvVol);
		// matrice des Inverses de Surfaces
		MatCreate(PETSC_COMM_SELF, & InvSurface); 
		MatSetSizes(InvSurface, PETSC_DECIDE, PETSC_DECIDE, _Nmailles , _Nmailles );
		MatZeroEntries(InvSurface);

		// matrice DIVERGENCE (|K|div(u))
		MatCreate(PETSC_COMM_SELF, & B); 
		MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
		MatZeroEntries(B);
		// matrix GRADIENT (we will impose to be 0 on boundary faces so that u^n+1 = u^n at the boundary)
		MatCreate(PETSC_COMM_SELF, & Bt); 
		MatSetSizes(Bt, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
		MatZeroEntries(Bt);

		// matrice Btopo = {|K|div(u) sans |sigma| dans div(u)}-> pour définir facilement le laplacien à partir de l'opérateur divergence
		MatCreate(PETSC_COMM_SELF, & Btopo); 
		MatSetSizes(Btopo, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces ); 
		MatZeroEntries(Btopo);

		// matrice gradient de pression, qui contiendra les conditions aux limites de pression
		MatCreate(PETSC_COMM_SELF, & Btpressure); 
		MatSetSizes(Btpressure, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles ); 
		MatZeroEntries(Btpressure);
		
		// Assembly of matrices 
		for (int j=0; j<_Nfaces;j++){
			Face Fj = _mesh.getFace(j);
			bool _isBoundary=Fj.isBorder();
			std::vector< int > idCells = Fj.getCellsId();
			PetscScalar det, FaceArea;
			PetscInt IndexFace = _Nmailles + j;
		
			if (Fj.getNumberOfCells()==2 ){	// Fj is inside the domain and has two neighours (no junction)
				// compute the normal vector corresponding to face j : from Ctemp1 to Ctemp2
				Cell Ctemp1 = _mesh.getCell(idCells[0]);//origin of the normal vector
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
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
					det = Ctemp2.x() - Ctemp1.x();
					FaceArea = 1;
				} 
				else{
					for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
							if (j == Ctemp1.getFacesId()[l]){
								for (int idim = 0; idim < _Ndim; ++idim)
									_vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
								break;
							}
						}

					if (_Ndim = 2){
						std::vector< int > nodes =  Fj.getNodesId();
						Node vertex = _mesh.getNode( nodes[0] );
						// determinant of the vectors forming the diamond cell around the face sigma
						// TODO : vérifier déterminant
						det = (Ctemp1.x() - vertex.x() )* (Ctemp2.y() - vertex.y() ) - (Ctemp1.y() - vertex.y() )* (Ctemp2.x() - vertex.x() );
					}
					if (_Ndim = 3){	
						cout<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
						*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
						_runLogFile->close();
						throw CdmathException("WaveStaggered pas dispo en 3D, arret de calcul");				
					}
					
					//  TODO : vérifier orientation, ne vaut-il pas définir un vecteur n_sigma pour chaque face ?
					FaceArea = Fj.getMeasure();
					PetscScalar MinusFaceArea = -FaceArea;
					PetscScalar InvD_sigma = 1.0/PetscAbsReal(det);
					PetscScalar InvVol1 = 1/( Ctemp1.getMeasure()* Ctemp1.getNumberOfFaces());
					PetscScalar InvVol2 = 1/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());
					PetscScalar InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
					PetscScalar InvPerimeter2 = 1/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
					PetscScalar One=1;
					PetscScalar MinusOne=-1;


					MatSetValues(B, 1, &idCells[0], 1, &j, &FaceArea, ADD_VALUES ); 
					MatSetValues(B, 1, &idCells[1], 1, &j, &MinusFaceArea, ADD_VALUES );  
					MatSetValues(Btopo,1, &idCells[0], 1, &j, &One, ADD_VALUES); 
					MatSetValues(Btopo,1, &idCells[0], 1, &j, &MinusOne, ADD_VALUES ); 
					//sign minus because the exterior normal to Ctemp2 is the opposite of the Ctemp1's one

					MatSetValues(Btpressure, 1, &j, 1, &idCells[0], &FaceArea, ADD_VALUES ); 
					MatSetValues(Btpressure, 1, &j, 1, &idCells[1], &MinusFaceArea, ADD_VALUES ); 
					MatSetValues(Bt, 1, &j, 1, &idCells[0], &FaceArea, ADD_VALUES ); 
					MatSetValues(Bt, 1, &j, 1, &idCells[1], &MinusFaceArea, ADD_VALUES ); 

					MatSetValues(InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
					MatSetValues(InvSurface,1, &idCells[1],1, &idCells[1], &InvPerimeter2, ADD_VALUES );

					MatSetValues(InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1 , ADD_VALUES );
					MatSetValues(InvVol, 1, &idCells[1],1 ,&idCells[1], &InvVol2, ADD_VALUES );
					MatSetValues(InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	
				}		
			}
			else // boundary faces
			{
				Cell Cint = _mesh.getCell(idCells[0]);// TODO : vérifier que l'on obtient la cellule intérieure
				if (_Ndim = 1){
					_vec_normal[0] = 1;
					PetscScalar det = Fj.x() - Cint.x();
					FaceArea = 1;
				} 
				else{
					for(int l=0; l<Cint.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
							if (j == Cint.getFacesId()[l]){
								for (int idim = 0; idim < _Ndim; ++idim)
									_vec_normal[idim] = Cint.getNormalVector(l,idim);
								break;
							}
					}

					std::vector< int > nodes =  Fj.getNodesId();
					Node vertex1 = _mesh.getNode( nodes[0] );
					Node vertex2 = _mesh.getNode( nodes[1] );
					PetscScalar det = (Cint.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (vertex2.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
					// determinant of the vectors forming the interior half diamond cell around the face sigma
					//  TODO : 3D case
				}
				PetscScalar One=1;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar InvD_sigma = 1.0/PetscAbsReal(det);
				PetscScalar InvVol1 = 1/( Cint.getMeasure()* Cint.getNumberOfFaces());
				PetscScalar InvPerimeter1 = 1/( _perimeters(idCells[0])*Cint.getNumberOfFaces()  );
	
				MatSetValues(B, 1, &idCells[0], 1, &j, &FaceArea, ADD_VALUES ); //TODO : vérifier l'orientation
				MatSetValues(Btopo,1, &idCells[0], 1, &j, &One, ADD_VALUES); 
				MatSetValues(InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES ),
				MatSetValues(InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
				MatSetValues(InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	

				PetscScalar pInt;
				VecGetValues(_primitiveVars, 1, &idCells[0], &pInt);
				double pExt = _VV(idCells[0],j); // TODO : à modifier pour pouvoir imposer conditions aux limites
				MatSetValues(Bt, 1, &j, 1, &idCells[0], 0, INSERT_VALUES );  // TODO à vérifier 
				PetscScalar pressureGrad =  1 -FaceArea * pExt/pInt;
				MatSetValues(Btpressure, 1, &j, 1, &idCells[0], &pressureGrad, ADD_VALUES );  // TODO : orientation ok ?
			 
			}	
		}
		MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(Btopo,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Btopo, MAT_FINAL_ASSEMBLY);

		MatAssemblyBegin(Bt, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Bt, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(Btpressure, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Btpressure, MAT_FINAL_ASSEMBLY);

		MatAssemblyBegin(InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(InvVol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(InvVol, MAT_FINAL_ASSEMBLY);
		
		int *indices1 = new int[_Nmailles];
		std::iota(indices1, indices1 +_Nmailles, 0);
		int *indices2 = new int[_Nfaces];
		std::iota(indices2, indices2 +_Nfaces, _Nmailles);

		Mat Laplacian, GradDivTilde;
		MatMatMult(Btopo, Btpressure, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Laplacian);  
		MatMatMatMult(Bt,InvSurface, B , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &GradDivTilde); 

		PetscScalar minusdc = -_d*_c;
		PetscScalar minuskappa = -_kappa;
		PetscScalar invrho = 1.0/_rho;
		MatScale(Laplacian, minusdc );
		MatScale(B, invrho);
		MatScale(Bt, minuskappa);
		MatScale(GradDivTilde, minusdc );

		MatSetValuesBlockedLocal(_Q, _Nmailles, indices1, _Nmailles, indices1, Laplacian , INSERT_VALUES);  
		MatSetValuesBlockedLocal(_Q, _Nmailles, indices1, _Nfaces, indices2, B , INSERT_VALUES);
		MatSetValuesBlockedLocal(_Q, _Nfaces, indices2, _Nmailles, indices1, Bt, INSERT_VALUES);
		MatSetValuesBlockedLocal(_Q, _Nfaces, indices2, _Nfaces, indices2, GradDivTilde , INSERT_VALUES); // TODO : MATAssmebly for _Q ?

		Mat Prod;
		MatMatMult(InvVol, _Q, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_Q, SAME_NONZERO_PATTERN); // TODO : SAME_NONZERO_PATTERN ?

		if (_cfl > _d/2.0){
			cout << "cfl ="<< _cfl <<"is to high, cfl is updated to 0.9*_d/2" << endl;
			_cfl =  0.9 * _d/2.0;
		}
		VecAssemblyBegin(_b);
		MatMult(_Q,_old, _b); //TODO : _old = U^n ?
		VecAssemblyEnd(_b); 

		Vec V, W;
		MatGetDiagonal(InvVol,V);
		MatGetDiagonal(InvSurface, W);
		int *indices3 = new int[_Nmailles + _Nfaces];
		std::iota(indices3, indices3 +_Nmailles + _Nfaces, 0);
		int *indices4 = new int[_Nmailles];
		std::iota(indices4, indices4 +_Nmailles, 0);
		PetscScalar minInvSurf, maxInvVol;
		VecMax(W, indices3, &maxInvVol);
		VecMin(V, indices4, &minInvSurf);
		maxPerim = 1.0/minInvSurf;
		minCell = 1.0/maxInvVol;

		delete[] indices1, indices2, indices3, indices4;
		MatDestroy(& B); 
		MatDestroy(& Btopo); 
		MatDestroy(& Btpressure); 
		MatDestroy(& Bt); 
		MatDestroy(& InvVol); 
		MatDestroy(& InvSurface);
		MatDestroy(& Laplacian);
		MatDestroy(& GradDivTilde); 
	}
	if (_timeScheme == Explicit){
		VecAssemblyBegin(_b);
		MatMult(_Q,_old, _b); //TODO : _old = U^n ?
		VecAssemblyEnd(_b); 

	}
	stop=false;

	return _cfl * minCell / (maxPerim * _c) ; 

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

// void WaveStaggered::computeScaling(double maxvp)
// {}

void WaveStaggered::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{}

void WaveStaggered::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{}

void WaveStaggered::porosityGradientSourceVector()
{}

void WaveStaggered::jacobian(const int &j, string nameOfGroup,double * normale)
{}

Vector WaveStaggered::convectionFlux(Vector U,Vector V, Vector normale, double porosity)
{
	return U;
}

void WaveStaggered::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
{}

void WaveStaggered::getDensityDerivatives( double pressure, double temperature, double v2)
{}

void WaveStaggered::terminate(){ 
	VecDestroy(&_newtonVariation);
	VecDestroy(&_b);
	VecDestroy(&_primitiveVars);
	MatDestroy(& _Q); 
 
	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);
	for(int i=0;i<_nbPhases;i++)
		delete _fluides[i];

	// Destruction du solveur de Newton de PETSc
	if(_timeScheme == Implicit && _nonLinearSolver != Newton_SOLVERLAB)
		SNESDestroy(&_snes);
}

// void WaveStaggered::save(){
//     PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
//     *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

// 	string prim(_path+"/WaveStaggeredPrim_");///Results
// 	prim+=_fileName;

// 	PetscInt Ii;
// 	for (long i = 0; i < ; i++){
// 			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
// 		}
// 	}
// 	_Pressure.setTime(_time,_nbTimeStep);
// 	_Velocity.setTime(_time,_nbTimeStep);

// 	// create mesh and component info
// 	if (_nbTimeStep ==0){
// 		string prim_suppress ="rm -rf "+prim+"_*";
// 		string cons_suppress ="rm -rf "+cons+"_*";

// 		system(prim_suppress.c_str());//Nettoyage des précédents calculs identiques
// 		system(cons_suppress.c_str());//Nettoyage des précédents calculs identiques


// 		_VV.setInfoOnComponent(0,"Pressure_(Pa)");
// 		_VV.setInfoOnComponent(1,"Velocity_x_(m/s)");
// 		if (_Ndim>1)
// 			_VV.setInfoOnComponent(2,"Velocity_y_(m/s)");
// 		if (_Ndim>2)
// 			_VV.setInfoOnComponent(3,"Velocity_z_(m/s)");
// 		_VV.setInfoOnComponent(_nVar-1,"Temperature_(K)");

// 		switch(_saveFormat)
// 		{
// 		case VTK :
// 			_VV.writeVTK(prim);
// 			break;
// 		case MED :
// 			_VV.writeMED(prim);
// 			break;
// 		case CSV :
// 			_VV.writeCSV(prim);
// 			break;
// 		}
// 	}
// 	// do not create mesh
// 	else{
// 		switch(_saveFormat)
// 		{
// 		case VTK :
// 			_VV.writeVTK(prim,false);
// 			break;
// 		case MED :
// 			_VV.writeMED(prim,false);
// 			break;
// 		case CSV :
// 			_VV.writeCSV(prim);
// 			break;
// 		}
// 		if(_saveConservativeField){
// 			switch(_saveFormat)
// 			{
// 			case VTK :
// 				_UU.writeVTK(cons,false);
// 				break;
// 			case MED :
// 				_UU.writeMED(cons,false);
// 				break;
// 			case CSV :
// 				_UU.writeCSV(cons);
// 				break;
// 			}
// 		}
// 	}
// 	if(_saveVelocity){
// 		for (long i = 0; i < _Nmailles; i++){
// 				VecGetValues(_primitiveVars,1,&i,&_Vitesse(i));
// 			}
// 		}
// 		_Vitesse.setTime(_time,_nbTimeStep);
// 		if (_nbTimeStep ==0){
// 			_Vitesse.setInfoOnComponent(0,"Velocity_x_(m/s)");
// 			_Vitesse.setInfoOnComponent(1,"Velocity_y_(m/s)");
// 			_Vitesse.setInfoOnComponent(2,"Velocity_z_(m/s)");

// 			switch(_saveFormat)
// 			{
// 			case VTK :
// 				_Vitesse.writeVTK(prim+"_Velocity");
// 				break;
// 			case MED :
// 				_Vitesse.writeMED(prim+"_Velocity");
// 				break;
// 			case CSV :
// 				_Vitesse.writeCSV(prim+"_Velocity");
// 				break;
// 			}
// 		}
// 		else{
// 			switch(_saveFormat)
// 			{
// 			case VTK :
// 				_Vitesse.writeVTK(prim+"_Velocity",false);
// 				break;
// 			case MED :
// 				_Vitesse.writeMED(prim+"_Velocity",false);
// 				break;
// 			case CSV :
// 				_Vitesse.writeCSV(prim+"_Velocity");
// 				break;
// 			}
// 		}
// 	}
// 	if(_isStationary)
// 	{
// 		prim+="_Stat";
// 		cons+="_Stat";

// 		switch(_saveFormat)
// 		{
// 		case VTK :
// 			_VV.writeVTK(prim);
// 			break;
// 		case MED :
// 			_VV.writeMED(prim);
// 			break;
// 		case CSV :
// 			_VV.writeCSV(prim);
// 			break;
// 		}

// 		if(_saveConservativeField){
// 			switch(_saveFormat)
// 			{
// 			case VTK :
// 				_UU.writeVTK(cons);
// 				break;
// 			case MED :
// 				_UU.writeMED(cons);
// 				break;
// 			case CSV :
// 				_UU.writeCSV(cons);
// 				break;
// 			}
// 		}

// 		if(_saveVelocity){
// 			switch(_saveFormat)
// 			{
// 			case VTK :
// 				_Vitesse.writeVTK(prim+"_Velocity");
// 				break;
// 			case MED :
// 				_Vitesse.writeMED(prim+"_Velocity");
// 				break;
// 			case CSV :
// 				_Vitesse.writeCSV(prim+"_Velocity");
// 				break;
// 			}
// 		}
// 	}
// }
