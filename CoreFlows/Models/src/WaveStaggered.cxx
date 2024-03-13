/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include <numeric>

using namespace std;

WaveStaggered::WaveStaggered(int dim, double kappa, double rho, MPI_Comm comm):ProblemCoreFlows(comm){
	_Ndim=dim;
	_nVar = 2; // 2 equations TODO : sure ?
	_kappa = kappa;
	_rho = rho;
	_c = sqrt(kappa/rho);
	_saveVelocity=false; 
	_savePressure=false; 
	_facesBoundinit = false;
	_indexFacePeriodicSet = false;
}

std::map<int,double>  WaveStaggered::getboundaryPressure(){
		return _boundaryPressure;
}

void  WaveStaggered::setboundaryVelocity(std::map< int, double> BoundaryVelocity){
	if (_facesBoundinit == true){
		std::map<int,double>::iterator it;
		for( it= BoundaryVelocity.begin(); it != BoundaryVelocity.end(); it++){
			_Velocity( it->first ) = it->second; 
		}
	}
	else{
		*_runLogFile<<"WaveStaggered::setboundaryVelocity should be called after WaveStaggered::setInitialField(Velocity)"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered::setboundaryVelocity should be called after WaveStaggered::setInitialField(Velocity)");
	}
}

void  WaveStaggered::setboundaryPressure(std::map< int, double> BoundaryPressure){
	_boundaryPressure = BoundaryPressure;
}


void WaveStaggered::setInitialField(const Field &field)
{
	if(_Ndim != field.getSpaceDimension()){
		*_runLogFile<<"WaveStaggered::setInitialField: mesh has incorrect space dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered::setInitialField: mesh has incorrect space dimension");
	}
	if  (field.getName()  == "pressure" || field.getTypeOfField() == CELLS){  
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
		_facesBoundinit = true;

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

	_globalNbUnknowns = _Nmailles + _Nfaces; //Staggered discretisation : velocity is on faces

	if(!_initialDataSet)
	{
		*_runLogFile<<"!!!!!!!!WaveStaggered::initialize() set initial data first"<<endl;
		_runLogFile->close();
		throw CdmathException("!!!!!!!!WaveStaggered::initialize() set initial data first");
	}
	cout << "mesh dimension = "<<_Ndim <<endl;
	*_runLogFile << " spaceDim= "<<_Ndim <<endl;

	_d = 1/(2* sqrt(_neibMaxNbCells) );

	//Construction des champs primitifs initiaux comme avant dans ParaFlow
	double * initialFieldVelocity = new double[_Nfaces];
	for(int i =0; i<_Nfaces; i++)
		initialFieldVelocity[i]=_Velocity(i); 
		
	double * initialFieldPressure = new double[_Nmailles];
	for(int i =0; i<_Nmailles; i++)
		initialFieldPressure[i]=_Pressure(i); 
	/**********Petsc structures:  ****************/

	//creation des vecteurs
	VecCreate(PETSC_COMM_SELF, & _primitiveVars);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_primitiveVars, PETSC_DECIDE, _globalNbUnknowns);
	VecSetFromOptions(_primitiveVars);
	VecDuplicate(_primitiveVars, &_newtonVariation);//Newton variation Uk+1-Uk 
	VecDuplicate(_primitiveVars, &_b);//Right hand side of Newton method

	// transfer information de condition initial vers primitiveVars
	int *indices1 = new int[_Nmailles];
	int *indices2 = new int[_Nfaces];
	std::iota(indices1, indices1 + _Nmailles, 0);
	std::iota(indices2, indices2 + _Nfaces, _Nmailles);
	VecSetValues(_primitiveVars, _Nmailles , indices1, initialFieldPressure, INSERT_VALUES); 
	VecSetValues(_primitiveVars, _Nfaces, indices2, initialFieldVelocity, INSERT_VALUES); 
	VecAssemblyBegin(_primitiveVars);
	VecAssemblyEnd(_primitiveVars);

	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);
	MatZeroEntries(_A);

	// matrice des Inverses Volumes V^{-1}
	MatCreate(PETSC_COMM_SELF, &_InvVol); 
	MatSetSizes(_InvVol, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_InvVol);
	MatSetUp(_InvVol);
	MatZeroEntries(_InvVol);

	if(_system)
	{
		cout << "Variables primitives initiales : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}

	delete[] initialFieldVelocity, initialFieldPressure, indices1, indices2;

	createKSP();
	PetscPrintf(PETSC_COMM_WORLD,"SOLVERLAB Newton solver ");
	*_runLogFile << "SOLVERLAB Newton solver" << endl;
	_runLogFile->close();

	_initializedMemory=true;
	save();//save initial data
}



double WaveStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme

	//The matrices are assembled only in the first time step since linear problem
	if (_timeScheme == Explicit && _nbTimeStep == 0 ){ //TODO : pourquoi la solution exate n'évolue pas quand on enlève _nbTimeStep==0
		cout << "WaveStaggered::computeTimeStep : Début calcul matrice implicite et second membre"<<endl;
		cout << endl;
		Mat B, Bt, Laplacian, InvSurface;
		
		// matrice DIVERGENCE (|K|div(u))
		MatCreate(PETSC_COMM_SELF, & B); 
		MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
		MatSetFromOptions(B);
		MatSetUp(B);
		MatZeroEntries(B);
		
		// matrix GRADIENT (we will impose to be 0 on faces so that u^n+1 = u^n at the boundary)
		MatCreate(PETSC_COMM_SELF, & Bt); 
		MatSetSizes(Bt, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
		MatSetFromOptions(Bt);
		MatSetUp(Bt);
		MatZeroEntries(Bt);

		// matrix LAPLACIAN (without boundary terms)
		MatCreate(PETSC_COMM_SELF, & Laplacian); 
		MatSetSizes(Laplacian, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles ); 
		MatSetFromOptions(Laplacian);
		MatSetUp(Laplacian);
		MatZeroEntries(Laplacian);

		// Vector BoundaryTerms for Pressure
		VecCreate(PETSC_COMM_SELF, & _BoundaryTerms); 
		VecSetSizes(_BoundaryTerms, PETSC_DECIDE, _globalNbUnknowns); 
		VecSetFromOptions(_BoundaryTerms);
		VecSetUp(_BoundaryTerms);
		VecZeroEntries(_BoundaryTerms);

		// matrice des Inverses de Surfaces
		MatCreate(PETSC_COMM_SELF, & InvSurface); 
		MatSetSizes(InvSurface, PETSC_DECIDE, PETSC_DECIDE, _Nmailles , _Nmailles );
		MatSetFromOptions(InvSurface);
		MatSetUp(InvSurface);
		MatZeroEntries(InvSurface);

		// Assembly of matrices 
		for (int j=0; j<_Nfaces;j++){
			Face Fj = _mesh.getFace(j);
			bool _isBoundary=Fj.isBorder();
			std::vector< int > idCells = Fj.getCellsId();

			std::map<int,int>::iterator it;
			bool periodicFaceComputed, periodicFaceNotComputed;
			if (_indexFacePeriodicSet == true ){
				it=_indexFacePeriodicMap.find(j)  ; // if periodic
				periodicFaceComputed = (it->first == j ) ;
				periodicFaceNotComputed = (it->second == j ) ;
			}
			else{
				periodicFaceComputed = false;
				periodicFaceNotComputed = false;
			}
			
			PetscScalar det, FaceArea, InvD_sigma, InvPerimeter1, InvPerimeter2;
			PetscInt IndexFace = _Nmailles + j;
			if (Fj.getNumberOfCells()==2 || (periodicFaceComputed == true) ){	// Fj is inside the domain or is a boundary periodic face (computed)
				Cell Ctemp1 ,Ctemp2 ;
				Ctemp1 = _mesh.getCell(idCells[0]);
				if ( periodicFaceComputed == true){ 
					std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
					idCells.push_back( idCells_other_Fj[0]  );
					Ctemp2 = _mesh.getCell( idCells[1]);
				}
				else
					Ctemp2 = _mesh.getCell(idCells[1]);
				FaceArea = Fj.getMeasure();
				if (_Ndim == 1){
					det = Ctemp2.x() - Ctemp1.x();
					InvPerimeter1 = 1.0/Ctemp1.getNumberOfFaces();
					InvPerimeter2 = 1.0/Ctemp2.getNumberOfFaces();
				} 
				if (_Ndim ==2){
					std::vector<int> nodes =  Fj.getNodesId();
					Node vertex = _mesh.getNode( nodes[0] );
					// determinant of the vectors forming the diamond cell around the face sigma
					det = (Ctemp1.x() - vertex.x() )* (Ctemp2.y() - vertex.y() ) - (Ctemp1.y() - vertex.y() )* (Ctemp2.x() - vertex.x() );
					InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
					InvPerimeter2 = 1/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
				}
				if (_Ndim == 3){	
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
					*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
					_runLogFile->close();
					throw CdmathException("WaveStaggered pas dispo en 3D, arret de calcul");				
				}
		
				PetscScalar MinusFaceArea = -FaceArea;
				PetscScalar InvD_sigma = 1.0/PetscAbsReal(det);
				PetscScalar InvVol1 = 1/( Ctemp1.getMeasure()* Ctemp1.getNumberOfFaces());
				PetscScalar InvVol2 = 1/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());

				MatSetValues(B, 1, &idCells[0], 1, &j, &FaceArea, ADD_VALUES ); 
				MatSetValues(B, 1, &idCells[1], 1, &j, &MinusFaceArea, ADD_VALUES );  
				MatSetValues(Bt, 1, &j, 1, &idCells[0], &FaceArea, ADD_VALUES ); 
				MatSetValues(Bt, 1, &j, 1, &idCells[1], &MinusFaceArea, ADD_VALUES ); 

				MatSetValues(Laplacian, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea, ADD_VALUES ); 
				MatSetValues(Laplacian, 1, &idCells[0], 1, &idCells[1], &FaceArea, ADD_VALUES );  
				MatSetValues(Laplacian, 1, &idCells[1], 1, &idCells[1], &MinusFaceArea, ADD_VALUES ); 
				MatSetValues(Laplacian, 1, &idCells[1], 1, &idCells[0], &FaceArea, ADD_VALUES );  

				MatSetValues(InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
				MatSetValues(InvSurface,1, &idCells[1],1, &idCells[1], &InvPerimeter2, ADD_VALUES );
				MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1 , ADD_VALUES );
				MatSetValues(_InvVol, 1, &idCells[1],1 ,&idCells[1], &InvVol2, ADD_VALUES );
				MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	
			
						
			}
			else if (Fj.getNumberOfCells()==1 && (periodicFaceNotComputed == false) ) { //if boundary face and face index is different from periodic faces not computed 
				Cell Cint = _mesh.getCell(idCells[0]);
				FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				if (_Ndim == 1){
					PetscScalar det = Fj.x() - Cint.x();
					InvD_sigma = 2.0/Cint.getMeasure() ;
					InvPerimeter1 = 1.0/Cint.getNumberOfFaces() ;
				} 
				if (_Ndim == 2){
					std::vector< int > nodes =  Fj.getNodesId();
					Node vertex1 = _mesh.getNode( nodes[0] );
					Node vertex2 = _mesh.getNode( nodes[1] );
					PetscScalar det = (Cint.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Cint.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
					// determinant of the vectors forming the interior half diamond cell around the face sigma
					InvD_sigma = 1.0/PetscAbsReal(det);
					InvPerimeter1 = 1/Cint.getNumberOfFaces();	
					
				}
				double* vec_normal = new double[_Ndim];
				for (int l=0; l< Cint.getNumberOfFaces(); l++){
					if (j == Cint.getFacesId()[l]){
						for (int idim=0; idim < _Ndim; idim ++){
							vec_normal[idim] = Cint.getNormalVector(l,idim);
						}
					}
				}
				for ( int k=0; k < _Ndim; k++){
					if (vec_normal[k] < 0)
						FaceArea = -Fj.getMeasure();
				}
				delete[] vec_normal;
				PetscScalar InvVol1 = 1.0/(Cint.getMeasure()*Cint.getNumberOfFaces());
				
				MatSetValues(B, 1, &idCells[0], 1, &j, &FaceArea, ADD_VALUES ); 
				MatSetValues(InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES ),
				MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
				MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	

				std::map<int,double> boundaryPressure = getboundaryPressure(); 
				std::map<int,double>::iterator it = boundaryPressure.find(j);
				PetscScalar pExt =Fj.getMeasure()*boundaryPressure[it->first]; 
				VecSetValues(_BoundaryTerms, 1,&idCells[0], &pExt, ADD_VALUES ); 
				MatSetValues(Laplacian, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea, ADD_VALUES ); 
			 
			}	
		}
		MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(Bt, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Bt, MAT_FINAL_ASSEMBLY);

		MatAssemblyBegin(Laplacian,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Laplacian, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_BoundaryTerms);
		VecAssemblyEnd(_BoundaryTerms);

		MatAssemblyBegin(InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);

		Mat  GradDivTilde; 
		MatScale(Bt, -1.0);
		MatMatMatMult(Bt,InvSurface, B , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &GradDivTilde); 
		if (_verbose){
			_d = 1;
			_rho =1;
			_kappa = 1;
		}
		VecScale(_BoundaryTerms, _d*_c);
		MatScale(Laplacian, _d*_c );
		MatScale(B, -1.0/_rho);
		MatScale(Bt, -1.0*_kappa);
		MatScale(GradDivTilde, _d*_c);
		if (_verbose){
			cout<<"Vectueur BoundaryTerms "<<endl;
			VecView(_BoundaryTerms,PETSC_VIEWER_STDOUT_SELF);
			cout<<"Laplacian"<<endl;
			MatView(Laplacian,PETSC_VIEWER_STDOUT_SELF);
			MatView(B,PETSC_VIEWER_STDOUT_SELF);
			MatView(Bt,PETSC_VIEWER_STDOUT_SELF);
			MatView(GradDivTilde,PETSC_VIEWER_STDOUT_SELF);
		}
		
		
		// _A = (dc Laplacian  ;  -1/rho B         )
		//      (kappa B^t     ;  dc -B^t(1/|dK|) B ) 
		Mat G[4];
		G[0] = Laplacian;
		G[1] = B;
		G[2] = Bt;
		G[3] = GradDivTilde;
		MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
		Mat Prod;
		MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
		MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, SAME_NONZERO_PATTERN); 

		if (_cfl > _d/2.0 && _Ndim > 1){
			cout << "cfl = "<< _cfl <<" is to high, cfl is updated to _d/2 = "<< 0.99*_d/2 << endl; 
		 	_cfl =  0.99 * _d/2.0; //WARNING : cfl = _d/2.0 theoretical but proof leads to think that it is the double (cfl = _d)
		}
		if (_Ndim == 1){
			cout << "the explicit in 1D is stable with cfl = 0.4, cfl is updated "<< endl;
		 	_cfl =  0.4;
		}
		Vec V, W;
		PetscScalar minInvSurf, maxInvVol;

		// Minimum size of mesh volumes
		VecCreate(PETSC_COMM_SELF, & V);
		VecSetSizes(V, PETSC_DECIDE, _globalNbUnknowns);
		int *indices3 = new int[_globalNbUnknowns];
		std::iota(indices3, indices3 +_globalNbUnknowns, 0);
		VecSetFromOptions(V);
		MatGetDiagonal(_InvVol,V);
		VecMax(V, indices3, &maxInvVol);
		_minCell = 1.0/maxInvVol;

		//Maximum size of surfaces
		VecCreate(PETSC_COMM_SELF, & W);
		VecSetSizes(W, PETSC_DECIDE, _Nmailles);
		VecSetFromOptions(W);
		MatGetDiagonal(InvSurface, W);
		int *indices4 = new int[_Nmailles];
		std::iota(indices4, indices4 +_Nmailles, 0);
		VecMin(W, indices4, &minInvSurf);
		_maxPerim = 1.0/minInvSurf;
	
		delete[] indices3, indices4;
		VecDestroy(& V);
		VecDestroy(& W); 
		MatDestroy(& B); 
		MatDestroy(& Bt); 
		MatDestroy(& InvSurface);
		MatDestroy(& Laplacian);
		MatDestroy(& GradDivTilde); 
	}
	if (_timeScheme == Explicit){	
		Vec Prod2;
		VecDuplicate(_BoundaryTerms, &Prod2);
		MatMult(_InvVol, _BoundaryTerms, Prod2);  
		MatMult(_A,_primitiveVars, _b); 
		VecAXPY(_b,     1, Prod2);
	}
	return _cfl * _minCell / (_maxPerim * _c);

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

	//Change the relaxation coefficient to ease convergence
	double relaxation=1;
	VecAXPY(_primitiveVars, relaxation, _newtonVariation);//Vk+1=Vk+relaxation*deltaV

	return true;
}

void WaveStaggered::abortTimeStep(){
	_dt = 0;
}

void WaveStaggered::validateTimeStep()
{
	//Calcul de la variation Un+1-Un
	_erreur_rel= 0;
	double x, dx;
	for(int j=1; j<_globalNbUnknowns; j++){
		VecGetValues(_newtonVariation, 1, &j, &dx);
		VecGetValues(_primitiveVars, 1, &j, &x);
		if (fabs(x)< _precision)
		{
			if(_erreur_rel < fabs(dx))
				_erreur_rel = fabs(dx);
		}
		else if(_erreur_rel < fabs(dx/x))
			_erreur_rel = fabs(dx/x);
	}

	_isStationary =_erreur_rel <_precision;

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
		if (_timeScheme == Implicit)
			cout << "Matrice du système linéaire avant contribution delta t" << endl;
		if (_timeScheme == Explicit) {
			cout << "Matrice _A tel que _A = V^-1(dc Laplacian  ;  -1/rho B         )       "<<endl; 
			cout << "                            (kappa B^t     ;  dc -B^t(1/|dK|) B) : du second membre avant contribution delta t" << endl;
		}
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

void WaveStaggered::testConservation()
{
	double SUM, DELTA, x;
	double InvcellMeasure;
	SUM = 0;
	DELTA = 0;
	//TODO : seulement si Neumann au bord ?
	for(int j=0; j<_globalNbUnknowns; j++)
	{
		VecGetValues(_primitiveVars, 1, &j, &x);//on recupere la valeur du champ
		MatGetValues(_InvVol, 1, &j,1, &j, &InvcellMeasure);//on recupere la valeur du champ
		SUM += x/InvcellMeasure;
		VecGetValues(_newtonVariation, 1, &j, &x);//on recupere la variation du champ
		DELTA += x/InvcellMeasure;
	}
	if(fabs(SUM)>_precision)
		cout << SUM << ", variation relative: " << fabs(DELTA /SUM)  << endl;
	else
		cout << " a une somme quasi nulle,  variation absolue: " << fabs(DELTA) << endl;
	
}

bool WaveStaggered::initTimeStep(double dt){
	_dt = dt;
	return _dt>0;//No need to call MatShift as the linear system matrix is filled at each Newton iteration (unlike linear problem)
}

void WaveStaggered::setVerticalPeriodicFaces(){
    for (int j=0;j<_mesh.getNumberOfFaces() ; j++){
        Face my_face=_mesh.getFace(j);
        int iface_perio=-1;
		double x=my_face.x();
		if (my_face.getNumberOfCells() ==1 && my_face.x()>0 && my_face.x()<1){ //TODO : dim =1
			/* if(_Ndim ==1){
				for (int iface=0;iface<_mesh.getNumberOfFaces() ; iface++){
					Face face_i=_mesh.getFace(iface);
					if (face_i.getNumberOfCells() ==1){
						if(iface!=j){
							iface_perio=iface;
							break;
						}
					}
				}
			} */
			if(_Ndim==2){
				for (int iface=0;iface<_mesh.getNumberOfFaces() ; iface++){
					Face face_i=_mesh.getFace(iface);
					double xi =face_i.x();
					if (face_i.getNumberOfCells() ==1 && iface !=j && ( abs(x-xi)<1e-3) ){ //TODO : pas générique quelle condition mettre pour ne pas compter face de bord
						bool empty = (_indexFacePeriodicMap.find(iface) == _indexFacePeriodicMap.end()) ;
						if (empty == true)
							_indexFacePeriodicMap[j]=iface;
					}
				}
			}
			else
				throw CdmathException("Mesh::setPeriodicFaces: Mesh dimension should be 2");		
			
			// TODO : do not fill if already in the map
		}
	}
	std::map<int,int>::iterator it;
	for (it = _indexFacePeriodicMap.begin(); it != _indexFacePeriodicMap.end(); it++ )
		cout << "_indexFacePeriodicMap["<<it->first<<"] = "<<it->second <<endl;
	_indexFacePeriodicSet = true;
}

vector<string> WaveStaggered::getInputFieldsNames()
{
	vector<string> result(1);
	
	result[0]="NOT DEFINED";
	return result;
}
void WaveStaggered::setInputField(const string& nameField, Field& inputField )
{}

double WaveStaggered::getTimeStep()
{
	return _dt;
}

void WaveStaggered::terminate(){ 
	VecDestroy(&_newtonVariation);
	VecDestroy(&_b);
	VecDestroy(&_primitiveVars);
	MatDestroy(& _A); 
	MatDestroy(&_InvVol); 
	VecDestroy(& _BoundaryTerms);	

 
	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);
}

void WaveStaggered::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/WaveStaggered_");///Results
	prim+=_fileName;

	if(_savePressure){
		for (int i = 0 ; i < _Nmailles  ; i++){
				VecGetValues(_primitiveVars,1,&i,&_Pressure(i));
			}
			
		_Pressure.setTime(_time,_nbTimeStep);
		if (_nbTimeStep ==0){
			_Pressure.setInfoOnComponent(0,"_Pressure (N/m²)");
			switch(_saveFormat)
			{
			case VTK :
				_Pressure.writeVTK(prim+"_Pressure");
				break;
			case MED :
				_Pressure.writeMED(prim+"_Pressure");
				break;
			case CSV :
				_Pressure.writeCSV(prim+"_Pressure");
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Pressure.writeVTK(prim+"_Pressure",false);
				break;
			case MED :
				_Pressure.writeMED(prim+"_Pressure",false);
				break;
			case CSV :
				_Pressure.writeCSV(prim+"_Pressure");
				break;
			}
		}
	}
	if(_saveVelocity){
		for (int i = 0 ; i < _Nfaces ; i++){
				std::map<int,int>::iterator it = _indexFacePeriodicMap.find(i); // TODO ; ne va pas fonctionner car cherche seulement dans les clés
				int k = ( (it->second == i ) && (_indexFacePeriodicSet == true) )? it->first : i; // in periodic k stays i if it has been computed by scheme and takes the value of its matched face 
				int I= _Nmailles + k;
				VecGetValues(_primitiveVars,1,&I,&_Velocity(i));
			}
			
		_Velocity.setTime(_time,_nbTimeStep);
		if (_nbTimeStep ==0){
			_Velocity.setInfoOnComponent(0,"Velocity . n_sigma_(m/s)");
			switch(_saveFormat)
			{
			case VTK :
				_Velocity.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Velocity.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Velocity.writeCSV(prim+"_Velocity");
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Velocity.writeVTK(prim+"_Velocity",false);
				break;
			case MED :
				_Velocity.writeMED(prim+"_Velocity",false);
				break;
			case CSV :
				_Velocity.writeCSV(prim+"_Velocity");
				break;
			}
		}
	}
	if(_isStationary)
	{
		prim+="_Stat";

		if(_saveVelocity){
			switch(_saveFormat)
			{
			case VTK :
				_Velocity.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Velocity.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Velocity.writeCSV(prim+"_Velocity");
				break;
			}
		}
		if(_savePressure){
			switch(_saveFormat)
			{
			case VTK :
				_Pressure.writeVTK(prim+"_Pressure");
				break;
			case MED :
				_Pressure.writeMED(prim+"_Pressure");
				break;
			case CSV :
				_Pressure.writeCSV(prim+"_Pressure");
				break;
			}
		}
	}
}
