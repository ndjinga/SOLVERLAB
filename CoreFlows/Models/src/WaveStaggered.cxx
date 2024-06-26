/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include <numeric>

using namespace std;

WaveStaggered::WaveStaggered(int dim, double kappa, double rho, MPI_Comm comm):ProblemCoreFlows(comm){
	_Ndim=dim;
	_nVar = 2; 
	_kappa = kappa;
	_rho = rho;
	_c = sqrt(kappa/rho);
	_saveVelocity=false; 
	_savePressure=false; 
	_facesBoundinit = false;
	_indexFacePeriodicSet = false;
	_vec_normal=NULL;
	if (_Ndim == 3){	
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered pas dispo en 3D, arret de calcul");				
	}
				
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

void WaveStaggered::setWallBoundIndex(int j ){
	_indexWallBoundFaceSet.push_back(j);
	_isWall = true;
}

void WaveStaggered::DisplayVelocity(){
for (int j = 0; j < _Velocity_at_Cells.getNumberOfElements(); j++) {
	for (int i=0; i< _Velocity_at_Cells.getNumberOfComponents(); i++)
		cout << "Velocity at Cells compenent["<< i<<", elem "<<j <<"]"<< _Velocity_at_Cells(j,i) <<endl;
	}
}

void WaveStaggered::setOrientation(int j,std::vector<double> vec_normal_sigma){
	for (int idim = 0; idim < _Ndim; ++idim)
		_vec_sigma[j].push_back(vec_normal_sigma[idim]);
}

double WaveStaggered::getOrientation(int j, Cell Cint){
	std::map<int, std::vector<double>  >::iterator it = _vec_sigma.find(j);
	double *vec =new double [_Ndim];
			
	for(int l=0; l<Cint.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
		if (j == Cint.getFacesId()[l]){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(l,idim);
			}
		}
	double dotprod = 0;
	double orien;
	for (int idim = 0; idim < _Ndim; ++idim)
		dotprod += vec[idim] * it->second[idim]; 

	if (dotprod > 0)
		orien = 1;
	else if (dotprod < 0)
		orien = -1;
	delete []vec;
	return orien;
}

void WaveStaggered::setExactVelocityInterpolate(const Field &Interpolate){
	_ExactVelocityInftyInterpolate = Interpolate;

	_ExactVelocityInftyInterpolate.setName("_ExactVelocityInftyInterpolate");
	_time=_ExactVelocityInftyInterpolate.getTime();
	_mesh=_ExactVelocityInftyInterpolate.getMesh();
	_ExactVelocityInftyInterpolate.setInfoOnComponent(0,"_ExactVelocityInftyInterpolate_x(m/s)");
	_ExactVelocityInftyInterpolate.setInfoOnComponent(1,"_ExactVelocityInftyInterpolate_y(m/s)");
	string prim(_path+"/WaveStaggered_");///Results
	prim+=_fileName;
	switch(_saveFormat)
	{
	case VTK :
		_ExactVelocityInftyInterpolate.writeVTK(prim+"_ExactVelocityInftyInterpolate");
		break;
	}
}


void WaveStaggered::setExactVelocityFieldAtCells(const Field &atCells){

	_ExactVelocityInftyAtCells = atCells;

	_ExactVelocityInftyAtCells.setName("_ExactVelocityInftyAtCells");
	_time=_ExactVelocityInftyAtCells.getTime();
	_mesh=_ExactVelocityInftyAtCells.getMesh();
	_ExactVelocityInftyAtCells.setInfoOnComponent(0,"ExactVelocityInfty_x(m/s)");
	_ExactVelocityInftyAtCells.setInfoOnComponent(1,"ExactVelocityInfty_y(m/s)");
	string prim(_path+"/WaveStaggered_");///Results
	prim+=_fileName;
	switch(_saveFormat)
	{
	case VTK :
		_ExactVelocityInftyAtCells.writeVTK(prim+"_ExactVelocityInftyAtCells");
		break;
	case MED :
		_ExactVelocityInftyAtCells.writeMED(prim+"_ExactVelocityInftyAtCells");
		break;
	case CSV :
		_ExactVelocityInftyAtCells.writeCSV(prim+"_ExactVelocityInftyAtCells	");
		break;
	}

}

double WaveStaggered::ErrorL2VelocityInfty(const Field &ExactVelocityInfty){
	double error =0;
	for (int j=0; j < _Nfaces; j++){
		PetscInt I = _Nmailles + j;
		double InvD_sigma;
		MatGetValues(_InvVol, 1, &I,1, &I, &InvD_sigma);
		double Dsigma = 1/InvD_sigma;
		error += Dsigma * (_Velocity(j) - ExactVelocityInfty(j))*(_Velocity(j) - ExactVelocityInfty(j));
	}
	return error;
}

void WaveStaggered::ErrorRelativeVelocityInfty(const Field &ExactVelocityInfty){
	double max = 0.1;
	for (int j=0; j < _Nfaces; j++){
		Face Fj = _mesh.getFace(j);
		double error =0;
		if ( abs(_Velocity(j) - ExactVelocityInfty(j)) > 1e-10)
			error = abs(_Velocity(j) - ExactVelocityInfty(j))/abs(ExactVelocityInfty(j));
		else 
			error = abs(_Velocity(j) - ExactVelocityInfty(j));
		if (max < error)
			max = error;
	}
	cout << "max = " << max << endl;
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
	_vec_normal = new double[_Ndim];

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

	// matrice DIVERGENCE (|K|div(u))
	MatCreate(PETSC_COMM_SELF, & _B); 
	MatSetSizes(_B, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_B);
	MatSetUp(_B);
	MatZeroEntries(_B);
	
	// matrix GRADIENT (we will impose to _Be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _Bt); 
	MatSetSizes(_Bt, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_Bt);
	MatSetUp(_Bt);
	MatZeroEntries(_Bt);

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
	if (_timeScheme == Explicit ){ 
		if ( _nbTimeStep == 0 ){
			cout << "WaveStaggered::computeTimeStep : Début calcul matrice implicite et second membre"<<endl;
			cout << endl;
			Mat Laplacian, InvSurface;
			
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
				Cell Ctemp1 = _mesh.getCell(idCells[0]);
				double orien = getOrientation(j,Ctemp1);
				

				// Metrics
				PetscScalar orientedFaceArea = orien * Fj.getMeasure();
				PetscScalar orientedMinusFaceArea = -orientedFaceArea;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				PetscScalar det, InvD_sigma, InvPerimeter1, InvPerimeter2;
				PetscInt IndexFace = _Nmailles + j;
				PetscScalar InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());

				//Is the face periodic face ? If yes will it be seen by the scheme or is it the "other" face ?
				std::map<int,int>::iterator it;
				bool periodicFaceComputed, periodicFaceNotComputed;
				if (_indexFacePeriodicSet == true ){ // if periodic 
					it = _indexFacePeriodicMap.find(j);
					periodicFaceComputed = (it != _indexFacePeriodicMap.end());
					std::map<int,int>::iterator it2 = _indexFacePeriodicMap.begin();
					while ( ( j !=it2->second) && (it2 !=_indexFacePeriodicMap.end() ) )
						it2++;
					periodicFaceNotComputed = (it2 !=  _indexFacePeriodicMap.end());
				}
				else{
					periodicFaceComputed = false;
					periodicFaceNotComputed = false;
				}			
				
				if (Fj.getNumberOfCells()==2 || (periodicFaceComputed == true) ){	// Fj is inside the domain or is a boundary periodic face (computed)
					if ( periodicFaceComputed == true){ 
						std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
						idCells.push_back( idCells_other_Fj[0]  );
					}
					Cell Ctemp2 = _mesh.getCell(idCells[1]);
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
				
					InvD_sigma = 1.0/PetscAbsReal(det);
					PetscScalar InvVol2 = 1/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());

					MatSetValues(_B, 1, &idCells[0], 1, &j, &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_B, 1, &idCells[1], 1, &j, &orientedMinusFaceArea, ADD_VALUES );  
					MatSetValues(_Bt, 1, &j, 1, &idCells[0], &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_Bt, 1, &j, 1, &idCells[1], &orientedMinusFaceArea, ADD_VALUES ); 

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
					if (_Ndim == 1){
						InvD_sigma = 2.0/Ctemp1.getMeasure() ;
						InvPerimeter1 = 1/Ctemp1.getNumberOfFaces();
					} 
					if (_Ndim == 2){
						std::vector< int > nodes =  Fj.getNodesId();
						Node vertex1 = _mesh.getNode( nodes[0] );
						Node vertex2 = _mesh.getNode( nodes[1] );
						det = (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
						// determinant of the vectors forming the interior half diamond cell around the face sigma
						InvD_sigma = 1.0/PetscAbsReal(det);	
						InvPerimeter1 = 1/Ctemp1.getNumberOfFaces(); //TODO ?? pourquoi pas pareil que face intérieure ?InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
					}
					MatSetValues(_B, 1, &idCells[0], 1, &j, &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 
					MatSetValues(Laplacian, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea, ADD_VALUES );

					//Is the face a wall boundarycondition face
					PetscScalar pExt, pInt;
					if (std::find(_indexWallBoundFaceSet.begin(), _indexWallBoundFaceSet.end(), j)!=_indexWallBoundFaceSet.end()){
						VecGetValues(_primitiveVars,1,&idCells[0],&pInt);
						pExt =  Fj.getMeasure()*pInt; //pExt = pin so (grad p)_j = 0
					}
					else{ //Imposed boundaryconditions
						std::map<int,double> boundaryPressure = getboundaryPressure(); 
						std::map<int,double>::iterator it = boundaryPressure.find(j);
						pExt = Fj.getMeasure()*boundaryPressure[it->first]; 
					}
					
					VecSetValues(_BoundaryTerms, 1,&idCells[0], &pExt, INSERT_VALUES );
				}	
			}
			MatAssemblyBegin(_B,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_B, MAT_FINAL_ASSEMBLY);
			MatAssemblyBegin(_Bt, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_Bt, MAT_FINAL_ASSEMBLY);

			MatAssemblyBegin(Laplacian,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(Laplacian, MAT_FINAL_ASSEMBLY);
			VecAssemblyBegin(_BoundaryTerms);
			VecAssemblyEnd(_BoundaryTerms);
			VecScale(_BoundaryTerms, _d * _c);

			MatAssemblyBegin(InvSurface, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(InvSurface, MAT_FINAL_ASSEMBLY);
			MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);

			Mat  GradDivTilde; 
			MatScale(_Bt, -1.0);
			MatMatMatMult(_Bt,InvSurface, _B , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &GradDivTilde); 
			MatScale(Laplacian, _d*_c );
			MatScale(_B, -1.0/_rho);
			MatScale(_Bt, -1.0*_kappa);
			MatScale(GradDivTilde, _d*_c);
			
			
			// _A = (dc Laplacian  ;  -1/rho B         )
			//      (kappa B^t     ;  dc -B^t(1/|dK|) B ) 
			Mat G[4];
			G[0] = Laplacian;
			G[1] = _B;
			G[2] = _Bt;
			G[3] = GradDivTilde;
			MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
			Mat Prod;
			MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
			MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
			MatCopy(Prod,_A, SAME_NONZERO_PATTERN); 

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
			MatDestroy(& InvSurface);
			MatDestroy(& Laplacian);
			MatDestroy(& GradDivTilde); 
		}
		if (_isWall && _nbTimeStep >0){	
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				if (Fj.getNumberOfCells()==1) { //if boundary face 
					//Is the face a wall boundarycondition face
					PetscScalar pExt, pInt;
					if (std::find(_indexWallBoundFaceSet.begin(), _indexWallBoundFaceSet.end(), j)!=_indexWallBoundFaceSet.end()){
						std::vector< int > idCells = Fj.getCellsId();
						VecGetValues(_primitiveVars,1,&idCells[0],&pInt);
						pExt = _d * _c * Fj.getMeasure()*pInt; //pExt = pin so (grad p)_j = 0
						VecSetValues(_BoundaryTerms, 1,&idCells[0], &pExt, INSERT_VALUES );
					} 
				}	
			}
			VecAssemblyBegin(_BoundaryTerms);
			VecAssemblyEnd(_BoundaryTerms);
		}
		Vec Prod2;
		VecDuplicate(_BoundaryTerms, &Prod2);
		MatMult(_InvVol, _BoundaryTerms, Prod2);  
		MatMult(_A,_primitiveVars, _b); 
		VecAXPY(_b,     1, Prod2);
	}

	ComputeEnergyAtTimeT();
	return _cfl * _minCell / (_maxPerim * _c);
}

void WaveStaggered::ComputeEnergyAtTimeT(){
	double E = 0;
	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		PetscInt I = _Nmailles + j;
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		PetscScalar InvD_sigma, InvCell1measure, InvCell2measure, pressure_in, pressure_out, velocity;
		
		if (Fj.getNumberOfCells()==2  ){	// Fj is inside the domain or is a boundary periodic face (computed)
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			MatGetValues(_InvVol, 1, &I,1, &I, &InvD_sigma);
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell1measure );
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell2measure );
			VecGetValues(_primitiveVars, 1, &idCells[0], &pressure_in );
			VecGetValues(_primitiveVars, 1, &idCells[1], &pressure_out );
			VecGetValues(_primitiveVars, 1, &I, &velocity );

			double pressure_int=  1/(InvCell1measure*Ctemp1.getNumberOfFaces()) * (pressure_in)*(pressure_in) ;
			double pressure_ext=  1/(InvCell2measure*Ctemp2.getNumberOfFaces()) * (pressure_out)*(pressure_out);
			double velocity_part = 1/(InvD_sigma) * (velocity)*(velocity);
			E += pressure_int + pressure_ext + velocity_part ;
						
		}
		else if (Fj.getNumberOfCells()==1 ) { //if boundary face and face index is different from periodic faces not computed 	
			MatGetValues(_InvVol,1,&I, 1, &I,&InvD_sigma);
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell1measure );
			VecGetValues(_primitiveVars, 1, &idCells[0], &pressure_in );
			VecGetValues(_primitiveVars, 1, &I, &velocity );

			double pressure_part_cellint=  1/(InvCell1measure*Ctemp1.getNumberOfFaces()) * (pressure_in)*(pressure_in) ; 
			double velocity_part = 1/(InvD_sigma) * (velocity)*(velocity);
			E += pressure_part_cellint + velocity_part ;
		}	
	}
	_Energy.push_back(E);
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

void WaveStaggered::validateTimeStep()
{
	//Calcul de la variation Un+1-Un
	_erreur_rel= 0;
	double x, dx;
	for(int j=0; j<_globalNbUnknowns; j++){
		VecGetValues(_newtonVariation, 1, &j, &dx);
		VecGetValues(_primitiveVars, 1, &j, &x);
		if (fabs(x)< _precision)
		{
			if(_erreur_rel < fabs(dx)){
				_erreur_rel = fabs(dx);
			}
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


bool WaveStaggered::initTimeStep(double dt){
	_dt = dt;
	return _dt>0;//No need to call MatShift as the linear system matrix is filled at each Newton iteration (unlike linear problem)
}

void WaveStaggered::abortTimeStep(){
	_dt = 0;
}

void WaveStaggered::setVerticalPeriodicFaces(){
    for (int j=0;j<_mesh.getNumberOfFaces() ; j++){
        Face my_face=_mesh.getFace(j);
        int iface_perio=-1;
		double x=my_face.x();
		if (my_face.getNumberOfCells() ==1 && my_face.x()>0 && my_face.x()<1){ //TODO : dim =1
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
		}
	}
	_indexFacePeriodicSet = true;
}

void WaveStaggered::setHorizontalPeriodicFaces(){
    for (int j=0;j<_mesh.getNumberOfFaces() ; j++){
        Face my_face=_mesh.getFace(j);
        int iface_perio=-1;
		double y=my_face.y();
		if (my_face.getNumberOfCells() ==1 && my_face.y()>0 && my_face.y()<1){ //TODO : dim =1 & : pas générique ; quelle condition mettre pour ne pas compter face de bord
			if(_Ndim==2){
				for (int iface=0;iface<_mesh.getNumberOfFaces() ; iface++){
					Face face_i=_mesh.getFace(iface);
					double yi =face_i.y();
					if (face_i.getNumberOfCells() ==1 && iface !=j && ( abs(y-yi)<1e-3) ){ 
						bool empty = (_indexFacePeriodicMap.find(iface) == _indexFacePeriodicMap.end()) ;
						if (empty == true)
							_indexFacePeriodicMap[j]=iface;
					}
				}
			}
			else
				throw CdmathException("Mesh::setPeriodicFaces: Mesh dimension should be 2");		
		}
	}
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
	delete[]_vec_normal;
	VecDestroy(&_newtonVariation);
	VecDestroy(&_b);
	VecDestroy(&_primitiveVars);
	MatDestroy(& _A); 
	MatDestroy(& _B); 
	MatDestroy(& _Bt); 
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
	if(_saveVelocity  ){ 
		Field _Velocity_at_Cells("Velocity at cells results", CELLS, _mesh,3);
		Field  _DivVelocity("velocity divergence", CELLS, _mesh, 1);

		_Velocity_at_Cells.setTime(_time,_nbTimeStep);
		_DivVelocity.setTime(_time,_nbTimeStep);
		for (int l=0; l < _Nmailles ; l++){
			_DivVelocity(l) =0;
			for (int k=0; k< 3; k++){
				_Velocity_at_Cells(l, k) =0;
			}
		}
		
		for (int i = 0 ; i < _Nfaces ; i++){
			bool periodicFaceNotComputed;
			std::map<int,int>::iterator it2 = _indexFacePeriodicMap.begin();
			while ( ( i !=it2->second) && (it2 != _indexFacePeriodicMap.end() ) )
				it2++;
			periodicFaceNotComputed = (it2 !=  _indexFacePeriodicMap.end());
			int k = (periodicFaceNotComputed ==true) && (_indexFacePeriodicSet == true) ? it2->first : i; // in periodic k stays i, if it has been computed by scheme and takes the value of its matched face 
			int I= _Nmailles + k;
			VecGetValues(_primitiveVars,1,&I,&_Velocity(i));

			Face Fj = _mesh.getFace(i);
			std::vector< int > idCells = Fj.getCellsId();
			Cell Ctemp1 = _mesh.getCell(idCells[0]); //origin of the normal vector
		
			if (_Ndim >1 ){
				bool found = false;
				for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
					if (i == Ctemp1.getFacesId()[l]){
						found = true;
						for (int idim = 0; idim < _Ndim; ++idim)
							_vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
					}
				}
				assert(found);
			}

			double orien1 = getOrientation(i,Ctemp1);
			PetscScalar det, detL, detR, D_sigmaL, D_sigmaR;
			/* cout << " \n 2) Ctemp1.x() =  "<< Ctemp1.x() <<" Ctemp1.y() = "<< Ctemp1.y() <<" Fj.x()==" << Fj.x() << " Fj.y() ="<<Fj.y() <<endl;
			cout << " Velocity( "<< i<< " )= " <<  _Velocity(i)  << endl;
			for (int k=0; k< 2; k++)
				cout << "normal ["<< k <<"] = " << _vec_normal[k] << endl; */

			if (Fj.getNumberOfCells() ==2){
				Cell Ctemp2 = _mesh.getCell(idCells[1]); 
				/* if (_Ndim ==2){
					std::vector<int> nodes =  Fj.getNodesId();
					Node vertex = _mesh.getNode( nodes[0] );
					Node vertex2 = _mesh.getNode( nodes[1]);
					detL = (vertex.x()-vertex2.x())*(Ctemp1.y() - vertex2.y() ) - (vertex.y()-vertex2.y())*(Ctemp1.x() - vertex.x() ); 
					detR = (vertex.x()-vertex2.x())*(Ctemp2.y() - vertex2.y() ) - (vertex.y()-vertex2.y())*(Ctemp2.x() - vertex.x() ) ;
					D_sigmaL= PetscAbsReal(detL)/2.0;
					D_sigmaR= PetscAbsReal(detR)/2.0;
				} */
				
				for (int k=0; k< 2; k++){ //TODO : cas _ndim = 1 !
					_Velocity_at_Cells(idCells[0], k) +=  _Velocity(i) * _vec_normal[k]/4.0; 
					_Velocity_at_Cells(idCells[1], k) +=  _Velocity(i) * _vec_normal[k]/4.0; 
				}
				_DivVelocity( idCells[0]) += Fj.getMeasure() * orien1 * _Velocity(i)/(Ctemp1.getMeasure());
				_DivVelocity( idCells[1]) -= Fj.getMeasure() * orien1 * _Velocity(i)/(Ctemp2.getMeasure());
			}
			else if (Fj.getNumberOfCells() ==1){
				/* if (_Ndim ==2){	
					std::vector<int> nodes =  Fj.getNodesId();
					Node vertex = _mesh.getNode( nodes[0] );
					Node vertex2 = _mesh.getNode( nodes[1]);
					detL = (vertex.x()-vertex2.x())*(Ctemp1.y() - vertex2.y() ) - (vertex.y()-vertex2.y())*(Ctemp1.x() - vertex.x() );
					D_sigmaL= PetscAbsReal(detL)/2.0;
				} */
				for (int k=0; k< 2; k++){ //TODO : cas _ndim = 1 !
					_Velocity_at_Cells(idCells[0], k) +=  _Velocity(i) * _vec_normal[k]/4.0; 
					} 
				_DivVelocity( idCells[0]) += Fj.getMeasure() * orien1 * _Velocity(i)/(Ctemp1.getMeasure());
			}
		}

		_Velocity.setTime(_time,_nbTimeStep);
		_Velocity_at_Cells.setTime(_time,_nbTimeStep);
		_DivVelocity.setTime(_time,_nbTimeStep);
		_Velocity.setInfoOnComponent(0,"Velocity . n_sigma_(m/s)");
		_Velocity_at_Cells.setInfoOnComponent(0,"Velocity at cells x_(m/s)");
		_Velocity_at_Cells.setInfoOnComponent(1,"Velocity at cells y_(m/s)");
		_DivVelocity.setInfoOnComponent(0,"divergence velocity (s^-1)");
	
		switch(_saveFormat)
		{
		case VTK :
			_Velocity_at_Cells.writeVTK(prim+"_Velocity at cells");
			_DivVelocity.writeVTK(prim+"Divergence Velocity");
			_Velocity.writeVTK(prim+"_Velocity");
			break;
		case MED :
			_Velocity.writeMED(prim+"_Velocity");
			break;
		case CSV :
			_Velocity.writeCSV(prim+"_Velocity");
			break;
		}

		if (_isStationary || _time == _timeMax){
			double boundaryIntegral =0;
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				if (Fj.getNumberOfCells() == 1){ 
					std::vector< int > idCells = Fj.getCellsId();
					Cell Ctemp1 = _mesh.getCell(idCells[0]);

					double u;
					int I = _Nmailles + j;
					VecGetValues(_primitiveVars, 1, &I, &u);
					double orien1 = getOrientation(j, Ctemp1);
					boundaryIntegral += Fj.getMeasure() * orien1 * u;
				}
			}
			double norm = 0;
			for (int i = 0; i < _Nmailles; i++){
				if (norm < fabs(_DivVelocity(i)))
					norm = fabs(_DivVelocity(i));	
			}
			cout << "max|div(u)|= "<< norm << " while /int_{/partial /Omega} u_b.n d/gamma = "<< boundaryIntegral <<endl;
		}
	}
}
