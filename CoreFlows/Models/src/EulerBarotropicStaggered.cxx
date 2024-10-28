/*
 * EulerBarotropicStaggered.cxx
 */

#include "EulerBarotropicStaggered.hxx"
#include "StiffenedGas.hxx"
#include <cassert>
#include <numeric>

using namespace std;

EulerBarotropicStaggered::EulerBarotropicStaggered(phaseTypeStaggered fluid, pressureEstimate pEstimate, int dim):WaveStaggered( dim, 1.0 ,1.0, MPI_COMM_WORLD ){
;	_Ndim=dim;
	_nVar = 2; 
	_saveVelocity=false; 
	_savePressure=false; 
	_facesBoundinit = false;
	_indexFacePeriodicSet = false;
	_vec_normal=NULL;
	_Time.push_back(0);

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		_Tref=300;
		_Pref=1e5;
		if(fluid==GasStaggered){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_compressibleFluid = new StiffenedGas(1.4,743,_Tref,2.22e5);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, e=2.22e5, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_compressibleFluid = new StiffenedGas(996,_Pref,_Tref,1.12e5,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C, e=1.12e5, c_v=4130
		}
	}
	else{//EOS at 155 bars and 618K 
		_Tref=618;
		_Pref=155e5;
		if(fluid==GasStaggered){
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_compressibleFluid = new StiffenedGas(102,_Pref,_Tref,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
		}
		else{
			cout<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			_compressibleFluid= new StiffenedGas(726.82,_Pref,_Tref,1.3e6, 971.,5454.);  //stiffened gas law for water at pressure 155 bar, and temperature 345°C
		}
	}
	_c = _compressibleFluid->vitesseSonEnthalpie(_compressibleFluid->getEnthalpy(_Tref, _compressibleFluid->getDensity(_Pref, _Tref)));
	//TODO découplage burgers _c=0
	//Save into the fluid list
	_fluides.resize(0);
	_fluides.push_back(_compressibleFluid);
	if (_Ndim == 3){	
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!EulerBarotropicStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!EulerBarotropicStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("EulerBarotropicStaggered pas dispo en 3D, arret de calcul");				
	}		
}


void EulerBarotropicStaggered::initialize(){
	double * initialFieldVelocity = new double[_Nfaces];
	double * initialFieldPressure = new double[_Nmailles];
	if (_mpi_rank == 0){
		cout<<"\n Initialising the Wave System model\n"<<endl;
		*_runLogFile<<"\n Initialising the Wave Sytem model\n"<<endl;

		_globalNbUnknowns = _Nmailles + _Nfaces; //Staggered discretisation : velocity is on faces

		if(!_initialDataSet)
		{
			*_runLogFile<<"!!!!!!!!EulerBarotropicStaggered::initialize() set initial data first"<<endl;
			_runLogFile->close();
			throw CdmathException("!!!!!!!!EulerBarotropicStaggered::initialize() set initial data first");
		}
		cout << "mesh dimension = "<<_Ndim <<endl;
		*_runLogFile << " spaceDim= "<<_Ndim <<endl;

		_vec_normal = new double[_Ndim];
		//Construction des champs primitifs initiaux comme avant dans ParaFlow
		
		for(int i =0; i<_Nfaces; i++)
			initialFieldVelocity[i]=_Velocity(i); 
			
		for(int i =0; i<_Nmailles; i++)
			initialFieldPressure[i]=_Pressure(i); 
	}

	/**********Petsc structures:  ****************/
	
	//creation des vecteurs
	VecCreate(PETSC_COMM_SELF, & _primitiveVars);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_primitiveVars, PETSC_DECIDE, _globalNbUnknowns);
	VecSetFromOptions(_primitiveVars);
	VecDuplicate(_primitiveVars, &_newtonVariation);//Newton variation Uk+1-Uk 
	VecDuplicate(_primitiveVars, &_b);//Right hand side of Newton method

	// transfer information de condition initial vers primitiveVars  
	if (_mpi_rank == 0){
		int *indices1 = new int[_Nmailles];
		int *indices2 = new int[_Nfaces];
		std::iota(indices1, indices1 + _Nmailles, 0);
		std::iota(indices2, indices2 + _Nfaces, _Nmailles);
		VecSetValues(_primitiveVars, _Nmailles , indices1, initialFieldPressure, INSERT_VALUES); 
		VecSetValues(_primitiveVars, _Nfaces, indices2, initialFieldVelocity, INSERT_VALUES);
		delete[] initialFieldVelocity;
		delete[] initialFieldPressure;
		delete[] indices1;
		delete[] indices2;
	} 
	VecAssemblyBegin(_primitiveVars);
	VecAssemblyEnd(_primitiveVars);

	VecCreate(PETSC_COMM_SELF, & _DualDensity);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_DualDensity, PETSC_DECIDE, _Nfaces);
	VecSetFromOptions(_DualDensity);
	VecZeroEntries(_DualDensity);


	//************ Time independent matrices ********** //
	// matrice des Inverses Volumes V^{-1}
	MatCreate(PETSC_COMM_SELF, &_InvVol); 
	MatSetSizes(_InvVol, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_InvVol);
	MatSetUp(_InvVol);
	MatZeroEntries(_InvVol);

	// matrice des Inverses Volumes Primaux V^{-1}
	MatCreate(PETSC_COMM_SELF, &_InvVolPrim); 
	MatSetSizes(_InvVolPrim, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles );
	MatSetFromOptions(_InvVolPrim);
	MatSetUp(_InvVolPrim);
	MatZeroEntries(_InvVolPrim);

	// matrice des Inverses Duaux V^{-1}
	MatCreate(PETSC_COMM_SELF, &_InvVolDual); 
	MatSetSizes(_InvVolDual, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces );
	MatSetFromOptions(_InvVolDual);
	MatSetUp(_InvVolDual);
	MatZeroEntries(_InvVolDual);

	// matrice des Inverses de Surfaces
	MatCreate(PETSC_COMM_SELF, & _InvSurface); 
	MatSetSizes(_InvSurface, PETSC_DECIDE, PETSC_DECIDE, _Nmailles , _Nmailles );
	MatSetFromOptions(_InvSurface);
	MatSetUp(_InvSurface);
	MatZeroEntries(_InvSurface);

	AssembleMetricsMatrices();
	
	// *************** Time dependent matrices ************ //
	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);
	MatZeroEntries(_A);

	// matrice DIVERGENCE (|K|div(u))
	MatCreate(PETSC_COMM_SELF, & _Div); 
	MatSetSizes(_Div, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_Div);
	MatSetUp(_Div);

	// matrice DIVERGENCE  rho U (|K|div(rho u))
	MatCreate(PETSC_COMM_SELF, & _DivRhoU); 
	MatSetSizes(_DivRhoU, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_DivRhoU);
	MatSetUp(_DivRhoU);

	// matrix LAPLACIAN Pressure (without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianPressure); 
	MatSetSizes(_LaplacianPressure, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles ); 
	MatSetFromOptions(_LaplacianPressure);
	MatSetUp(_LaplacianPressure);
	
	// matrix GRADIENT (we will impose to be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_DivTranspose, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_DivTranspose);
	MatSetUp(_DivTranspose);

	// matrix CONVECTION 
	MatCreate(PETSC_COMM_SELF, & _Conv); 
	MatSetSizes(_Conv, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces );
	MatSetFromOptions(_Conv);
	MatSetUp(_Conv);

	// matrix LAPLACIAN VeLOCITY(without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianVelocity); 
	MatSetSizes(_LaplacianVelocity, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces ); 
	MatSetFromOptions(_LaplacianVelocity);
	MatSetUp(_LaplacianVelocity);

	// Vector BoundaryTerms for Pressure
	VecCreate(PETSC_COMM_SELF, & _BoundaryTerms); 
	VecSetSizes(_BoundaryTerms, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_BoundaryTerms);
	VecSetUp(_BoundaryTerms);
	VecZeroEntries(_BoundaryTerms);
	

	_ConvectiveMax = 0;
	if(_system)
	{
		cout << "Variables primitives initiales : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}

	if(_mpi_size>1 && _mpi_rank == 0)
    	VecCreateSeq(PETSC_COMM_SELF, _globalNbUnknowns, &_primitiveVars_seq);//For saving results on proc 0
    VecScatterCreateToZero(_primitiveVars,&_scat,&_primitiveVars_seq);


	createKSP();
	PetscPrintf(PETSC_COMM_WORLD,"SOLVERLAB Newton solver ");
	*_runLogFile << "SOLVERLAB Newton solver" << endl;
	_runLogFile->close();

	_initializedMemory=true;
	UpdateDualDensity();
	save();//save initial data
}

std::vector<double>  EulerBarotropicStaggered::getTimeEvol(){
	return _Time;
}

void EulerBarotropicStaggered::AssembleMetricsMatrices(){
	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		PetscScalar det, InvD_sigma, InvPerimeter1, InvPerimeter2;	;
		PetscScalar InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;		

		if (IsInterior){
			std::map<int,int>::iterator it = _FacePeriodicMap.find(j);
			if ( it != _FacePeriodicMap.end()  ){ 
				std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
				idCells.push_back( idCells_other_Fj[0]  );
			}
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			PetscScalar InvVol2 = 1.0/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());
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
				InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
				InvPerimeter2 = 1.0/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
			}
			InvD_sigma = 1.0/PetscAbsReal(det);
			MatSetValues(_InvVolPrim, 1, &idCells[0],1 ,&idCells[0], &InvVol1 , ADD_VALUES );
			MatSetValues(_InvVolPrim, 1, &idCells[1],1 ,&idCells[1], &InvVol2, ADD_VALUES );
			MatSetValues(_InvVolDual, 1, &j, 1, &j,  &InvD_sigma, ADD_VALUES); 	
			MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
			MatSetValues(_InvSurface,1, &idCells[1],1, &idCells[1], &InvPerimeter2, ADD_VALUES );
		}
		else if (IsWallBound || IsSteggerBound ) { 
			if (_Ndim == 1){
				InvD_sigma = 2.0/Ctemp1.getMeasure() ;
				InvPerimeter1 = 1.0/Ctemp1.getNumberOfFaces();
			} 
			if (_Ndim == 2){
				std::vector< int > nodes =  Fj.getNodesId();
				Node vertex1 = _mesh.getNode( nodes[0] );
				Node vertex2 = _mesh.getNode( nodes[1] );
				det = (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
				// determinant of the vectors forming the interior half diamond cell around the face sigma
				InvD_sigma = 1.0/PetscAbsReal(det);	
				InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces() ); 
			}
			MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
			MatSetValues(_InvVolPrim, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
			MatSetValues(_InvVolDual, 1, &j, 1, &j,  &InvD_sigma, ADD_VALUES); 
		}
	}
	MatAssemblyBegin(_InvVolPrim,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvVolPrim, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(_InvVolDual,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvVolDual, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(_InvSurface,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
}




void EulerBarotropicStaggered::UpdateDualDensity(){
	VecZeroEntries(_DualDensity);

	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		PetscScalar detL, detR, D_sigmaL, D_sigmaR, rho_sigma,inv_rho_sigma, rhoL, rhoR;

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;		
		PetscInt I = _Nmailles + j;
		PetscScalar one = 1.0;
		if (IsInterior){
			std::map<int,int>::iterator it = _FacePeriodicMap.find(j);
			if ( it != _FacePeriodicMap.end()  ){ 
				std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
				idCells.push_back( idCells_other_Fj[0]  );
			}
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			if (_Ndim == 1){
				detL = Ctemp2.x() - Ctemp1.x()/2.0;
				detR = Ctemp2.x() - Ctemp1.x()/2.0;
			}
			if (_Ndim ==2){
				std::vector< int > nodes =  Fj.getNodesId();
				Node vertex1 = _mesh.getNode( nodes[0] );
				Node vertex2 = _mesh.getNode( nodes[1] );
				detL = (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() )/2.0;
				detR = (Ctemp2.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp2.y() - vertex1.y() )* (vertex2.x() - vertex1.x() )/2.0;
			}
			D_sigmaL = PetscAbsReal(detL);
			D_sigmaR = PetscAbsReal(detR);
			VecGetValues(_primitiveVars,1,&idCells[0],&rhoL);
			VecGetValues(_primitiveVars,1,&idCells[1],&rhoR);
			rho_sigma = (rhoL * D_sigmaL + rhoR*D_sigmaR)/(D_sigmaL +  D_sigmaR) ;
			//rho_sigma = 1.0; //TODO découplage advection
			inv_rho_sigma = 1.0/rho_sigma;
			VecSetValues(_DualDensity, 1, &j, &rho_sigma, INSERT_VALUES );
		}
	}
	VecAssemblyBegin(_DualDensity);
	VecAssemblyEnd(_DualDensity);
}

double EulerBarotropicStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme
	/************ Max rho, Max u *******************/
	PetscScalar rho,p, u;
	PetscInt zero=0;
	VecGetValues(_primitiveVars,1,&zero,&rho);
	VecGetValues(_primitiveVars,1,&_Nmailles,&u);
	_uMax = abs(u);
	_rhoMax = rho;
	for (int n=0; n <_Nmailles; n++){
		VecGetValues(_primitiveVars,1,&n,&rho);
		_rhoMax = max(_rhoMax, rho);
	}
	for (int n=0; n <_Nfaces; n++){
		PetscInt I = _Nmailles + n;
		VecGetValues(_primitiveVars,1,&I,&u);
		_uMax = max(_uMax, abs(u));
	}

	MatZeroEntries(_DivRhoU);
	MatZeroEntries(_LaplacianPressure);
	MatZeroEntries(_InvSurface);
	MatZeroEntries(_InvVol);

	MatZeroEntries(_Conv);
	MatZeroEntries(_Div); 
	MatZeroEntries(_DivTranspose); 
	MatZeroEntries(_LaplacianVelocity);
	
	if (_timeScheme == Explicit ){ 
		if (_mpi_rank ==0){
			// Assembly of matrices 
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				std::vector< int > idCells = Fj.getCellsId();
				Cell Ctemp1 = _mesh.getCell(idCells[0]);
				PetscScalar epsilon, ConvectiveFlux, rhoL, rhoR, u, absConvectiveFlux, MinusabsConvectiveFlux,det, InvD_sigma, InvPerimeter1, InvPerimeter2;	
				
				// Metrics
				PetscScalar orientedFaceArea = getOrientation(j,Ctemp1) * Fj.getMeasure();
				PetscScalar orientedMinusFaceArea = -orientedFaceArea;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				PetscInt IndexFace = _Nmailles + j;
				PetscScalar InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());

				bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
				bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
				bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;			
		
				if (IsInterior){
					std::map<int,int>::iterator it = _FacePeriodicMap.find(j);
					if ( it != _FacePeriodicMap.end()  ){ 
						std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
						idCells.push_back( idCells_other_Fj[0]  );
					}
					Cell Ctemp2 = _mesh.getCell(idCells[1]);

					/*************** Matrices related to metrics ***************/
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
						InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
						InvPerimeter2 = 1.0/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
					}
					InvD_sigma = 1.0/PetscAbsReal(det);
					PetscScalar InvVol2 = 1.0/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());
					MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1 , ADD_VALUES );
					MatSetValues(_InvVol, 1, &idCells[1],1 ,&idCells[1], &InvVol2, ADD_VALUES );
					MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	
					MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
					MatSetValues(_InvSurface,1, &idCells[1],1, &idCells[1], &InvPerimeter2, ADD_VALUES );

					/****************** Density conservation equation *********************/
					VecGetValues(_primitiveVars,1,&idCells[0],&rhoL);
					VecGetValues(_primitiveVars,1,&idCells[1],&rhoR);
					VecGetValues(_primitiveVars,1,&IndexFace,&u);

					PetscScalar orientedFaceArea_densityMean = orientedFaceArea * (rhoL + rhoR)/2.0;
					PetscScalar MinusorientedFaceArea_densityMean = -orientedFaceArea_densityMean;
					PetscScalar FaceArea_upwinding = (abs(u) + _c ) * FaceArea/2.0; 
					PetscScalar MinusFaceArea_upwinding = -FaceArea_upwinding;
					MatSetValues(_DivRhoU, 1, &idCells[0], 1, &j, &orientedFaceArea_densityMean, ADD_VALUES ); 
					MatSetValues(_DivRhoU, 1, &idCells[1], 1, &j, &MinusorientedFaceArea_densityMean, ADD_VALUES );  
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea_upwinding, ADD_VALUES ); 
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[1], &FaceArea_upwinding, ADD_VALUES );  
					MatSetValues(_LaplacianPressure, 1, &idCells[1], 1, &idCells[1], &MinusFaceArea_upwinding, ADD_VALUES ); 
					MatSetValues(_LaplacianPressure, 1, &idCells[1], 1, &idCells[0], &FaceArea_upwinding, ADD_VALUES );  

					/*************** Momentum conservation equation *****************/
					MatSetValues(_Div, 1, &idCells[0], 1, &j, &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_Div, 1, &idCells[1], 1, &j, &orientedMinusFaceArea, ADD_VALUES ); 
					MatSetValues(_DivTranspose, 1, &j, 1, &idCells[0], &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_DivTranspose, 1, &j, 1, &idCells[1], &orientedMinusFaceArea, ADD_VALUES ); 
					
					
					// Convective terms //
					PetscInt jepsilon, L, I;
					// Loop on half diamond cells
					for (int nei =0; nei <Fj.getNumberOfCells(); nei ++){
						Cell K = _mesh.getCell(Fj.getCellsId() [nei]);
						Point xb =  K.getBarryCenter();
						std::vector<int> idFaces = K.getFacesId();
						VecGetValues(_primitiveVars,1,&Fj.getCellsId() [nei],&rhoL);
						//Loop on epsilon the boundary of a dual cell in the fixed half diamond cell 
						for (int Nbepsilon = 0; Nbepsilon < Fj.getNumberOfNodes() ;Nbepsilon ++){ //TODO : works in 1D but is not clean
							ConvectiveFlux = 0 ;		
							Node xsigma = _mesh.getNode( Fj.getNodesId()[Nbepsilon] ); //TODO DIMENSION DEPENDANT
							// For fixed epsilon (and thus the node on sigma defining it ) Loop on the faces that are in the boundary of the cell containing the fixed half diamond cell	
							// Compute the flux through epsilon
							for (int f =0; f <K.getNumberOfFaces(); f ++){
								bool IsfInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),idFaces[f] ) != _InteriorFaceSet.end() ;
								bool IsfWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),idFaces[f] ) != _WallBoundFaceSet.end() ;
								bool IsfSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),idFaces[f] ) != _SteggerBoundFaceSet.end() ;			
								Face Facef = _mesh.getFace( idFaces[f] );
								std::vector< int> idCellsOfFacef =  Facef.getCellsId();
								I = _Nmailles + idFaces[f];
								VecGetValues(_primitiveVars,1,&I	,&u);
								double psif;
								if (_Ndim==1){
									if (Facef.x() < K.getBarryCenter().x())
										psif = -1.0/K.getNumberOfFaces();        
									else      		
										psif = 1.0/K.getNumberOfFaces();
								}
								//TODO : calculer 2d
								if (IsfInterior){												
									std::map<int,int>::iterator it = _FacePeriodicMap.find(f);
									if ( it != _FacePeriodicMap.end()  ){ 
										std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
										idCellsOfFacef.push_back( idCells_other_Fj[0]  );
									}
									VecGetValues(_primitiveVars,1,&idCellsOfFacef[1],&rhoR);	
								}
								else if (IsfWallBound ){			
									rhoR =  rhoL;
								}
								else if (IsfSteggerBound){ 
									std::map<int,double> boundaryPressure = getboundaryPressure(); 
									std::map<int,double>::iterator it = boundaryPressure.find( idFaces[f]);
									rhoR = boundaryPressure[it->first]; 
								} 
								/* rhoL =1.0; //TODO Découplage advection 
								rhoR =1.0;  */
								
								ConvectiveFlux += ( u *(rhoL + rhoR)/2.0  - (abs(u) + _c)* (rhoR - rhoL)/2.0* getOrientation(idFaces[f], K) )* psif ; // TODO  
								
								std::map<int,int>::iterator it = _FacePeriodicMap.begin();
								if (_Ndim == 1 && K.getFacesId()[f] != j){ // -> Search for the unique face that is not sigma that is in the boundary of K-> jepsilon will be the index of this cell
									while (it->second != K.getFacesId()[f] &&  it != _FacePeriodicMap.end() )
										it++;
									if (it == _FacePeriodicMap.end())
										jepsilon =  K.getFacesId()[f];
									else 
										jepsilon = it->first;
								}
								if (_Ndim == 2){ // -> Search for the unique face that is not sigma that has x_sigma has an extremity -> jepsilon will be the index of this cell
									std::vector< int > idNodes = Facef.getNodesId(); 
									if ( std::find(idNodes.begin(), idNodes.end(), Fj.getNodesId()[Nbepsilon])!=idNodes.end() && K.getFacesId()[f] != j){
										while (it->second != K.getFacesId()[f] &&  it != _FacePeriodicMap.end() )
											it++;
										if (it == _FacePeriodicMap.end())
											jepsilon =  K.getFacesId()[f];
										else 
											jepsilon = it->first;
									}
								}
							}
							epsilon = 1.0;
							if (_Ndim > 1)
								epsilon = (xb.x() - xsigma.x())*(xb.x() - xsigma.x()) + (xb.y() - xsigma.y() )*(xb.y() - xsigma.y());
							epsilon = sqrt(epsilon);
							_ConvectiveMax = max(abs(ConvectiveFlux), abs(_ConvectiveMax));
							ConvectiveFlux *= epsilon/2.0  ;
							MatSetValues(_Conv, 1, &j, 1, &j, &ConvectiveFlux, ADD_VALUES );  		
							MatSetValues(_Conv, 1, &j, 1, &jepsilon, &ConvectiveFlux, ADD_VALUES ); 
							absConvectiveFlux = epsilon *_rhoMax * _uMax/2.0; //abs(ConvectiveFlux); TODO que faire des ces termes  ?
							MinusabsConvectiveFlux = -epsilon * _rhoMax * _uMax/2.0; //-abs(ConvectiveFlux);
							MatSetValues(_LaplacianVelocity, 1, &j, 1, &j, &MinusabsConvectiveFlux, ADD_VALUES ); 
							MatSetValues(_LaplacianVelocity, 1, &j, 1, &jepsilon, &absConvectiveFlux, ADD_VALUES ); 
						}
					} 
				}
				else if (IsWallBound || IsSteggerBound ) { 
					/*************** Matrices related to metrics ***************/		
					if (_Ndim == 1){
						InvD_sigma = 2.0/Ctemp1.getMeasure() ;
						InvPerimeter1 = 1.0/Ctemp1.getNumberOfFaces();
					} 
					if (_Ndim == 2){
						std::vector< int > nodes =  Fj.getNodesId();
						Node vertex1 = _mesh.getNode( nodes[0] );
						Node vertex2 = _mesh.getNode( nodes[1] );
						det = (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
						// determinant of the vectors forming the interior half diamond cell around the face sigma
						InvD_sigma = 1.0/PetscAbsReal(det);	
						InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces() ); 
					}
					MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 
					
					/****************** Density conservation equation *********************/
					PetscScalar rhoExt, rhoInt, u;
					VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);
					VecGetValues(_primitiveVars,1,&IndexFace,&u);
					
					//Is the face a wall boundarycondition face
					if (IsWallBound){	
						rhoExt =  rhoInt;
					}
					else if (IsSteggerBound){ //Imposed boundaryconditions
						std::map<int,double> boundaryPressure = getboundaryPressure(); 
						std::map<int,double>::iterator it = boundaryPressure.find(j);
						rhoExt = boundaryPressure[it->first]; 
					}
					PetscScalar orientedFaceArea_densityMean = orientedFaceArea * (rhoInt + rhoExt)/2.0;
					MatSetValues(_DivRhoU, 1, &idCells[0], 1, &j, &orientedFaceArea_densityMean, ADD_VALUES ); 
					MatSetValues(_Div, 1, &idCells[0], 1, &j, &orientedFaceArea, ADD_VALUES ); 
					
					PetscScalar MinusFaceArea_upwinding = -((abs(u) + _c ) * FaceArea)/2.0; 
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea_upwinding, ADD_VALUES );
					PetscScalar boundterm = -rhoExt*MinusFaceArea_upwinding;
					VecSetValues(_BoundaryTerms, 1,&idCells[0], &boundterm, INSERT_VALUES );
				}	
			}
		}
		MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);

		MatAssemblyBegin(_DivRhoU,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_DivRhoU, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_LaplacianPressure,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_LaplacianPressure, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_BoundaryTerms);
		VecAssemblyEnd(_BoundaryTerms);

		MatAssemblyBegin(_DivTranspose, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_DivTranspose, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_Div, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_Div, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_Conv, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_Conv, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_LaplacianVelocity, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_LaplacianVelocity, MAT_FINAL_ASSEMBLY);		

		
		
		/***********Assembling the matrix _A such that : *************/
		// _A = ( (u+c) _Laplacian   ;                 -Div(\rhomean .)                  )
		//      (0_FacesxCells   ; -Conv+ c (\tilde rho) MinusGrad 1/|dK| Div + LaplacianVelocity	)

		MatMatMatMult(_DivTranspose,_InvSurface, _Div , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &_GradDivTilde); // -grad (inv_Surf) Div													
		MatScale(_GradDivTilde, -1.0 * _c*_rhoMax) ; 										// -(-grad (inv_Surf) Div) = grad (inv_Surf) Div
		MatScale(_DivRhoU, -1.0);	
		//MatScale(_DivRhoU, 0);	//TODO découplage burgers
		MatAXPY(_GradDivTilde, 1, _LaplacianVelocity, UNKNOWN_NONZERO_PATTERN);
		MatAXPY(_GradDivTilde, -1, _Conv, UNKNOWN_NONZERO_PATTERN); 

		Mat G[4], ZeroNfaces_Ncells;
		MatCreate(PETSC_COMM_SELF, &ZeroNfaces_Ncells); 
		MatSetSizes(ZeroNfaces_Ncells, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
		MatSetFromOptions(ZeroNfaces_Ncells);
		MatSetUp(ZeroNfaces_Ncells);
		MatZeroEntries(ZeroNfaces_Ncells);
		MatAssemblyBegin(ZeroNfaces_Ncells, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(ZeroNfaces_Ncells, MAT_FINAL_ASSEMBLY);	//TODO dupliquer grad dans zeroNfaces_ncells 
		G[0] = _LaplacianPressure;
		G[1] = _DivRhoU;
		G[2] = ZeroNfaces_Ncells ;
		G[3] = _GradDivTilde;
		MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
		Mat Prod; 
		MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
		MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, UNKNOWN_NONZERO_PATTERN); 						
		MatMult(_A,_primitiveVars, _b); 
		
		/***********Adding boundary terms and pressure gradient**************/
		Vec Temporary1,Temporary2, Pressure, GradPressure;
		VecCreate(PETSC_COMM_SELF, & Temporary1);
		VecCreate(PETSC_COMM_SELF, & Temporary2);
		VecSetSizes(Temporary1, PETSC_DECIDE, _globalNbUnknowns);
		VecSetSizes(Temporary2, PETSC_DECIDE, _Nfaces);
		VecSetFromOptions(Temporary1);
		VecSetFromOptions(Temporary2);
		//add boundary terms to AU^n//
		MatMult(_InvVol, _BoundaryTerms, Temporary1); 
		VecAXPY(_b,     1, Temporary1); 
		//Extract pressure to multiply by _DivTranspose //
		VecCreate(PETSC_COMM_SELF, & Pressure); 
		VecSetSizes(Pressure, PETSC_DECIDE, _Nmailles); 
		VecSetFromOptions(Pressure);
		VecSetUp(Pressure);
		for (int i=0; i <_Nmailles; i++){
			VecGetValues(_primitiveVars,1,&i,&rho);
			p =  _compressibleFluid->getPressureFromEnthalpy(_compressibleFluid->getEnthalpy(_Tref,rho) ,rho);											
			VecSetValues(Pressure, 1,&i, &p, INSERT_VALUES );
		}														
		//add pressure gradient to AU^n + Boundterms//
		MatMult(_DivTranspose, Pressure, Temporary2); 						
		double *Product = new double[_Nfaces]; 								//TODO can we avoid creating a pointer only to insert in the bigger Vector ?
		PetscScalar gradp;
		for (int i=0; i <_Nfaces; i++){
			VecGetValues(Temporary2,1,&i,&gradp);
			Product[i] = gradp; // TODO Découplage gradp;
		}
		VecDuplicate(_primitiveVars, &GradPressure);
		VecZeroEntries(GradPressure);
		int *indices2 = new int[_Nfaces];
		std::iota(indices2, indices2 + _Nfaces, _Nmailles);
		VecSetValues(GradPressure, _Nfaces, indices2, Product, INSERT_VALUES);	
		MatMult(_InvVol, GradPressure, Temporary1);
		VecAXPY(_b,     1, Temporary1); 									
																			
		VecDestroy(& Temporary1);
		VecDestroy(& Temporary2);
		MatDestroy(& Prod); 
		MatDestroy(&ZeroNfaces_Ncells);
		VecDestroy(& GradPressure);
		VecDestroy(& Pressure);	
		delete[] indices2;
		

	}
	if (_nbTimeStep == 0)
		ComputeMinCellMaxPerim(); 
	double numax = _Ndim*2;	//TODO triangles ?
	double dt = _cfl * _minCell / (_maxPerim * max(_uMax,_c) ) ;//max(_uMax,_c)* 8 ); //* numax TODO
	double PreviousTime = _Time.back();
	_Time.push_back(PreviousTime+ dt);
	return dt; 
}

bool EulerBarotropicStaggered::iterateTimeStep(bool &converged)
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

	UpdateDualDensity();
	PetscScalar rhou, rho_sigma, u;
	for (int f=0; f< _Nfaces; f++){
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),f) != _InteriorFaceSet.end() ;
		if (IsInterior){
			PetscInt I = _Nmailles +f;
			VecGetValues(_primitiveVars,1,&I,&u);
			VecGetValues(_DualDensity,1,&f, &rho_sigma);
			rhou = rho_sigma * u;
			VecSetValues(_primitiveVars,1,&I, &rhou, INSERT_VALUES);
		}
	}
	VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk+relaxation*deltaV
	UpdateDualDensity();
	for (int f=0; f< _Nfaces; f++){
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),f) != _InteriorFaceSet.end() ;
		if (IsInterior){
			PetscInt I = _Nmailles +f;
			VecGetValues(_primitiveVars, 1,&I, &rhou);
			VecGetValues(_DualDensity, 1,&f, &rho_sigma);
			u = rhou / rho_sigma ;
			VecSetValues(_primitiveVars,1,&I, &u, INSERT_VALUES);
		}
	}
	return converged;
}

void EulerBarotropicStaggered::terminate(){ 
	MatDestroy(& _LaplacianVelocity);
	MatDestroy(& _Conv); 
	MatDestroy(& _DivRhoU); 

	/* // 	PCDestroy(_pc);
	KSPDestroy(&_ksp);

	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq); */

	WaveStaggered::terminate();
	
}


