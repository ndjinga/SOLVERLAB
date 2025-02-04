/*
 * EulerBarotropicStaggered.cxx
 */

#include "EulerBarotropicStaggered.hxx"
#include "StiffenedGas.hxx"
#include "BarotropicLaw.hxx"
#include <cassert>
#include <numeric>

using namespace std;

EulerBarotropicStaggered::EulerBarotropicStaggered(phaseTypeStaggered fluid, pressureEstimate pEstimate, double a, double gamma, int dim):WaveStaggered( dim, 1.0 ,1.0, MPI_COMM_WORLD ){
;	_Ndim=dim;
	_nVar = 2; 
	_saveVelocity=false; 
	_savePressure=false; 
	_facesBoundinit = false;
	_indexFacePeriodicSet = false;
	_vec_normal=NULL;
	_Time.push_back(0);
	_compressibleFluid = new BarotropicLaw(a, gamma);
	_fluides.resize(0);
	_fluides.push_back(_compressibleFluid);
	if (_Ndim == 3){	
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!EulerBarotropicStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!EulerBarotropicStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("EulerBarotropicStaggered pas dispo en 3D, arret de calcul");				
	}		
}

std::vector<double>  EulerBarotropicStaggered::getTimeEvol(){
	return _Time;
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

	// matrice des Inverses de Surfaces
	MatCreate(PETSC_COMM_SELF, & _InvSurface); 
	MatSetSizes(_InvSurface, PETSC_DECIDE, PETSC_DECIDE, _Nmailles , _Nmailles );
	MatSetFromOptions(_InvSurface);
	MatSetUp(_InvSurface);
	MatZeroEntries(_InvSurface);

	AssembleMetricsMatrices();
	
	// matrice DIVERGENCE (|K|div(u))
	MatCreate(PETSC_COMM_SELF, & _Div); 
	MatSetSizes(_Div, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_Div);
	MatSetUp(_Div);

	// matrix minus GRADIENT (we will impose to be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_DivTranspose, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_DivTranspose);
	MatSetUp(_DivTranspose);

	// *************** Time dependent matrices ************ //
	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);
	MatZeroEntries(_A);

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

	// Vector BoundaryTerms for Pressure
	VecCreate(PETSC_COMM_SELF, & _BoundaryTerms); 
	VecSetSizes(_BoundaryTerms, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_BoundaryTerms);
	VecSetUp(_BoundaryTerms);
	VecZeroEntries(_BoundaryTerms);
	
	// Vector CONVECTION 
	VecCreate(PETSC_COMM_SELF, & _Conv); 
	VecSetSizes(_Conv, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_Conv);
	VecSetUp(_Conv);

	// matrix LAPLACIAN VeLOCITY(without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianVelocity); 
	MatSetSizes(_LaplacianVelocity, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces ); 
	MatSetFromOptions(_LaplacianVelocity);
	MatSetUp(_LaplacianVelocity);


	if(_system){
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


void EulerBarotropicStaggered::AssembleMetricsMatrices(){
	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		std::vector< int > NodesFj =  Fj.getNodesId();
 		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		std::vector< int > idFaces = Ctemp1.getFacesId();
		PetscScalar det, InvD_sigma, InvPerimeter1, InvPerimeter2, deltaY	;
		PetscInt IndexFace = _Nmailles + j;
		PetscScalar InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;	

		for (int f =0; f <Ctemp1.getNumberOfFaces(); f ++){			
			Face Facef = _mesh.getFace( idFaces[f] );
			std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
			if ( ( (std::find(idNodesOfFacef.begin(), idNodesOfFacef.end(), NodesFj[0])!=idNodesOfFacef.end()) || (std::find(idNodesOfFacef.begin(), idNodesOfFacef.end(), NodesFj[1])!=idNodesOfFacef.end()) ) && (idFaces[f] != j) ){
				deltaY = Facef.getMeasure();
			}	
		}

		if (IsInterior){
			std::map<int,int>::iterator it = _FacePeriodicMap.find(j);
			if ( it != _FacePeriodicMap.end()  ){ 
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
				det = (Ctemp1.x() - vertex.x() )* (Ctemp2.y() - vertex.y() ) - (Ctemp1.y() - vertex.y() )* (Ctemp2.x() - vertex.x() );
				det = Fj.getMeasure() *deltaY; //TODO only cartesian 
				InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
				InvPerimeter2 = 1.0/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
			}
			InvD_sigma =1.0/(PetscAbsReal(det));
			PetscScalar InvVol2 = 1.0/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());
			MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1 , ADD_VALUES );
			MatSetValues(_InvVol, 1, &idCells[1],1 ,&idCells[1], &InvVol2, ADD_VALUES );
			MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 	
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
				det = Fj.getMeasure() *deltaY /2.0 ;// TODO only cartesian 
				InvD_sigma = 1.0/(PetscAbsReal(det)); 	
				InvPerimeter1 = 1.0/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces() ); 
			}
			MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
			MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
			MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES);  
		}
	}
	MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(_InvSurface,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
}




void EulerBarotropicStaggered::UpdateDualDensity(){
	VecZeroEntries(_DualDensity);

	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		PetscScalar detL, detR, D_sigmaL, D_sigmaR, rho_sigma,inv_rho_sigma, rhoL, rhoR, usigma;

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
				detL = (Ctemp2.x() - Ctemp1.x())/2.0;
				detR = (Ctemp2.x() - Ctemp1.x())/2.0;
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
			VecGetValues(_primitiveVars,1,&j,&usigma);
			rho_sigma = (rhoL * D_sigmaL + rhoR*D_sigmaR)/(D_sigmaL +  D_sigmaR) ;
			
			//rho_sigma = 1.0; //Découplage advectionrho_sigma
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
	_c =  _compressibleFluid->vitesseSon(_rhoMax); 	 // Découplage Burgers : c =0
	
	MatZeroEntries(_DivRhoU);
	MatZeroEntries(_LaplacianPressure);
	VecZeroEntries(_Conv);
	MatZeroEntries(_Div); 
	MatZeroEntries(_DivTranspose); 
	MatZeroEntries(_LaplacianVelocity);
	
	if (_timeScheme == Explicit ){ 
		if (_mpi_rank ==0){
			// Assembly of matrices 
			for (int j=0; j<_Nfaces; j++){	 
				Face Fj = _mesh.getFace(j);
				std::vector< int > idCells = Fj.getCellsId();
				Cell Ctemp1 = _mesh.getCell(idCells[0]);
				PetscScalar Convection_integral, rhoL, rhoR, u, u_sigma, det, InvD_sigma, InvPerimeter1, InvPerimeter2;	
				
				PetscScalar orientedFaceArea = getOrientation(j,Ctemp1) * Fj.getMeasure();
				PetscScalar orientedMinusFaceArea = -orientedFaceArea;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				PetscInt IndexFace = _Nmailles + j;
				//PetscScalar InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());

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

					/****************** Density conservation equation *********************/
					VecGetValues(_primitiveVars,1,&idCells[0],&rhoL);
					VecGetValues(_primitiveVars,1,&idCells[1],&rhoR);
					VecGetValues(_primitiveVars,1,&IndexFace,&u_sigma);

					PetscScalar orientedFaceArea_densityMean = orientedFaceArea * (rhoL + rhoR)/2.0 ;
					PetscScalar MinusorientedFaceArea_densityMean = -orientedFaceArea_densityMean;
					PetscScalar FaceArea_upwinding = (abs(u)+ _c ) * FaceArea/2.0; 
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
					
					// Convective terms (WARNING !!!!! is not computed in 3 space dimensions)
					PetscInt I;
					Convection_integral= 0;
					for (int nei =0; nei <2; nei ++){ // we are in the case where an interior face has two neighbours
						Cell K = _mesh.getCell(idCells[nei]); 
						Point BarycenterK =  K.getBarryCenter();
						std::vector<int> idFaces = K.getFacesId();
						VecGetValues(_primitiveVars,1,&idCells[nei],&rho);

						//If  the face is periodic, recover the geometric informations of the associated periodic cell 
						std::vector< int > NodesFj = ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[nei]    ? _mesh.getFace(_FacePeriodicMap.find(j)->second ).getNodesId() :  Fj.getNodesId(); 
						Point BarycenterFj =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[nei] 	 ? _mesh.getFace(_FacePeriodicMap.find(j)->second ).getBarryCenter() :  Fj.getBarryCenter();
						
						for (int f =0; f <K.getNumberOfFaces(); f ++){
							bool IsfInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),idFaces[f] ) != _InteriorFaceSet.end() ;
							bool IsfWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),idFaces[f] ) != _WallBoundFaceSet.end() ;
							bool IsfSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),idFaces[f] ) != _SteggerBoundFaceSet.end();

							Face Facef = _mesh.getFace( idFaces[f] );
							std::vector< int> idCellsOfFacef =  Facef.getCellsId();
							std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
							I = _Nmailles + idFaces[f];
							VecGetValues(_primitiveVars,1,&I	,&u);

							if (_Ndim==1)
								Convection_integral -= getOrientation(j, K) * u*u * rho/2.0; 
							if (_Ndim == 2){
								if ( ( (std::find(idNodesOfFacef.begin(), idNodesOfFacef.end(), NodesFj[0])==idNodesOfFacef.end()) && (std::find(idNodesOfFacef.begin(), idNodesOfFacef.end(), NodesFj[1])==idNodesOfFacef.end()) ) || (idFaces[f] ==j)  ){
									Convection_integral -=  getOrientation(j, K) * Fj.getMeasure()* u*u * rho/2.0 ;
								}
							}

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
						}
					} 
					VecSetValues(_Conv, 1, &IndexFace, &Convection_integral, ADD_VALUES ); 
				}
				else if (IsWallBound || IsSteggerBound ) { 
					/****************** Density conservation equation *********************/
					PetscScalar rhoExt, rhoInt;
					VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);
					//VecGetValues(_primitiveVars,1,&IndexFace,&u);
					
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
					
					PetscScalar MinusFaceArea_upwinding = -( abs(u)+_c ) * FaceArea/2.0; 
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea_upwinding, ADD_VALUES );
					PetscScalar boundterm = -rhoExt*MinusFaceArea_upwinding;
					VecSetValues(_BoundaryTerms, 1,&idCells[0], &boundterm, INSERT_VALUES );
				}	
			}
		}

		MatAssemblyBegin(_DivRhoU,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_DivRhoU, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_LaplacianPressure,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_LaplacianPressure, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_BoundaryTerms);
		VecAssemblyEnd(_BoundaryTerms);

		MatAssemblyBegin(_DivTranspose, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_DivTranspose, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_Div, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_Div, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_Conv);
		VecAssemblyEnd(_Conv);	
		
		/***********Assembling the matrix _A such that : *************/
		// _A = ( (u+c) _Laplacian   ;                 -Div(\rhomean .)                  )
		//      (0_FacesxCells   ; -Conv+ c (\tilde rho) MinusGrad 1/|dK| Div + LaplacianVelocity	)

		MatMatMatMult(_DivTranspose,_InvSurface, _Div , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &_GradDivTilde); // -grad (inv_Surf) Div// -(-grad (inv_Surf) Div) = grad (inv_Surf) Div
		MatScale(_DivRhoU, -1.0);	
		VecScale(_Conv, -1.0);
		MatScale(_GradDivTilde, -1. *( _c*_rhoMax + _uMax)/2.0 );
		//MatScale(_DivRhoU, 0);	// Découplage burgers

		Mat G[4], ZeroNfaces_Ncells, Prod;
		MatDuplicate(_DivTranspose, MAT_DO_NOT_COPY_VALUES, &ZeroNfaces_Ncells); 	
		MatZeroEntries(ZeroNfaces_Ncells);
		G[0] = _LaplacianPressure;
		G[1] = _DivRhoU;
		G[2] = ZeroNfaces_Ncells ;
		G[3] = _GradDivTilde; 
		MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
		MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
		MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, UNKNOWN_NONZERO_PATTERN); 						
		MatMult(_A,_primitiveVars, _b); 
		
		/***********Adding boundary terms to density equaiton, convection term and pressure gradient to momentum equation **************/
		Vec Temporary1, Temporary2, Pressure, GradPressure;
		VecCreate(PETSC_COMM_SELF, & Temporary1);
		VecCreate(PETSC_COMM_SELF, & Temporary2);
		VecSetSizes(Temporary1, PETSC_DECIDE, _globalNbUnknowns);
		VecSetSizes(Temporary2, PETSC_DECIDE, _Nfaces);
		VecSetFromOptions(Temporary1);
		VecSetFromOptions(Temporary2);
		//add boundary terms : AU^n + V^{-1}Boundterms //
		MatMult(_InvVol, _BoundaryTerms, Temporary1); 
		VecAXPY(_b,     1, Temporary1); 
		//add convection terms  : AU^n + V^{-1}Boundterms + V^{-1}_Conv//
		MatMult(_InvVol, _Conv, Temporary1); 
		VecAXPY(_b,     1, Temporary1); 
		//Extract pressure to multiply by _DivTranspose = -grad //
		VecCreate(PETSC_COMM_SELF, & Pressure); 
		VecSetSizes(Pressure, PETSC_DECIDE, _Nmailles); 
		VecSetFromOptions(Pressure);
		VecSetUp(Pressure);
		for (int i=0; i <_Nmailles; i++){
			VecGetValues(_primitiveVars,1,&i,&rho);
			p =  _compressibleFluid->getPressure(rho); 				
			VecSetValues(Pressure, 1,&i, &p, INSERT_VALUES );
		}														
		//add pressure gradient : AU^n+ V^{-1}Boundterms - V^{-1}_Conv - V^{-1}GradPressure//
		MatMult(_DivTranspose, Pressure, Temporary2); 						
		double *Product = new double[_Nfaces]; //Can we avoid creating a pointer only to insert in the bigger Vector ?
		PetscScalar gradp;
		for (int i=0; i <_Nfaces; i++){
			VecGetValues(Temporary2,1,&i,&gradp);
			Product[i] = gradp; // Découplage gradp=0;
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
	double numax = _Ndim*2;	
	double dt = _cfl * _minCell / (_maxPerim * (_uMax + _c) * 2 ) ;
	double PreviousTime = _Time.back();
	_Time.push_back(PreviousTime+ dt);
	
	return dt; 
}


//################ Raviar-Thomas ##############//
std::vector<double> EulerBarotropicStaggered::ReferenceBasisFunctionRaviartThomas(int i, Point Xhat){
	std::vector<double> Psihat(_Ndim);
	if (_Ndim ==2){
		if (i ==0){
			Psihat[0] = (Xhat.x() + 1)/2.0;;
			Psihat[1] = 0;
		}
		else if (i ==1){
			Psihat[0] = 0;
			Psihat[1] = (Xhat.y() + 1)/2.0;
		}
		else if (i ==2){	
			Psihat[0]=  -(Xhat.x() - 1)/2.0;
			Psihat[1] = 0;
		}
		else if (i ==3){
			Psihat[0] = 0;
			Psihat[1] = -(Xhat.y() - 1)/2.0;	
		}
	}
	return Psihat;
}

std::vector<double>  EulerBarotropicStaggered::Gradient_ReferenceBasisFunctionRaviartThomas(int i){
	std::vector<double>  Gradient_Psihat(_Ndim*_Ndim); 
	if (_Ndim ==2){
		if (i ==0){
			Gradient_Psihat[0] = 1/2.0;
			Gradient_Psihat[1] = 0;
			Gradient_Psihat[2] = 0;
			Gradient_Psihat[3] = 0;
		}
		else if (i ==1){
			Gradient_Psihat[0] = 0;
			Gradient_Psihat[1] = 0;
			Gradient_Psihat[2] = 0;
			Gradient_Psihat[3] = 1/2.0;

		}
		else if (i ==2){	
			Gradient_Psihat[0] = -1/2.0;
			Gradient_Psihat[1] = 0;
			Gradient_Psihat[2] = 0;
			Gradient_Psihat[3] = 0;
		}
		else if (i ==3){
			Gradient_Psihat[0] = 0;
			Gradient_Psihat[1] = 0;
			Gradient_Psihat[2] = 0;
			Gradient_Psihat[3] = -1/2.0;
		}
	}
	return Gradient_Psihat;	
}

Point EulerBarotropicStaggered::xToxhat(Cell K, Point X, std::vector<Node> K_Nodes){
	Point Xhat;
	Point Xb = K.getBarryCenter();
	if (_Ndim ==2){
		if (X .x() == Xb.x() && X.y() == Xb.y() ) //if X is the center of gravity
			Xhat = Point(0,0,0);
		
		//Nodes related
		else if (X.x() == K_Nodes[0].x() && X.y() == K_Nodes[0].y() ) // if X is the Node X_1
			Xhat = Point(-1,-1,0);
		else if (X.x() == K_Nodes[1].x() && X.y() == K_Nodes[1].y() ) // if X is the Node X_2
			Xhat = Point(1,-1,0);
		else if (X.x() == K_Nodes[2].x() && X.y() == K_Nodes[2].y() ) // if X is the Node X_3
			Xhat = Point(1,1,0);
		else if (X.x() == K_Nodes[3].x() && X.y() == K_Nodes[3].y() ) // if X is the Node X_4
			Xhat = Point(-1,1,0);

		//Faces related
		else if (X.x() == ( K_Nodes[0].x() +  K_Nodes[1].x())/2.0 && X.y() == ( K_Nodes[0].y() +  K_Nodes[1].y())/2.0  ) // if X is the middle point of the first face
			Xhat = Point(0,-1,0);
		else if (X.x() == ( K_Nodes[1].x() +  K_Nodes[2].x())/2.0 && X.y() == ( K_Nodes[1].y() +  K_Nodes[2].y())/2.0 )  // if X is the middle point of the second face
			Xhat = Point(1,0,0);
		else if (X.x() == ( K_Nodes[2].x() +  K_Nodes[3].x())/2.0 && X.y() == ( K_Nodes[2].y() +  K_Nodes[3].y())/2.0 )  // if X is the middle point of the third face
			Xhat = Point(0,1,0);
		else if (X.x() == ( K_Nodes[3].x() +  K_Nodes[0].x())/2.0 && X.y() == ( K_Nodes[3].y() +  K_Nodes[0].y())/2.0 )  // if X is the middle point of the fourth face
			Xhat = Point(-1,0,0);
	}	
	return Xhat;
}

std::vector<double> EulerBarotropicStaggered::JacobianTransfor_K_X(Point X, std::vector<Node> K_Nodes){
	std::vector<double> JacobianTransfor_K(_Ndim * _Ndim);
	if (_Ndim ==2){
	//first column : should only depend on Xhat_x
	JacobianTransfor_K[0] = 1/4.0 *( (K_Nodes[0].x()  - K_Nodes[1].x())*(1 - X.y()) + (K_Nodes[3].x()  - K_Nodes[2].x())*(1 + X.y())  );
	JacobianTransfor_K[2] = 1/4.0 *( (K_Nodes[0].y()  - K_Nodes[1].y())*(1 - X.y()) + (K_Nodes[3].y()  - K_Nodes[2].y())*(1 + X.y())  );
	//second column: should only depend on Xhat_y
	JacobianTransfor_K[1] = 1/4.0 *( (K_Nodes[3].x()  - K_Nodes[0].x())*(1 - X.x()) + (K_Nodes[2].x()  - K_Nodes[1].x())*(1 + X.x())  );
	JacobianTransfor_K[3] = 1/4.0 *( (K_Nodes[3].y()  - K_Nodes[0].y())*(1 - X.x()) + (K_Nodes[2].y()  - K_Nodes[1].y())*(1 + X.x())  );
	double absJacobian =  abs(JacobianTransfor_K[0] * JacobianTransfor_K[3] - JacobianTransfor_K[2]* JacobianTransfor_K[1] ); 
	for (int k =0; k< _Ndim * _Ndim; k++)
		JacobianTransfor_K[k] /= absJacobian;
	}
	return JacobianTransfor_K;
}

std::vector<double> EulerBarotropicStaggered::PhysicalBasisFunctionRaviartThomas(Cell K, Face Facef,int f, Point X){
	std::vector<double> PhysicalPsif_in_X(_Ndim, 0.0);
	std::vector<int> idNodesofCellK = K.getNodesId();
	std::vector<Node> K_Nodes;
	for (int i=0; i < idNodesofCellK.size(); i++)
		K_Nodes.push_back(_mesh.getNode(idNodesofCellK[i]) );

	Point Xhat_f =  xToxhat(K,  Facef.getBarryCenter(), K_Nodes) ; 
	Point Xhat = xToxhat(K, X, K_Nodes); 
	std::vector<double>  JacobianTransfor_K_Xhatf = JacobianTransfor_K_X(Xhat_f, K_Nodes );
	std::vector<double>  JacobianTransfor_K_Xhat = JacobianTransfor_K_X(Xhat, K_Nodes );
	
	for (int  j = 0; j < K.getNumberOfFaces(); j ++){
		std::vector<double> ReferencePsij = ReferenceBasisFunctionRaviartThomas(j, Xhat_f);
		std::vector<double> PhysicalPsij(_Ndim);
		for (int k =0; k < _Ndim ; k++){
			PhysicalPsij[k] = 0;
			for (int l =0; l < _Ndim ; l++)
				PhysicalPsij[k] += JacobianTransfor_K_Xhatf[k* _Ndim +l] * ReferencePsij[l]* Facef.getMeasure() /2.0 ;
		}
		double cdot =0;
		std::map<int, std::vector<double>  >::iterator it = _vec_sigma.find(f);
		for (int e=0; e <_Ndim; e++)
			cdot += PhysicalPsij[e] * it->second[e] 	; 
		if (abs(cdot) !=0){
			std::vector<double> ReferencePsif_in_X = ReferenceBasisFunctionRaviartThomas(j, Xhat );
			for (int k =0; k < _Ndim ; k++){	
				for (int l =0; l < _Ndim ; l++)
					PhysicalPsif_in_X[k] += JacobianTransfor_K_Xhat[k*_Ndim + l] * ReferencePsif_in_X[l] * Facef.getMeasure()/2.0 ;//TODO why 1/2.0 ?
			}
		}
	}

	return PhysicalPsif_in_X;
}


//TODO ; for now only on affine meshes and need to add B^{-t}
std::vector<double> EulerBarotropicStaggered::Gradient_PhysicalBasisFunctionRaviartThomas(Cell K, Face Facef, int f, Point X){
	std::vector<double> Gradient_PhysicalPsif_in_X(_Ndim* _Ndim);
	for (int k =0; k < _Ndim* _Ndim ; k++)
		Gradient_PhysicalPsif_in_X[k+ k*_Ndim ] = 0;

	std::vector<int> idNodesofCellK = K.getNodesId();
	std::vector<Node> K_Nodes;
	for (int i=0; i < idNodesofCellK.size(); i++)
		K_Nodes.push_back(_mesh.getNode(idNodesofCellK[i]) );
	
	Point Xhat_f =  xToxhat(K,  Facef.getBarryCenter(), K_Nodes) ; 
	Point Xhat = xToxhat(K, X, K_Nodes); 
	std::vector<double>  JacobianTransfor_K_Xhatf = JacobianTransfor_K_X(Xhat_f, K_Nodes );
	std::vector<double>  JacobianTransfor_K_Xhat = JacobianTransfor_K_X(Xhat, K_Nodes ); 

	for (int  j = 0; j < K.getNumberOfFaces(); j ++){
		std::vector<double> ReferencePsij = ReferenceBasisFunctionRaviartThomas(j, Xhat_f);
		std::vector<double> PhysicalPsij(_Ndim);
		for (int k =0; k < _Ndim ; k++){
			PhysicalPsij[k] = 0;
			for (int l =0; l < _Ndim ; l++)
				PhysicalPsij[k] += JacobianTransfor_K_Xhatf[k*_Ndim   + l] * ReferencePsij[l]* Facef.getMeasure() ;
		}
		double cdot =0;
		std::map<int, std::vector<double>  >::iterator it = _vec_sigma.find(f);
		for (int e=0; e <_Ndim; e++)
			cdot += PhysicalPsij[e] * it->second[e] 	; 
		if (abs(cdot) !=0){
			std::vector<double>  Gradient_ReferencePsif_in_X = Gradient_ReferenceBasisFunctionRaviartThomas(j);
			for (int i =0; i < _Ndim ; i++){
				for (int j =0; j < _Ndim ; j++){
					for (int l =0; l < _Ndim ; l++)
						Gradient_PhysicalPsif_in_X[i*_Ndim + j] += JacobianTransfor_K_Xhat[i*_Ndim  +  l] * Gradient_ReferencePsif_in_X[l*_Ndim + j] * Facef.getMeasure()/2.0 ;//TODO why 1/2.0 ?
				}
			}		
		}
	}
	return Gradient_PhysicalPsif_in_X;
}


std::vector<double> EulerBarotropicStaggered::VelocityRaviartThomas_at_point_X(Cell K, Point X){
	std::vector<double> VelocityRT(_Ndim,0.0);
	std::vector< int > idFaces = K.getFacesId();
	PetscInt IndexFace;
	PetscScalar u;
	for (int f=0; f< K.getNumberOfFaces(); f++){
		PetscInt IndexFace = _Nmailles + idFaces[f];
		std::vector<double> Psif = PhysicalBasisFunctionRaviartThomas(K,_mesh.getFace(idFaces[f]), idFaces[f], X);
		VecGetValues(_primitiveVars,1,&IndexFace,&u);
		for (int k=0; k< _Ndim; k++)
			VelocityRT[k] += u * Psif[k];
	}
	return VelocityRT;
 }


 std::vector<double> EulerBarotropicStaggered::TensorProduct(std::vector<double> &u, std::vector<double> &v){
	std::vector<double> tensorproduct(u.size() * v.size());
	for (int i=0; i< _Ndim; i++){
		for (int j=0; j< _Ndim; j++){
			tensorproduct[i*_Ndim + j] = u[i] * v[j];
		}
	}
	return tensorproduct;
 }

double EulerBarotropicStaggered::Contraction(std::vector<double> &u, std::vector<double> &v){
	double contraction;
	for (int i=0; i< _Ndim; i++){
		for (int j=0; j< _Ndim; j++){
			contraction += u[i*_Ndim + j]* v[i*_Ndim + j];
		}
	}
	return contraction;
 }


bool EulerBarotropicStaggered::iterateTimeStep(bool &converged){
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

	UpdateDualDensity(); // \rho^n_K -> \rho^n_\sigma
	PetscScalar rhou, rho_sigma, u;
	for (int f=0; f< _Nfaces; f++){
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),f) != _InteriorFaceSet.end() ;
		if (IsInterior){
			PetscInt I = _Nmailles +f;
			VecGetValues(_primitiveVars,1,&I,&u);
			VecGetValues(_DualDensity,1,&f, &rho_sigma);
			rhou = rho_sigma * u;
			VecSetValues(_primitiveVars,1,&I, &rhou, INSERT_VALUES); // (\rho^n, u^n) -> (\rho^n, \rho u^n)
		}
	}
	VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk+relaxation*deltaV : (\rho^{n+1}, \rho u^{n+1}) =  Delta(q) + \rho^n u^n
	UpdateDualDensity(); // \rho^{n+1}_K -> \rho^{n+1}_\sigma
	for (int f=0; f< _Nfaces; f++){
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),f) != _InteriorFaceSet.end() ;
		if (IsInterior){
			PetscInt I = _Nmailles +f;
			VecGetValues(_primitiveVars, 1,&I, &rhou);
			VecGetValues(_DualDensity, 1,&f, &rho_sigma);
			u = rhou / rho_sigma ;
			VecSetValues(_primitiveVars,1,&I, &u, INSERT_VALUES); // (\rho^{n+1}, \rho^{n+1}u^{n+1}) ->(\rho^{n+1}, u^{n+1})
		}
	}
	return converged;
}

void EulerBarotropicStaggered::terminate(){ 
	MatDestroy(& _LaplacianVelocity);
	VecDestroy(& _Conv); 
	MatDestroy(& _DivRhoU); 

	/* // 	PCDestroy(_pc);
	KSPDestroy(&_ksp);

	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq); */

	WaveStaggered::terminate();
	
}


