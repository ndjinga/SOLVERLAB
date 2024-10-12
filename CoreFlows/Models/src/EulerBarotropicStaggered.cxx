/*
 * EulerBarotropicStaggered.cxx
 */

#include "EulerBarotropicStaggered.hxx"
#include "StiffenedGas.hxx"
#include <cassert>
#include <numeric>

using namespace std;

EulerBarotropicStaggered::EulerBarotropicStaggered(phaseType fluid, pressureEstimate pEstimate, int dim):WaveStaggered( dim, 1.0 ,1.0, MPI_COMM_WORLD ){
;	_Ndim=dim;
	_nVar = 2; 
	_fluides.resize(1);
	_saveVelocity=false; 
	_savePressure=false; 
	_facesBoundinit = false;
	_indexFacePeriodicSet = false;
	_vec_normal=NULL;

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		_Tref=300;
		_Pref=1e5;
		if(fluid==Gas){
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
		if(fluid==Gas){
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

	//Save into the fluid list
	_fluides.resize(1,_compressibleFluid);
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
		int *indices1 = new int[_Nmailles]; //TODO : peut-on faire cela en parallèle ?
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
	
	// *************** Time dependent matrices ************ //
	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);

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
	
	// matrix GRADIENT (we will impose to _Be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_DivTranspose, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_DivTranspose);
	MatSetUp(_DivTranspose);

	// matrix CONVECTION 
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_Conv, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces );
	MatSetFromOptions(_Conv);
	MatSetUp(_Conv);

	// matrice Grad Div
	MatCreate(PETSC_COMM_SELF, & _InvSurface); 
	MatSetSizes(_GradDivTilde, PETSC_DECIDE, PETSC_DECIDE, _Nfaces , _Nfaces );
	MatSetFromOptions(_GradDivTilde);
	MatSetUp(_GradDivTilde);
	MatZeroEntries(_GradDivTilde);

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
	MatZeroEntries(_A);
	MatZeroEntries(_Div);
	
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
	save();//save initial data
}

//TODO is it sufficient to put virtual in front of WaveStagg::Computetimestep

double EulerBarotropicStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme

	MatZeroEntries(_DivRhoU);
	MatZeroEntries(_GradDivTilde);
	MatZeroEntries(_LaplacianPressure);
	MatZeroEntries(_LaplacianVelocity);
	MatZeroEntries(_Conv);
	MatZeroEntries(_DivTranspose); //TODO : vérifier que cela remet bien à zéro les matrices

	if (_timeScheme == Explicit ){ 
		cout << "EulerBarotropicStaggered::computeTimeStep : Début calcul matrice implicite et second membre"<<endl;
		
		if (_mpi_rank ==0){
			// Assembly of matrices 
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				std::vector< int > idCells = Fj.getCellsId();
				Cell Ctemp1 = _mesh.getCell(idCells[0]);
				PetscScalar epsilon, ConvectiveFlux, rhoL, rhoR, u, absConvectiveFlux, MinusabsConvectiveFlux	;
				
				// Metrics
				PetscScalar orientedFaceArea = getOrientation(j,Ctemp1) * Fj.getMeasure();
				PetscScalar orientedMinusFaceArea = -orientedFaceArea;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				PetscScalar det, InvD_sigma, InvPerimeter1, InvPerimeter2;
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
						InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
						InvPerimeter2 = 1/(_perimeters(idCells[1])*Ctemp2.getNumberOfFaces()  );
					}
					InvD_sigma = 1.0/PetscAbsReal(det);
					PetscScalar InvVol2 = 1/( Ctemp2.getMeasure()* Ctemp2.getNumberOfFaces());
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
					PetscScalar MinusorientedFaceArea = - orientedFaceArea;
					PetscScalar FaceArea_upwinding = (abs(u) + _c) * FaceArea/2.0;
					PetscScalar MinusFaceArea_upwinding = -FaceArea_upwinding;
					MatSetValues(_DivRhoU, 1, &idCells[0], 1, &j, &orientedFaceArea_densityMean, ADD_VALUES ); 
					MatSetValues(_DivRhoU, 1, &idCells[1], 1, &j, &MinusorientedFaceArea_densityMean, ADD_VALUES );  
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea_upwinding, ADD_VALUES ); 
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[1], &FaceArea_upwinding, ADD_VALUES );  
					MatSetValues(_LaplacianPressure, 1, &idCells[1], 1, &idCells[1], &MinusFaceArea_upwinding, ADD_VALUES ); 
					MatSetValues(_LaplacianPressure, 1, &idCells[1], 1, &idCells[0], &FaceArea_upwinding, ADD_VALUES );  

					/*************** Momentum conservation equation *****************/
					MatSetValues(_DivTranspose, 1, &j, 1, &idCells[0], &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_DivTranspose, 1, &j, 1, &idCells[1], &orientedMinusFaceArea, ADD_VALUES ); 
					MatSetValues(_Div, 1, &idCells[0], 1, &j, &orientedFaceArea, ADD_VALUES ); 
					MatSetValues(_Div, 1, &idCells[1], 1, &j, &MinusorientedFaceArea, ADD_VALUES ); 
					
					// Convective terms //
					PetscInt jepsilon, L, I;
					// Loop on half diamond cells
					for (int nei =0; nei <Fj.getNumberOfCells(); nei ++){
						Cell K = _mesh.getCell(Fj.getCellsId() [nei]);
						Point xb =  K.getBarryCenter();
						std::vector<int> idFaces = K.getFacesId();
						VecGetValues(_primitiveVars,1,&Fj.getCellsId() [nei],&rhoL);
						//Loop on epsilon the boundary of a dual cell in the fixed half diamond cell 
						for (int Nbepsilon = 0; Nbepsilon < Fj.getNumberOfNodes() ;Nbepsilon ){ 
							ConvectiveFlux = 0 ;		
							Node xsigma = _mesh.getNode( Fj.getNodesId()[Nbepsilon] );
							// For fixed epsilon (and thus the node on sigma defining it ) Loop on the faces that are in the boundary of the cell containing the fixed half diamond cell	
							for (int f =0; f <K.getNumberOfFaces(); f ++){
								bool IsfInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),idFaces[f] ) != _InteriorFaceSet.end() ;
								bool IsfWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),idFaces[f] ) != _WallBoundFaceSet.end() ;
								bool IsfSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),idFaces[f] ) != _SteggerBoundFaceSet.end() ;			
			
								// Compute the flux through epsilon
								Face Facef = _mesh.getFace( idFaces[f] );
								std::vector< int> idCellsOfFacef =  Facef.getCellsId();
								I = _Nmailles + idFaces[f];
								VecGetValues(_primitiveVars,1,&I	,&u);
								PetscScalar psif = 1/K.getNumberOfFaces();               		//TODO à calculer
								if (IsfInterior){												
									std::map<int,int>::iterator it = _FacePeriodicMap.find(f);
									if ( it != _FacePeriodicMap.end()  ){ 
										std::vector< int > idCells_other_Fj =  _mesh.getFace(it->second).getCellsId();
										idCellsOfFacef.push_back( idCells_other_Fj[0]  );
									}
									VecGetValues(_primitiveVars,1,&idCellsOfFacef[1],&rhoR);	
								}
								else if (IsfWallBound )					
										rhoR =  rhoL;
								else if (IsfSteggerBound){ 
										std::map<int,double> boundaryPressure = getboundaryPressure(); 
										std::map<int,double>::iterator it = boundaryPressure.find(j);
										rhoR = boundaryPressure[it->first]; 
								} 
								ConvectiveFlux += ( u *(rhoL + rhoR)/2.0 - (abs(u) + _c)* (rhoR - rhoL) )* psif;
								std::vector< int > idNodes = Facef.getNodesId();
								// -> Search for the unique face that is not sigma that has x_sigma has an extremity -> jepsilon will be the index of this cell
								if ( std::find(idNodes.begin(), idNodes.end(), Fj.getNodesId()[Nbepsilon])!=idNodes.end() && K.getFacesId()[f] != j) 	
									jepsilon =  K.getFacesId()[f];
							}
							
							epsilon = (xb.x() - xsigma.x())*(xb.x() - xsigma.x());
							if (_Ndim > 1){
								epsilon += (xb.y() - xsigma.y() )*(xb.y() - xsigma.y());
								if (_Ndim > 2)
									epsilon += (xb.z() - xsigma.z())*(xb.z() - xsigma.z());
							}
							epsilon = sqrt(epsilon);
							
							ConvectiveFlux *= epsilon/2.0  ;
							MatSetValues(_Conv, 1, &j, 1, &j, &ConvectiveFlux, ADD_VALUES );  		// TODO Should we pass twice ? 
							MatSetValues(_Conv, 1, &j, 1, &jepsilon, &ConvectiveFlux, ADD_VALUES ); 
							absConvectiveFlux = abs(ConvectiveFlux);
							MinusabsConvectiveFlux = -abs(ConvectiveFlux);
							MatSetValues(_LaplacianVelocity, 1, &j, 1, &j, &MinusabsConvectiveFlux, ADD_VALUES ); 
							MatSetValues(_LaplacianVelocity, 1, &j, 1, &jepsilon, &absConvectiveFlux, ADD_VALUES ); 
						}
					} 
				}
				else if (IsWallBound || IsSteggerBound ) { 
					/*************** Matrices related to metrics ***************/		
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
						InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces() ); //TODO ?? pourquoi pas pareil que face intérieure ?InvPerimeter1 = 1/( _perimeters(idCells[0])*Ctemp1.getNumberOfFaces()  );
					}
					MatSetValues(_InvSurface,1, &idCells[0],1, &idCells[0], &InvPerimeter1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &idCells[0],1 ,&idCells[0], &InvVol1, ADD_VALUES );
					MatSetValues(_InvVol, 1, &IndexFace, 1, &IndexFace,  &InvD_sigma, ADD_VALUES); 
					
					/****************** Density conservation equation *********************/
					PetscScalar rhoExt, rhoInt;
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
					PetscScalar MinusFaceArea_upwinding = -(abs(u) + _c) * FaceArea/2.0;
					MatSetValues(_LaplacianPressure, 1, &idCells[0], 1, &idCells[0], &MinusFaceArea_upwinding, ADD_VALUES );
					PetscScalar boundterm = -rhoExt*MinusFaceArea_upwinding; // TODO ??
					VecSetValues(_BoundaryTerms, 1,&idCells[0], &boundterm, INSERT_VALUES );
				}	
			}
		}
		MatAssemblyBegin(_InvSurface, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
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
		MatAssemblyBegin(_Conv, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_Conv, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_LaplacianVelocity, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_LaplacianVelocity, MAT_FINAL_ASSEMBLY);


		// Minimum size of mesh volumes
		int *indices1 = new int[_Nmailles]; 
		std::iota(indices1, indices1+_Nmailles, 0);
		VecMax(_primitiveVars, indices1, &_rhoMax);
		int *indices2 = new int[_Nfaces]; 
		std::iota(indices2, indices2+_Nfaces, _Nmailles );
		VecMax(_primitiveVars, indices2, &_uMax);
		delete[] indices1;
		delete[] indices2;

		MatScale(_DivTranspose, -1.0);
		MatMatMatMult(_DivTranspose,_InvSurface, _Div , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &_GradDivTilde); 
		MatScale(_DivRhoU, -1.0);
		MatScale(_DivTranspose, -1.0);
		MatScale(_GradDivTilde, _c*_rhoMax) ; 
		
		// _A = ( (u+c) _Laplacian   ;  -Div(\rho .)         )
		//      (-Conv  ;  -c (\tilde rho) MinusGrad (1/|dK| Div ) 
		Mat G[4];
		G[0] = _LaplacianPressure;
		G[1] = _DivRhoU;
		G[2] = _Conv ;
		G[3] = _GradDivTilde;
		// TODO ajouter minus grad _DivTranspose _compressibleFluid->getPressure() ???
		// TODO ajouter Laplacian velocity

		MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
		Mat Prod;
		MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
		MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, SAME_NONZERO_PATTERN); //TODO Coment sortir les redondances avec Prod

		Vec Prod2;
		VecDuplicate(_BoundaryTerms, &Prod2);
		MatMult(_InvVol, _BoundaryTerms, Prod2);  
		MatMult(_A,_primitiveVars, _b); 
		VecAXPY(_b,     1, Prod2);
		VecDestroy(& Prod2);	
	}
	if (_nbTimeStep == 0)
		ComputeMinCellMaxPerim();

	return _cfl * _minCell / (_maxPerim * max(_uMax,_c) );
}


void EulerBarotropicStaggered::terminate(){ 
	delete[]_vec_normal;
	WaveStaggered::terminate();
	MatDestroy(& _LaplacianVelocity);
	MatDestroy(& _Conv); 
	MatDestroy(& _DivRhoU); 
	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);

	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq);
}


