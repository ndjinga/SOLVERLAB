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
	delete[] initialFieldVelocity;
	delete[] initialFieldPressure;
	delete[] indices1;
	delete[] indices2;
	
	VecAssemblyBegin(_primitiveVars);
	VecAssemblyEnd(_primitiveVars);

	VecCreate(PETSC_COMM_SELF, & _DualDensity);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_DualDensity, PETSC_DECIDE, _Nfaces);
	VecSetFromOptions(_DualDensity);
	VecZeroEntries(_DualDensity);

	VecCreate(PETSC_COMM_SELF, & _velocityVec);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_velocityVec, PETSC_DECIDE, _Nfaces);
	VecSetFromOptions(_velocityVec);
	VecZeroEntries(_velocityVec);

	//************ Time independent matrices ********** //
	// matrice des Inverses Volumes V^{-1}
	MatCreate(PETSC_COMM_SELF, &_InvVol); 
	MatSetSizes(_InvVol, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_InvVol);
	MatSetUp(_InvVol);
	MatZeroEntries(_InvVol);

	// matrix minus GRADIENT (we will impose to be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_DivTranspose, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_DivTranspose);
	MatSetUp(_DivTranspose);

	// matrice des Inverses de Surfaces
	MatCreate(PETSC_COMM_SELF, & _InvSurface); 
	MatSetSizes(_InvSurface, PETSC_DECIDE, PETSC_DECIDE, _Nmailles , _Nmailles );
	MatSetFromOptions(_InvSurface);
	MatSetUp(_InvSurface);
	MatZeroEntries(_InvSurface);

	// matrice DIVERGENCE (|K|div(u))
	MatCreate(PETSC_COMM_SELF, & _Div); 
	MatSetSizes(_Div, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_Div);
	MatSetUp(_Div);

	// *************** Time dependent matrices ************ //
	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);

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

	// matrix LAPLACIAN Pressure (without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianPressure); 
	MatSetSizes(_LaplacianPressure, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles ); 
	MatSetFromOptions(_LaplacianPressure);
	MatSetUp(_LaplacianPressure);

	AssembleMetricsMatrices();

	if(_system){
		cout << "Variables primitives initiales : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}

	/* if(_mpi_size>1 && _mpi_rank == 0)
    	VecCreateSeq(PETSC_COMM_SELF, _globalNbUnknowns, &_primitiveVars_seq);//For saving results on proc 0
    VecScatterCreateToZero(_primitiveVars,&_scat,&_primitiveVars_seq); */

	createKSP();
	PetscPrintf(PETSC_COMM_WORLD,"SOLVERLAB Newton solver ");
	*_runLogFile << "SOLVERLAB Newton solver" << endl;
	_runLogFile->close();

	_initializedMemory=true;
	UpdateDualDensity();
	save();//save initial data
}



void EulerBarotropicStaggered::UpdateDualDensity(){
	VecZeroEntries(_DualDensity);

	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		double D_sigmaL, D_sigmaR, rho_sigma, rhoL, rhoR ;

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;		
		PetscInt I = _Nmailles + j;
		VecGetValues(_primitiveVars,1,&idCells[0],&rhoL);

		std::vector< int > nodes =  Fj.getNodesId();
		Node vertex1,vertex2;
		vertex1  = _mesh.getNode( nodes[0] );
		if (_Ndim == 2) vertex2 = _mesh.getNode( nodes[1] );
		D_sigmaL = (_Ndim == 2) ? fabs( (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() ) ) /2.0 : fabs(vertex1.x() - Ctemp1.x())/2.0;
				
		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  )
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			D_sigmaR = (_Ndim == 2) ? fabs( (Ctemp2.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp2.y() - vertex1.y() )* (vertex2.x() - vertex1.x() ) ) /2.0 :  fabs(Ctemp2.x() - vertex1.x())/2.0;
		
			VecGetValues(_primitiveVars,1,&idCells[1],&rhoR);
			VecSetValue(_DualDensity, j, (rhoL * D_sigmaL + rhoR*D_sigmaR)/(D_sigmaL +  D_sigmaR), INSERT_VALUES );
		}
		else if (IsSteggerBound)  VecSetValue(_DualDensity, j, (rhoL + getboundaryPressure().find(j)->second )/2.0, INSERT_VALUES );
		else if (IsWallBound)      VecSetValue(_DualDensity, j, rhoL , INSERT_VALUES );

	}
	VecAssemblyBegin(_DualDensity);
	VecAssemblyEnd(_DualDensity);
}

double EulerBarotropicStaggered::getOrientationNode(int n, int j) {
	assert(_Ndim == 2); 
	std:vector<double> n_sigma_perp(_Ndim, 0.0), XnMinusXsigma(_Ndim, 0.0);
	n_sigma_perp[0] = -_vec_sigma.find(j)->second[1];
	n_sigma_perp[1] =  _vec_sigma.find(j)->second[0];
	XnMinusXsigma[0] = _mesh.getNode(n).x() - _mesh.getFace(j).getBarryCenter().x() ;
	XnMinusXsigma[1] = _mesh.getNode(n).y() - _mesh.getFace(j).getBarryCenter().y() ;
	
	double dotprod = 0;
	if ( _mesh.getFace(j).getNodesId()[0] == n  ||  _mesh.getFace(j).getNodesId()[1] == n ){
		for (int idim = 0; idim < _Ndim; ++idim)
			dotprod += n_sigma_perp[idim] * XnMinusXsigma[idim] * 2.0/_mesh.getFace(j).getMeasure(); 
	}
	return dotprod ; 
}

double EulerBarotropicStaggered::computeTimeStep(bool & stop){

	Rhomax_Umax_Cmax();
	VecZeroEntries(_Conv);
	VecZeroEntries(_BoundaryTerms);
	MatZeroEntries(_A);
	PetscScalar rho,p, u, Convection, rhoInt, rhoExt, q;
	

	if (_timeScheme == Explicit ){ 
		for (int j=0; j<_Nfaces ; j++){	 
			Face Fj = _mesh.getFace(j);
			std::vector< int > idCells = Fj.getCellsId();
			Cell Ctemp1 = _mesh.getCell(idCells[0]);
		
			PetscInt IndexFace = _Nmailles + j;
			PetscInt I;
			VecGetValues(_primitiveVars,1,&IndexFace,&q);

			bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
			bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
			bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;			

			if (IsInterior){
				if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  ) 
					idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j) ->second).getCellsId()[0]  );
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
		
				// -DivRhoU_{idcell[1], j}, -DivRhoU_{idcell[0], j}, (IndexFace = ncells + j)
				MatSetValue(_A, idCells[0], IndexFace,  -getOrientation(j,Ctemp1) * Fj.getMeasure() , ADD_VALUES ); 
				MatSetValue(_A, idCells[1], IndexFace,  -getOrientation(j,Ctemp2) * Fj.getMeasure() , ADD_VALUES );  
				// LaplacianPressure
				MatSetValue(_A, idCells[0], idCells[0], -(abs( _Velocity(j) )+ _c ) * Fj.getMeasure()/2.0, ADD_VALUES ); 
				MatSetValue(_A, idCells[0], idCells[1],  (abs( _Velocity(j) )+ _c ) * Fj.getMeasure()/2.0, ADD_VALUES );  
				MatSetValue(_A, idCells[1], idCells[1], -(abs( _Velocity(j) )+ _c ) * Fj.getMeasure()/2.0, ADD_VALUES ); 
				MatSetValue(_A, idCells[1], idCells[0],  (abs( _Velocity(j) )+ _c ) * Fj.getMeasure()/2.0, ADD_VALUES );  

				/*************** Momentum conservation equation *****************/
				if (_nbTimeStep == 0){
					MatSetValue(_DivTranspose, j, idCells[0], getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES );
					MatSetValue(_DivTranspose, j, idCells[1], getOrientation(j,Ctemp2) * Fj.getMeasure(), ADD_VALUES ); 
				}
				/*************(-1) x CONVECTION (is not computed in 3 space dimensions)********************/ 
				Convection= 0; 
				std::vector<Cell> Support_j;
				Support_j.push_back(_mesh.getCell(idCells[0]));
				Support_j.push_back(_mesh.getCell(idCells[1]));

				for (int nei =0; nei <2; nei ++){ // we are in the case where an interior face has two neighbours
					Cell K = _mesh.getCell(idCells[nei]); 
					std::vector<int> idFaces = K.getFacesId();
					VecGetValues(_primitiveVars,1,&idCells[nei],&rho);
					//If  the face is periodic and K isn't the direct neighbour of Fj, recover the geometric informations of the associated periodic cell 
					Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[nei] ? _mesh.getFace(_FacePeriodicMap.find(j)->second ):  Fj;
					
					for (int f =0; f <K.getNumberOfFaces(); f ++){
						bool IsfInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),idFaces[f] ) != _InteriorFaceSet.end() ;
						bool IsfWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),idFaces[f] ) != _WallBoundFaceSet.end() ;
						bool IsfSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),idFaces[f] ) != _SteggerBoundFaceSet.end();
						std::map<int,int>::iterator it = _FacePeriodicMap.begin();
						while ( ( idFaces[f] !=it->second) && (it !=_FacePeriodicMap.end() ) )it++;

						Face Facef = _mesh.getFace( idFaces[f] );
						std::vector< int> idCellsOfFacef =  Facef.getCellsId();
						std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
						I = _Nmailles + idFaces[f];
						VecGetValues(_primitiveVars,1,&I	,&q);
						double gradiv = -( _c + _uMax)/2.0 * Fj_physical.getMeasure() * getOrientation(j,K) *Facef.getMeasure()* getOrientation(idFaces[f], K) /( (_Ndim==2 )? _perimeters[idCells[nei]] : 1.0);
						MatSetValue(_A, IndexFace, I, gradiv, ADD_VALUES ); 

						//************* _Ndim-dimensional terms *************//
						std::vector<Node> K_Nodes;
						for (int i=0; i < K.getNodesId().size(); i++)
							K_Nodes.push_back(_mesh.getNode(K.getNodesId()[i]) );
						
						std::vector< Point > IntegrationNodes(Facef.getNumberOfNodes());
						std::vector< double > Weights(Facef.getNumberOfNodes());
						for (int i =0; i <Facef.getNumberOfNodes(); i++){
							IntegrationNodes[i] = _mesh.getNode( Facef.getNodesId()[i] ).getPoint();
							//In quads it's 1/8 (instead of 1/4 given by trapezoid formula given on reference elem) since the loop on the face then on the faces'nodes implies that the integral is computed twice 
							//TODO -> to improve & what about triangles
							Weights[i] = 1.0/(K.getNumberOfFaces() * _Ndim ) ; 
						}

						for (int inteNode=0; inteNode <IntegrationNodes.size(); inteNode ++){
							std::vector<double> MomentumRT_in_Xf = MomentumRaviartThomas_at_point_X(K, idCells[nei], IntegrationNodes[inteNode] ); 
							std::vector<double> qtensorielq = TensorProduct(MomentumRT_in_Xf, MomentumRT_in_Xf);
							std::vector<double> GradientPsi_j_in_Xf = Gradient_PhysicalBasisFunctionRaviartThomas(K, idCells[nei], Support_j, Fj_physical,j, IntegrationNodes[inteNode]  ); 
							Convection +=  Weights[inteNode] * fabs( det( JacobianTransfor_K_X( xToxhat(K, IntegrationNodes[inteNode]  , K_Nodes), K_Nodes) ) ) * 1.0/rho * Contraction(qtensorielq, GradientPsi_j_in_Xf); 	
						}
							
						//************* (_Ndim-1)-dimensional terms *************//
						if (IsfInterior){											
							if ( _FacePeriodicMap.find(idFaces[f]) != _FacePeriodicMap.end()  )
								idCellsOfFacef.push_back( _mesh.getFace( _FacePeriodicMap.find(idFaces[f])->second ).getCellsId()[0]  );
							VecGetValues(_primitiveVars,1,&idCellsOfFacef[1],&rhoExt);	
						}
						else if (it != _FacePeriodicMap.end())
							idCellsOfFacef.push_back( _mesh.getFace( it->first ).getCellsId()[0]  );
						else if (IsWallBound) rhoExt = rhoInt;
						else if (IsSteggerBound) rhoExt = getboundaryPressure().find(j)->second;
						
						for (int inteNode=0; inteNode <IntegrationNodes.size(); inteNode ++){
							std::vector<double> jumpPsi(_Ndim, 0.0);
							std::vector<double> meanRhoU(_Ndim, 0.0);
							
							for (int ncell =0; ncell < idCellsOfFacef.size(); ncell ++){
								Cell Neibourg_of_f = _mesh.getCell(idCellsOfFacef[ncell]);
								Face Facef_physical = (_FacePeriodicMap.find(idFaces[f]) != _FacePeriodicMap.end()) && idCellsOfFacef[ncell] != idCells[nei] ?  _mesh.getFace(_FacePeriodicMap.find(idFaces[f])->second ):  Facef;
								Face Facej_physical_2 = (_FacePeriodicMap.find(j) != _FacePeriodicMap.end()) && idCellsOfFacef[ncell] != idCells[0] ?  _mesh.getFace(_FacePeriodicMap.find(j)->second ):  Fj_physical;
								std::vector<double> Psi_j_in_Xf = PhysicalBasisFunctionRaviartThomas(Neibourg_of_f, idCellsOfFacef[ncell], Support_j, Facej_physical_2,j, IntegrationNodes[inteNode] ); 
								std::vector<double> velocityRT_in_Xf = MomentumRaviartThomas_at_point_X(Neibourg_of_f,idCellsOfFacef[ncell], IntegrationNodes[inteNode]  ); 
								for (int ndim =0; ndim < _Ndim; ndim ++ ){
									jumpPsi[ndim] += getOrientation( idFaces[f],Neibourg_of_f) *  Psi_j_in_Xf[ndim];  
									meanRhoU[ndim] += velocityRT_in_Xf[ndim]/2.0; //TODO what about outside ? 
								}
							}
							double dotprod =0;
							//The integral on the face j is computed twice because of choice of implementation, so divide by 2 to recover consistency
							double factor = ((idFaces[f] == j) || it->first == j ) ? 1.0/2.0 : 1.0;	
							for (int ndim =0; ndim < _Ndim; ndim ++ )
								dotprod  += 1.0/IntegrationNodes.size() * q * Facef.getMeasure() * jumpPsi[ndim] * meanRhoU[ndim] * 2.0/(rho +rhoExt) * factor ; 	
							Convection -= dotprod; 		
						}
					}
				} 
				VecSetValue(_Conv, IndexFace, Convection, ADD_VALUES );

				// *************** Grad^perp tilde (Grad^per)^*  ************//
				/* if (_Ndim ==2 ){
					double Inv_Dsigma, length_sigma_perp,  Inv_Df, length_f_perp;
					MatGetValues(_InvVol, 1, &IndexFace, 1, &IndexFace, &Inv_Dsigma);
					length_sigma_perp = 1/(Inv_Dsigma * Fj.getMeasure());
					//TODO periodic
					std::vector<int> idNodesSigma = Fj.getNodesId();
					for (int n=0; n < idNodesSigma.size() ; n++){
						std::vector<int> idFacesn = _mesh.getNode( idNodesSigma[n] ).getFacesId();
						double measureBoundaryDualVolumOfN =0;
						for (int f=0; f< idFacesn.size(); f++){
							Face FaceOfNodeN =  _mesh.getFace( idFacesn[f] );
							for (int cell=0; cell < FaceOfNodeN.getNumberOfCells(); cell++){
								Cell K = _mesh.getCell( FaceOfNodeN.getCellsId()[cell] );
								measureBoundaryDualVolumOfN += sqrt( (K.x() -  FaceOfNodeN.x())*(K.x() -  FaceOfNodeN.x()) + (K.y() -  FaceOfNodeN.y())*(K.y() -  FaceOfNodeN.y()));
							}
						}
						for (int f=0; f< idFacesn.size(); f++){
							int I_f = _Nmailles + idFacesn[f];
							MatGetValues(_InvVol, 1, &I_f, 1, &I_f, &Inv_Df);
							length_f_perp = 1/(Inv_Df * _mesh.getFace(idFacesn[f]).getMeasure());
							double rotrot = _rhoMax * _uMax/2.0 *length_sigma_perp * getOrientationNode(idNodesSigma[n],j) * length_f_perp * getOrientationNode(idNodesSigma[n], idFacesn[f])/measureBoundaryDualVolumOfN; 
							MatSetValue(_A, IndexFace, I_f, rotrot, ADD_VALUES ); 

						}
					}
				} */
			}
			else if (IsWallBound || IsSteggerBound ) { 
				/****************** Density conservation equation *********************/
				VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);		
				MatSetValue(_A, idCells[0], IndexFace, -getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES ); 
				if (IsSteggerBound){
					MatSetValue(_A, idCells[0], idCells[0], -( abs(_Velocity(j) )+_c ) * Fj.getMeasure()/2.0, ADD_VALUES );
					VecSetValue(_BoundaryTerms, idCells[0], ( abs(_Velocity(j) )+_c ) * Fj.getMeasure()/2.0 * getboundaryPressure().find(j)->second, ADD_VALUES );
				} 
			}	
		}

		VecAssemblyBegin(_BoundaryTerms);
		VecAssemblyEnd(_BoundaryTerms);
		MatAssemblyBegin(_DivTranspose, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_DivTranspose, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_Conv);
		VecAssemblyEnd(_Conv);	
		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

		Vec Au;
		VecDuplicate(_b, &Au);
		MatMult(_A,_primitiveVars, Au);
		MatMult(_InvVol,Au, _b); 

		/***********Adding boundary terms to density equation, convection term and pressure gradient to momentum equation **************/
		//TODO : créer vecteur dans intialize() et les supprimer dans terminate()
		Vec Temporary1, Temporary2, Pressure, GradPressure;
		VecCreate(PETSC_COMM_SELF, & Temporary2);
		VecDuplicate(_BoundaryTerms, & Temporary1);
		VecSetSizes(Temporary2, PETSC_DECIDE, _Nfaces);
		VecSetFromOptions(Temporary2);
		//add boundary terms : AU^n + V^{-1}Boundterms //
		MatMult(_InvVol, _BoundaryTerms, Temporary1); 
		VecAXPY(_b,     1, Temporary1); 
		//add convection terms  : AU^n + V^{-1}Boundterms - V^{-1}_Conv//
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
		double *Product = new double[_Nfaces]; //TODO  can we avoid creating a pointer only to insert in the bigger Vector ?
		PetscScalar gradp;
		for (int i=0; i <_Nfaces; i++){
			VecGetValues(Temporary2,1,&i,&gradp);
			Product[i] = gradp; 
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
		VecDestroy(& Au); 
		VecDestroy(& GradPressure);
		VecDestroy(& Pressure);	
		delete[] indices2;
		delete[] Product;
		
	}
	if (_nbTimeStep == 0)
		ComputeMinCellMaxPerim(); 
	double numax = _Ndim*2;	
	
	double dt = _cfl * _minCell / (_maxPerim * (_uMax + _c) * 2 ) ;
	double PreviousTime = _Time.back();
	_Time.push_back(PreviousTime+ dt);
	
	return dt; 
}


//*********************Gradient Raviart-Thomas********************//
std::vector<double> EulerBarotropicStaggered::HessienneTransfo(const int component, const std::vector<Node> &K_Nodes){
	std::vector<double> HessienneTransfo(_Ndim * _Ndim, 0.0);
	if (_Ndim ==2){ 
		if (K_Nodes.size()== 4){
			if (component == 0){
				HessienneTransfo[1] =  K_Nodes[3].x()  - K_Nodes[2].x() - K_Nodes[0].x()  + K_Nodes[1].x() ;
				HessienneTransfo[2] =  K_Nodes[3].x()  - K_Nodes[2].x() - K_Nodes[0].x()  + K_Nodes[1].x() ;
			}
			else if (component == 1){
				HessienneTransfo[1] = K_Nodes[3].y()  - K_Nodes[2].y() - K_Nodes[0].y()  + K_Nodes[1].y() ;
				HessienneTransfo[2] = K_Nodes[3].y()  - K_Nodes[2].y() - K_Nodes[0].y()  + K_Nodes[1].y() ;
			}
					
		}
	}
	return HessienneTransfo;
}

std::vector<double>  EulerBarotropicStaggered::Gradient_ReferenceBasisFunctionRaviartThomas(int i, const std::vector<Node> &K_Nodes, const Point & Xhat ){
	std::vector<double>  Gradient_Psihat(_Ndim*_Ndim, 0.0); 
	if (_Ndim==1)  Gradient_Psihat[0] = 1;
	else if (_Ndim ==2){
		if (K_Nodes.size() == 4){
			if (i ==0 || i ==2 ) Gradient_Psihat[0] = 1;
			else if (i ==1 || i ==3 ) Gradient_Psihat[3] = 1;
			double areaKleft = fabs( (K_Nodes[0].x()-K_Nodes[1].x())* (K_Nodes[2].y()-K_Nodes[1].y()) - (K_Nodes[0].y()-K_Nodes[1].y())* (K_Nodes[2].x()-K_Nodes[1].x())  )/2.0;
			double areaKright = fabs( (K_Nodes[2].x()-K_Nodes[3].x())* (K_Nodes[0].y()-K_Nodes[3].y()) - (K_Nodes[2].y()-K_Nodes[3].y())* (K_Nodes[0].x()-K_Nodes[3].x())  )/2.0;
			double distortedQuadsX = ((K_Nodes[0].x()-K_Nodes[1].x()) * (K_Nodes[3].y()-K_Nodes[2].y()) - (K_Nodes[0].y()-K_Nodes[1].y())* ( K_Nodes[3].x()-K_Nodes[2].x() ) )/(2.0 *(areaKleft+ areaKright) );
			double distortedQuadsY = ((K_Nodes[1].x()-K_Nodes[2].x()) * (K_Nodes[3].y()-K_Nodes[0].y()) - (K_Nodes[1].y()-K_Nodes[2].y())* ( K_Nodes[3].x()-K_Nodes[0].x() ) )/(2.0 *(areaKleft+ areaKright) );
			Gradient_Psihat[0] += distortedQuadsX * (2 * Xhat.x() - 1 );
			Gradient_Psihat[3] += distortedQuadsY * (2 * Xhat.y() - 1 );
		}
		else if (K_Nodes.size() == 3){
			if (i ==0 || i ==2 ){
				Gradient_Psihat[0] =1;
				Gradient_Psihat[3] =1;
			}
			else if (i==1) {
				Gradient_Psihat[0] =1;
				Gradient_Psihat[3] =1;
			}
		}
	}
	return Gradient_Psihat;	
}


std::vector<double> EulerBarotropicStaggered::Gradient_PhysicalBasisFunctionRaviartThomas(Cell K, int idcell, std::vector<Cell> Support, Face Facej, int j, Point X){
	std::vector<double> xf;
	xf.push_back(X.x());
	if (_Ndim ==2) xf.push_back(X.y());
	std::vector<double> Gradient_PhysicalPsif_in_X(_Ndim* _Ndim, 0.0);
	if (_nbTimeStep == 0 ){ 
		std::vector<Node> K_Nodes;
		for (int i=0; i < K.getNodesId().size(); i++)
			K_Nodes.push_back(_mesh.getNode(K.getNodesId()[i]) );
		std::vector<double>  JacobianTransfor_K_Xhat = JacobianTransfor_K_X( xToxhat(K, X, K_Nodes), K_Nodes ); 
		double J = fabs( det(JacobianTransfor_K_Xhat) );
		std::vector<double>  JacobianTransfor_K_Xhat_Inversed = Inverse( JacobianTransfor_K_Xhat );
		
		std::vector<double> derivTrace(_Ndim, 0.0);
		for (int edim =0; edim <_Ndim; edim ++){
			for (int m =0; m<_Ndim; m++){
				for (int u =0; u<_Ndim; u++){
					for (int l =0; l<_Ndim; l++)
						derivTrace[edim] +=JacobianTransfor_K_Xhat_Inversed[m*_Ndim + u ]* JacobianTransfor_K_Xhat_Inversed[l*_Ndim + edim ]*HessienneTransfo(u,K_Nodes)[l*_Ndim + m];
				}
			}
		}
	
		bool K_is_in_Support = false;
		for (const auto &cell: Support){
			if (K.getBarryCenter().x() == cell.getBarryCenter().x() && K.getBarryCenter().y() == cell.getBarryCenter().y())
				K_is_in_Support =true;
		}

		for (int  m = 0; m < K.getNumberOfFaces(); m ++){
			if (FindlocalBasis(m, Facej, j, K,  K_Nodes) && K_is_in_Support == true ){
				std::vector<double>  Gradient_ReferencePsif_in_X = Gradient_ReferenceBasisFunctionRaviartThomas(m, K_Nodes, xToxhat(K, X, K_Nodes));
				for (int i =0; i < _Ndim ; i++){
					for (int k =0; k < _Ndim ; k++){
						for (int l =0; l < _Ndim ; l++){
							for (int t =0; t < _Ndim ; t++){
								Gradient_PhysicalPsif_in_X[i*_Ndim + k] +=  getOrientation(j,K) * Facej.getMeasure()*1.0/fabs(J) *(
									JacobianTransfor_K_Xhat[i*_Ndim  +  l] * Gradient_ReferencePsif_in_X[l*_Ndim + t] * JacobianTransfor_K_Xhat_Inversed[t*_Ndim + k] 
									+ JacobianTransfor_K_Xhat_Inversed[l *_Ndim + k ] *  HessienneTransfo(i,K_Nodes)[l*_Ndim +t] * ReferenceBasisFunctionRaviartThomas(m, xToxhat(K, X, K_Nodes),  K_Nodes )[t] ); 
							}
							Gradient_PhysicalPsif_in_X[i*_Ndim + k] +=   -getOrientation(j,K) * Facej.getMeasure()*1.0/fabs(J) *( 
								JacobianTransfor_K_Xhat[i *_Ndim + l]*ReferenceBasisFunctionRaviartThomas(m, xToxhat(K, X, K_Nodes),  K_Nodes )[l] * derivTrace[k] );
						}
					}
				}		
			}
		}
		bool already_there =false;
		if (_GradientPhysicalPsif.find(idcell) != _GradientPhysicalPsif.end() ){
			if (_GradientPhysicalPsif.find(idcell)->second.find(j) != _GradientPhysicalPsif.find(idcell)->second.end() ){
				for (const auto & e : _GradientPhysicalPsif.find(idcell)->second.find(j)->second ){
					if ( (_Ndim==1 && e.first[0] == X.x()) || (_Ndim ==2  && e.first[0] == X.x() && e.first[1] == X.y() )){
						already_there = true;
					}
				}	
			}
		}
		if (already_there == false)
			_GradientPhysicalPsif[idcell][j].push_back( make_pair(xf, Gradient_PhysicalPsif_in_X) );
		
	}
	else{
		for (const auto & e : _GradientPhysicalPsif.find(idcell)->second.find(j)->second ){
			if ( (_Ndim==1 && e.first[0] == X.x()) || (_Ndim ==2  && e.first[0] == X.x() && e.first[1] == X.y()) ){
				for (int k= 0; k < Gradient_PhysicalPsif_in_X.size(); k++){
					Gradient_PhysicalPsif_in_X[k]= e.second[k];
				}
			}

		}
	}
	return Gradient_PhysicalPsif_in_X;
}


std::vector<double> EulerBarotropicStaggered::MomentumRaviartThomas_at_point_X(Cell K,int idcell, Point X){
	std::vector<double> VelocityRT(_Ndim,0.0);
	std::vector< int > idFaces = K.getFacesId();
	std::vector<Cell> Support ;
	Support.push_back(K);
	PetscInt IndexFace;
	PetscScalar q;
	for (int f=0; f< K.getNumberOfFaces(); f++){
		PetscInt IndexFace = _Nmailles + idFaces[f];
		std::vector<double> Psif = PhysicalBasisFunctionRaviartThomas(K, idcell, Support, _mesh.getFace(idFaces[f]), idFaces[f], X);
		VecGetValues(_primitiveVars,1,&IndexFace,&q);
		for (int k=0; k< _Ndim; k++)  VelocityRT[k] += q * Psif[k]; 
	}
	return VelocityRT;
 }

 std::vector<double> EulerBarotropicStaggered::TensorProduct(std::vector<double> &u, std::vector<double> &v){
	std::vector<double> tensorproduct(u.size() * v.size(),0.0);
	for (int i=0; i< _Ndim; i++){
		for (int j=0; j< _Ndim; j++) tensorproduct[i*_Ndim + j] = u[i] * v[j];
	}
	return tensorproduct;
 }

double EulerBarotropicStaggered::Contraction(std::vector<double> &u, std::vector<double> &v){
	double contraction =0;
	for (int i=0; i< _Ndim; i++){
		for (int j=0; j< _Ndim; j++) contraction += u[i*_Ndim + j]* v[i*_Ndim + j];
	}
	return contraction;
 }

 std::vector<double> EulerBarotropicStaggered::Inverse(std::vector<double> &u){
	std::vector<double> Inverse( u.size() );
	if (_Ndim ==1){
		Inverse[0] = 1.0/u[0]; 
	}
	else if (_Ndim ==2 ){
		double detu = (u[0] * u[3] - u[2] * u[1]); 
		Inverse[0] = u[3]/detu; 
		Inverse[1] = -u[1]/detu;
		Inverse[2] = -u[2]/detu;
		Inverse[3] = u[0]/detu;
	}
	return Inverse;
 }

void EulerBarotropicStaggered::Rhomax_Umax_Cmax(){
	/************ Max rho, Max u *******************/
	double rho;
	PetscInt zero=0;
	VecGetValues(_primitiveVars,1,&zero,&rho);
	_uMax = _Velocity(0);
	_rhoMax = rho;
	for (int n=0; n <_Nmailles; n++){
		VecGetValues(_primitiveVars,1,&n,&rho);
		_rhoMax = max(_rhoMax, rho);
	}
	for (int n=0; n <_Nfaces; n++){
		PetscInt I = _Nmailles + n;
		_uMax = max(_uMax, abs(_Velocity(n)));
	}
	_c =  _compressibleFluid->vitesseSon(_rhoMax); 	 
}


void EulerBarotropicStaggered::computeNewtonVariation(){
	if(_timeScheme == Explicit){	
		VecCopy(_b,_newtonVariation);
		VecScale(_newtonVariation, _dt);
	}
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

	PetscScalar q, rho_sigma, u;
	VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk + deltaV 
	UpdateDualDensity(); // \rho^{n+1}_K -> \rho^{n+1}_\sigma 

	for (int f=0; f< _Nfaces; f++){
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),f ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),f ) != _SteggerBoundFaceSet.end() ;	
		PetscInt I = _Nmailles +f;
		if (IsWallBound) VecSetValue(_primitiveVars, I, 0, INSERT_VALUES); 
		else if (IsSteggerBound){
			VecGetValues(_DualDensity, 1,&f, &rho_sigma);
			VecSetValue(_primitiveVars, I, _Velocity(f) * rho_sigma , INSERT_VALUES); 
		}
	}
	return converged;
}

void EulerBarotropicStaggered::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/");///Results
	prim+=_fileName;

	if(_mpi_size>1){
        VecScatterBegin(_scat,_primitiveVars,_primitiveVars_seq,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(  _scat,_primitiveVars,_primitiveVars_seq,INSERT_VALUES,SCATTER_FORWARD);
    }
	if(_savePressure){
		for (int i = 0 ; i < _Nmailles  ; i++){
				if (_mpi_size > 1)
					VecGetValues(_primitiveVars_seq,1,&i,&_Pressure(i));
				else 
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
			std::map<int,int>::iterator it2 = _FacePeriodicMap.begin();
			while ( ( i !=it2->second) && (it2 != _FacePeriodicMap.end() ) )
				it2++;
			periodicFaceNotComputed = (it2 !=  _FacePeriodicMap.end());
			int j = (periodicFaceNotComputed ==true) && (_indexFacePeriodicSet == true) ? it2->first : i; // in periodic k stays i, if it has been computed by scheme and takes the value of its matched face 
			int I= _Nmailles + j;
			double rho_sigma, q;
			VecGetValues(_DualDensity, 1,&j, &rho_sigma);

			if (_mpi_size > 1)
				VecGetValues(_primitiveVars_seq,1,&I,&q);
			else
				VecGetValues(_primitiveVars,1,&I,&q);
			_Velocity(j) = q / rho_sigma;

			Face Fj = _mesh.getFace(i);
			std::vector< int > idCells = Fj.getCellsId();
			Cell Ctemp1 = _mesh.getCell(idCells[0]);
			double orien1 = getOrientation(i,Ctemp1);
			std::vector<double> M1(_Ndim), M2(_Ndim);
			Point xK = Ctemp1.getBarryCenter();
			Point xsigma = Fj.getBarryCenter();

			M1[0] =  orien1 * Fj.getMeasure()*(xsigma.x()- xK.x());
			if (_Ndim >1) M1[1] = orien1 * Fj.getMeasure()*(xsigma.y()- xK.y());
			
			if (Fj.getNumberOfCells() == 2){
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
				double orien2 = getOrientation(i,Ctemp2);
				Point xK = Ctemp2.getBarryCenter();

				M2[0] = orien2 * Fj.getMeasure()*(xsigma.x()- xK.x());
				if (_Ndim >1) M2[1] = orien2 * Fj.getMeasure()*(xsigma.y()- xK.y());
			
				for (int k=0; k< _Ndim; k++){
					_Velocity_at_Cells(idCells[0], k) += _Velocity(i) * M1[k]/Ctemp1.getMeasure(); 
					_Velocity_at_Cells(idCells[1], k) += _Velocity(i) * M2[k]/Ctemp2.getMeasure(); 
				}
				_DivVelocity( idCells[0]) += orien1 * Fj.getMeasure() * _Velocity(i)/Ctemp1.getMeasure();
				_DivVelocity( idCells[1]) += orien2 * Fj.getMeasure() * _Velocity(i)/Ctemp2.getMeasure();

			}
			else if  (Fj.getNumberOfCells() == 1){
				for (int k=0; k< _Ndim; k++){
					_Velocity_at_Cells(idCells[0], k) += _Velocity(i) * M1[k]/Ctemp1.getMeasure(); 
				}
				_DivVelocity( idCells[0]) += orien1 * Fj.getMeasure() * _Velocity(i)/Ctemp1.getMeasure();
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
			_Velocity_at_Cells.writeVTK(prim+"_VelocityAtCells");
			_Velocity.writeVTK(prim+"_Velocity");
			_DivVelocity.writeVTK(prim+"DivVelocity");
			break;
		case MED :
			_Velocity.writeMED(prim+"_Velocity");
			break;
		case CSV :
			_Velocity.writeCSV(prim+"_Velocity");
			_Velocity_at_Cells.writeCSV(prim+"_VelocityAtCells");
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

void EulerBarotropicStaggered::terminate(){ 
	VecDestroy(& _Conv); 
	VecDestroy(& _DualDensity);
	delete _compressibleFluid;
	_compressibleFluid = nullptr;
	delete[]_vec_normal;
	VecDestroy(&_newtonVariation);
	VecDestroy(&_primitiveVars);
	VecDestroy(&_b);
	MatDestroy(& _A); 

	MatDestroy(&_InvVol); 
	MatDestroy(& _InvSurface);

	MatDestroy(& _Div); 
	MatDestroy(&_LaplacianPressure);
	VecDestroy(& _BoundaryTerms);	

	MatDestroy(& _DivTranspose); 
	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);

	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq);
	
}


