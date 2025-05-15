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

	// matrix LAPLACIAN Pressure (without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianPressure); 
	MatSetSizes(_LaplacianPressure, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles ); 
	MatSetFromOptions(_LaplacianPressure);
	MatSetUp(_LaplacianPressure);

	// *************** Time dependent matrices ************ //
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

	AssembleMetricsMatrices();

	/* if(_mpi_size>1 && _mpi_rank == 0)
    	VecCreateSeq(PETSC_COMM_SELF, _globalNbUnknowns, &_primitiveVars_seq);//For saving results on proc 0
    VecScatterCreateToZero(_primitiveVars,&_scat,&_primitiveVars_seq); */

	createKSP();
	PetscPrintf(PETSC_COMM_WORLD,"SOLVERLAB Newton solver ");
	*_runLogFile << "SOLVERLAB Newton solver" << endl;
	_runLogFile->close();

	_initializedMemory=true;
	
	
	
	VecCreate(PETSC_COMM_SELF, & _DualDensity);//Current primitive variables at Newton iteration k between time steps n and n+1
	VecSetSizes(_DualDensity, PETSC_DECIDE, _Nfaces);
	VecSetFromOptions(_DualDensity);
	VecZeroEntries(_DualDensity);

	UpdateDualDensity();	

	if (_timeScheme == Implicit) {
		MatCreate(PETSC_COMM_SELF, & _JacobianMatrix); 
		MatSetSizes(_JacobianMatrix, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
		MatSetFromOptions(_JacobianMatrix);
		MatSetUp(_JacobianMatrix);
	}

	//************ Time independent matrices ********** //

	// matrix minus GRADIENT (we will impose to be 0 on faces so that u^n+1 = u^n at the _Boundary)
	// Vector CONVECTION 
	VecCreate(PETSC_COMM_SELF, & _Conv); 
	VecSetSizes(_Conv, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_Conv);
	VecSetUp(_Conv);
	// Vector Pressure gradient
	VecCreate(PETSC_COMM_SELF, & _GradPressure); 
	VecSetSizes(_GradPressure, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_GradPressure);
	VecSetUp(_GradPressure);
	VecZeroEntries(_GradPressure);

	save();//save initial data

}


std::vector<double>  EulerBarotropicStaggered::getTimeEvol(){
	return _Time;
}

void  EulerBarotropicStaggered::setboundaryVelocityVector(int j,  std::vector<double>  boundaryVelocityVector){
	for (int idim = 0; idim < _Ndim; ++idim)
		_boundaryVelocityVector[j].push_back(boundaryVelocityVector[idim]);
}
std::map<int,std::vector<double> >  EulerBarotropicStaggered::getboundaryVelocityVector() const {
		return _boundaryVelocityVector;
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
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  )
				VecSetValue(_DualDensity, _FacePeriodicMap.find(j)->second, (rhoL * D_sigmaL + rhoR*D_sigmaR)/(D_sigmaL +  D_sigmaR), INSERT_VALUES );
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
	VecZeroEntries(_b);
	VecZeroEntries(_BoundaryTerms);
	VecZeroEntries(_GradPressure);
	MatZeroEntries(_A);
	if (_timeScheme == Implicit)
		MatZeroEntries(_JacobianMatrix);
	PetscScalar rho,p, u, Convection, rhoInt, rhoExt, q, WaveVelocity;
	
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
		
		//****************** Convection *****************// 
		Convection= 0; 
		std::vector<Cell> Support_j;
		Support_j.push_back(_mesh.getCell(idCells[0]));
		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  ) 
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j) ->second).getCellsId()[0]  );
			Support_j.push_back(_mesh.getCell(idCells[1]));
		}
		
		for (int nei =0; nei <Support_j.size(); nei ++){ 
			Cell K = _mesh.getCell(idCells[nei]); 
			std::vector<int> idFaces = K.getFacesId();
			VecGetValues(_primitiveVars,1,&idCells[nei],&rho);
			//If  the face is periodic and K isn't the direct neighbour of Fj, recover the geometric informations of the associated periodic cell 
			Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[nei] ? _mesh.getFace(_FacePeriodicMap.find(j)->second ):  Fj;
			
			for (int f =0; f <K.getNumberOfFaces(); f ++){
				std::map<int,int>::iterator it = _FacePeriodicMap.begin();
				while ( ( idFaces[f] !=it->second) && (it !=_FacePeriodicMap.end() ) )it++;
				bool IsfPeriodicTwin = (it != _FacePeriodicMap.end());
				bool IsfInterior     = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),idFaces[f] ) != _InteriorFaceSet.end() ;
				bool IsfSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),idFaces[f] ) != _SteggerBoundFaceSet.end();
				bool IsfWall         = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),idFaces[f] ) != _WallBoundFaceSet.end();
				
				Face Facef = _mesh.getFace( idFaces[f] );
				std::vector< int> idCellsOfFacef =  Facef.getCellsId();
				std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
				I = _Nmailles + idFaces[f];
				VecGetValues(_primitiveVars,1,&I	,&q);
				WaveVelocity = (_timeScheme == Implicit ) ? _uMax/2.0 : ( _c + _uMax)/2.0 ;
				// gradDiv //
				double gradiv = - WaveVelocity * Fj_physical.getMeasure() * getOrientation(j,K) *Facef.getMeasure()* getOrientation(idFaces[f], K) /( (_Ndim==2 )? _perimeters[idCells[nei]] : 1.0);
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
				cout <<j<< " Convection = " << Convection<< endl;
 				
				//************* (_Ndim-1)-dimensional terms *************//
				std::vector<double> rhoMean;
				rhoMean.push_back( rho );
				if (IsfInterior){											
					if ( _FacePeriodicMap.find(idFaces[f]) != _FacePeriodicMap.end()  )
						idCellsOfFacef.push_back( _mesh.getFace( _FacePeriodicMap.find(idFaces[f])->second ).getCellsId()[0]  );
					VecGetValues(_primitiveVars,1,&idCellsOfFacef[1],&rhoExt);	
					rhoMean.push_back( rhoExt ) ;
				}
				if (IsfPeriodicTwin){   
					idCellsOfFacef.push_back( _mesh.getFace( it->first ).getCellsId()[0]  ); 
					VecGetValues(_primitiveVars,1,&idCellsOfFacef[1],&rhoExt);	
					rhoMean.push_back( rhoExt );
				}
				if (IsfInterior || IsfPeriodicTwin){
					for (int inteNode=0; inteNode <IntegrationNodes.size(); inteNode ++){
						std::vector<double> jumpPsi(_Ndim, 0.0);
						std::vector<double> meanRhoU(_Ndim, 0.0);
						for (int ncell =0; ncell < idCellsOfFacef.size(); ncell ++){
							Cell Neibourg_of_f    = _mesh.getCell(idCellsOfFacef[ncell]);
							Face Facef_physical   =   (_FacePeriodicMap.find(idFaces[f]) != _FacePeriodicMap.end()) && idCellsOfFacef[ncell] != idCells[nei] ?  _mesh.getFace(_FacePeriodicMap.find(idFaces[f])->second ):  Facef;
							Face Facej_physical_on_K_ncell =   (_FacePeriodicMap.find(j) != _FacePeriodicMap.end()) && idCellsOfFacef[ncell] != idCells[nei] ?  _mesh.getFace(_FacePeriodicMap.find(j)->second )         :  Fj_physical;
							for (int i =0; i <Facef.getNumberOfNodes(); i++)
								IntegrationNodes[i] = _mesh.getNode( Facef_physical.getNodesId()[i] ).getPoint();
							std::vector<double> Psi_j_in_Xf = PhysicalBasisFunctionRaviartThomas(Neibourg_of_f, idCellsOfFacef[ncell], Support_j, Facej_physical_on_K_ncell,j, IntegrationNodes[inteNode] ); 
							std::vector<double> MomenentumRT_in_Xf = MomentumRaviartThomas_at_point_X(Neibourg_of_f,idCellsOfFacef[ncell], IntegrationNodes[inteNode]  ); 
							for (int ndim =0; ndim < _Ndim; ndim ++ ){
								jumpPsi[ndim] += Psi_j_in_Xf[ndim]*getOrientation( idFaces[f],Neibourg_of_f)  ; //* ( ((it != _FacePeriodicMap.end()) && it->first == j) || IsfInterior ? getOrientation( idFaces[f],Neibourg_of_f) : 1.0 ) ;  
								meanRhoU[ndim] += MomenentumRT_in_Xf[ndim]/idCellsOfFacef.size() * 1.0/rhoMean[ncell];
							}
						}
						double dotprod =0;
						//The integral on the face j is computed twice because of choice of implementation, so divide by 2 to recover consistency
						for (int ndim =0; ndim < _Ndim; ndim ++ )
							dotprod  += jumpPsi[ndim] * meanRhoU[ndim] * ( ((idFaces[f] == j) || ((it != _FacePeriodicMap.end()) && it->first == j) ) ? 1.0/2.0 : 1.0 ) ; 	
						Convection -= Facef.getMeasure() * q* getOrientation(idFaces[f], _mesh.getCell( idCellsOfFacef[0]) ) * dotprod /IntegrationNodes.size(); 	
					}
				}
				else if (IsfWall && (idFaces[f] == j ))
					Convection -= _compressibleFluid->vitesseSon(rho) * Facef.getMeasure() * q;
				else if (IsfSteggerBound){
					Cell interiorCell = _mesh.getCell( idCellsOfFacef[0]);
					double u_b = getboundaryVelocity().find(idFaces[f])->second * getOrientation(idFaces[f], interiorCell ); //TODO idfaces f and its interior cell
					double c_b = _compressibleFluid->vitesseSon(getboundaryPressure().find(idFaces[f])->second);
					double lambdaPlus  =  u_b + c_b;
					double lambdaMinus =  u_b - c_b;
					VecGetValues(_primitiveVars,1,&idCellsOfFacef[0],&rhoInt);	
					VecGetValues(_primitiveVars,1,&I,&q);	
					std::vector<double> lambdaPlusVector(_Ndim), lambdaMinusVector(_Ndim);
					for (int ndim =0; ndim < _Ndim; ndim ++ ){
						lambdaPlusVector[ndim] =  getboundaryVelocityVector().find(idFaces[f])->second[ndim]  + c_b * _vec_sigma.find(idFaces[f])->second[ndim]  * getOrientation(idFaces[f], interiorCell);
						lambdaMinusVector[ndim] = getboundaryVelocityVector().find(idFaces[f])->second[ndim]  - c_b * _vec_sigma.find(idFaces[f])->second[ndim]  * getOrientation(idFaces[f], interiorCell) ;
					}
					std::vector<double> Psi_j_in_Xf = PhysicalBasisFunctionRaviartThomas(interiorCell, idCellsOfFacef[0], Support_j, Fj ,j, Facef.getBarryCenter() );
					double dotprod = 0;
					for (int ndim =0; ndim < _Ndim; ndim ++ ){
						dotprod += ( ( (1 - u_b/c_b)* rhoInt +  getOrientation(idFaces[f], interiorCell ) * q/c_b )* lambdaPlusVector[ndim] * lambdaPlus/2.0
									+  getboundaryPressure().find(idFaces[f])->second * lambdaMinusVector[ndim] * lambdaMinus /2.0  ) * Psi_j_in_Xf[ndim];
					}
					Convection -= Facef.getMeasure() * dotprod ; 
				}	
			}
		} 
		cout <<j<< "Fj.x() = "<< Fj.x() <<" Convection before adding it = " << Convection<< endl;
		VecSetValue(_Conv, IndexFace, Convection, ADD_VALUES );

		// Density equation 
		if (IsInterior){
			/* if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  ) 
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j) ->second).getCellsId()[0]  ); */
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
	
			// -DivRhoU_{idcell[1], j}, -DivRhoU_{idcell[0], j}, (IndexFace = ncells + j)
			//TODO 
			MatSetValue(_A, idCells[0], IndexFace,  -getOrientation(j,Ctemp1) * Fj.getMeasure() , ADD_VALUES ); 
			MatSetValue(_A, idCells[1], IndexFace,  -getOrientation(j,Ctemp2) * Fj.getMeasure() , ADD_VALUES );  	

			// LaplacianPressure
			double WaveVelocity = (_timeScheme == Implicit ) ? abs( _Velocity(j) )/2.0 + _c :  ( abs( _Velocity(j) ) + _c )/2.0 ; 
			MatSetValue(_A, idCells[0], idCells[0], - WaveVelocity * Fj.getMeasure()/2.0, ADD_VALUES ); 
			MatSetValue(_A, idCells[0], idCells[1],   WaveVelocity * Fj.getMeasure()/2.0, ADD_VALUES );  
			MatSetValue(_A, idCells[1], idCells[1], - WaveVelocity * Fj.getMeasure()/2.0, ADD_VALUES ); 
			MatSetValue(_A, idCells[1], idCells[0],   WaveVelocity * Fj.getMeasure()/2.0, ADD_VALUES );  

			if (_timeScheme == Implicit){
				double rhoL, rhoR;
				VecGetValues(_primitiveVars,1,&idCells[0],&rhoL);
				VecGetValues(_primitiveVars,1,&idCells[1],&rhoR);
				MatSetValue(_JacobianMatrix, idCells[0], idCells[0],  _c * Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_JacobianMatrix, idCells[0], idCells[1], -_c * Fj.getMeasure(), ADD_VALUES );  
				MatSetValue(_JacobianMatrix, idCells[1], idCells[1],  _c * Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_JacobianMatrix, idCells[1], idCells[0], -_c * Fj.getMeasure(), ADD_VALUES );  
				MatSetValue(_JacobianMatrix, idCells[0], IndexFace,  getOrientation(j,Ctemp1) * Fj.getMeasure() , ADD_VALUES ); 
				MatSetValue(_JacobianMatrix, idCells[1], IndexFace,  getOrientation(j,Ctemp2) * Fj.getMeasure() , ADD_VALUES );  
				MatSetValue(_JacobianMatrix, IndexFace, idCells[0],  -pow(_compressibleFluid->vitesseSon(rhoL), 2) * getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES );
				MatSetValue(_JacobianMatrix, IndexFace, idCells[1],  -pow(_compressibleFluid->vitesseSon(rhoR), 2) * getOrientation(j,Ctemp2) * Fj.getMeasure(), ADD_VALUES ); 
			}
			
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
		else if (IsSteggerBound ) { 
			/****************** Density conservation equation *********************/
			VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);	
			double u_b = getboundaryVelocity().find(j)->second * getOrientation(j,Ctemp1);
			double c_b = _compressibleFluid->vitesseSon(getboundaryPressure().find(j)->second) ;
			double lambdaPlus =  u_b + c_b;
			double lambdaMinus = u_b - c_b;
			MatSetValue(_A, idCells[0], idCells[0], -Fj.getMeasure()*lambdaPlus/2.0 * (1 - u_b/c_b)  , ADD_VALUES );
			MatSetValue(_A, idCells[0], IndexFace,  -Fj.getMeasure()*getOrientation(j,Ctemp1)*lambdaPlus/(2.0 * c_b ), ADD_VALUES );
			VecSetValue(_BoundaryTerms, idCells[0], -Fj.getMeasure()*lambdaMinus * getboundaryPressure().find(j)->second/2.0 , ADD_VALUES );
		}	
	}

	VecAssemblyBegin(_BoundaryTerms);
	VecAssemblyEnd(_BoundaryTerms);
	VecAssemblyBegin(_Conv);
	VecAssemblyEnd(_Conv);	
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	if (_timeScheme == Implicit){
		MatAssemblyBegin(_JacobianMatrix, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_JacobianMatrix, MAT_FINAL_ASSEMBLY);
	}

	Vec Au;
	VecDuplicate(_b, &Au);
	MatMult(_A,_primitiveVars, Au); //TODO AU at initialize ?

	/***********Adding boundary terms to density equation, convection term and pressure gradient to momentum equation **************/
	// Add convection terms  : AU^n + Boundterms - _Conv
	VecAXPY(Au,     1, _BoundaryTerms); 
	VecAXPY(Au,     1, _Conv); 
	//Add pressure gradient : AU^n+ Boundterms - _Conv - GradPressure 
	for (int j =0; j <_Nfaces ; j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;
		PetscInt IndexFace = _Nmailles + j;		
		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  ) 
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j) ->second).getCellsId()[0]  );
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);
			VecGetValues(_primitiveVars,1,&idCells[1],&rhoExt);
			VecSetValue(_GradPressure, IndexFace, getOrientation(j,Ctemp1) * Fj.getMeasure() *_compressibleFluid->getPressure(rhoInt) , ADD_VALUES );
			VecSetValue(_GradPressure, IndexFace, getOrientation(j,Ctemp2) * Fj.getMeasure() *_compressibleFluid->getPressure(rhoExt) , ADD_VALUES );
		}
		if (IsSteggerBound){
			VecGetValues(_primitiveVars,1,&idCells[0],&rhoInt);
			VecSetValue(_GradPressure, IndexFace, 2 * getOrientation(j,Ctemp1) * Fj.getMeasure() *_compressibleFluid->getPressure(rhoInt) , ADD_VALUES ); //TODO
		}
	}
	VecAXPY(Au,     1, _GradPressure); 
	//Multiply by inverse volums V^{-1}( AU^n + Boundterms - _Conv - GradPressure  )
	MatMult(_InvVol, Au, _b);

	if (_nbTimeStep == 0)	ComputeMinCellMaxPerim(); 
	double dt ;
	if (_timeScheme == Implicit){
		Mat Prod;
		MatMatMult(_InvVol, _JacobianMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, DIFFERENT_NONZERO_PATTERN); 
		MatDestroy(& Prod);
		//TODO create Prod at init in WaveStaggered

		dt =   _minCell / (_maxPerim * _uMax  * 2 ) ;
		MatScale(_A,  dt);
		MatShift(_A,  1);
		
	}
	else if (_timeScheme == Explicit){
		dt = _cfl * _minCell / (_maxPerim * (_uMax + _c) * 2 ) ;
	}
	VecScale(_b, dt);

	_Time.push_back(_Time.back() + dt);	
	VecDestroy(& Au); 

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
	if(_timeScheme == Explicit)
		VecCopy(_b,_newtonVariation); //DELTA U = _b = delta t*  V^{-1}Au + delta t * V^{-1}_Boundterms
	else if (_timeScheme == Implicit){
		#if PETSC_VERSION_GREATER_3_5
        KSPSetOperators(_ksp, _A, _A);
#else
        KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif
        if(_conditionNumber)  KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
		KSPSolve(_ksp, _b, _newtonVariation);
        KSPConvergedReason reason;
        KSPGetConvergedReason(_ksp,&reason);
        KSPGetIterationNumber(_ksp, &_PetscIts);
        double residu;
        KSPGetResidualNorm(_ksp,&residu);

        if (reason!=2 and reason!=3) {
            PetscPrintf(PETSC_COMM_WORLD,"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc.\n");
            PetscPrintf(PETSC_COMM_WORLD,"!!!!!!!!!!!!! Itérations maximales %d atteintes, résidu = %1.2e, précision demandée= %1.2e.\n",_maxPetscIts,residu,_precision);
            PetscPrintf(PETSC_COMM_WORLD,"Solver used %s, preconditioner %s, Final number of iteration = %d.\n",_ksptype,_pctype,_PetscIts);
            if(_mpi_rank==0){//Avoid redundant printing 
                *_runLogFile<<"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc."<<endl;
                *_runLogFile<<"!!!!!!!!!!!!! Itérations maximales "<<_maxPetscIts<<" atteintes, résidu="<<residu<<", précision demandée= "<<_precision<<endl;
                *_runLogFile<<"Solver used "<<  _ksptype<<", preconditioner "<<_pctype<<", Final number of iteration= "<<_PetscIts<<endl;
                _runLogFile->close();
            }
            if( reason == -3)		cout<<"Maximum number of iterations "<<_maxPetscIts<<" reached"<<endl;
            else if( reason == -11) cout<<"!!!!!!! Construction of preconditioner failed !!!!!!"<<endl;
            else if( reason == -5)  cout<<"!!!!!!! Generic breakdown of the linear solver (Could be due to a singular matrix or preconditioner)!!!!!!"<<endl;
            else{
                cout<<"PETSc divergence reason  "<< reason <<endl;
                cout<<"Final iteration= "<<_PetscIts<<". Maximum allowed was " << _maxPetscIts<<endl;
            }
        }
        else if( _MaxIterLinearSolver < _PetscIts)  _MaxIterLinearSolver = _PetscIts;
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
	converged=true;
	VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk+relaxation*deltaV

	/* PetscScalar q, rho_sigma, u;
	UpdateDualDensity(); // \rho^{n+1}_K -> \rho^{n+1}_\sigma 
	for (int f=0; f< _Nfaces; f++){
		bool IsWallBound =    std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),f ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),f ) != _SteggerBoundFaceSet.end() ;	
		PetscInt I = _Nmailles +f;
		if (IsWallBound) VecSetValue(_primitiveVars, I, 0, INSERT_VALUES); 
		else if (IsSteggerBound){
			VecGetValues(_DualDensity, 1,&f, &rho_sigma);
			VecSetValue(_primitiveVars, I, _Velocity(f) * rho_sigma , INSERT_VALUES); 
		}
	}	 */
	return converged;//TODO not good
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
			while ( ( i !=it2->second) && (it2 != _FacePeriodicMap.end() ) ) it2++;
			int j = (_indexFacePeriodicSet == true) && (it2 !=  _FacePeriodicMap.end())  ? it2->first : i; // in periodic k stays i, if it has been computed by scheme and takes the value of its matched face 
			int I= _Nmailles + j;
			double rho_sigma, q;
			VecGetValues(_DualDensity, 1,&j, &rho_sigma);

			if (_mpi_size > 1)
				VecGetValues(_primitiveVars_seq,1,&I,&q);
			else
				VecGetValues(_primitiveVars,1,&I,&q);

			_Velocity(i) = q / rho_sigma;

			Face Fj = _mesh.getFace(i);
			std::vector< int > idCells = Fj.getCellsId();
			Cell Ctemp1 = _mesh.getCell(idCells[0]);
			double orien1 = getOrientation(j,Ctemp1);
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
	VecDestroy(&_newtonVariation);
	VecDestroy(&_primitiveVars);
	VecDestroy(&_b);
	MatDestroy(&_InvVol); 
	MatDestroy(& _InvSurface);
	MatDestroy(&_LaplacianPressure);
	VecDestroy(& _BoundaryTerms);
	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);
	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq);
	delete[]_vec_normal;


	if (_timeScheme == Implicit) MatDestroy(& _JacobianMatrix);
	VecDestroy(& _Conv); 
	VecDestroy(& _DualDensity);
	VecDestroy(& _GradPressure);	
	delete _compressibleFluid;
	_compressibleFluid = nullptr;
	
}


