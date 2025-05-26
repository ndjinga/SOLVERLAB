/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include <numeric>
#include <math.h>
#include <iomanip>

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
	_isWall =false ;
	_Time.push_back(0.0);
	_BasisFunctionAlreadyComputed = false;
	
	if (_Ndim == 3){	
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered pas dispo en 3D, arret de calcul");				
	}
				
}

void  WaveStaggered::setboundaryVelocity(std::map< int, double> BoundaryVelocity){
	/* if (_facesBoundinit == true){ 
		std::map<int,double>::iterator it;
		for( it= BoundaryVelocity.begin(); it != BoundaryVelocity.end(); it++){
			_Velocity( it->first ) = it->second; 
		}
	}
	else{
		*_runLogFile<<"WaveStaggered::setboundaryVelocity should be called after WaveStaggered::setInitialField(Velocity)"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered::setboundaryVelocity should be called after WaveStaggered::setInitialField(Velocity)");
	} */
	_boundaryVelocity = BoundaryVelocity;
}
void  WaveStaggered::setboundaryPressure(std::map< int, double> BoundaryPressure){
	_boundaryPressure = BoundaryPressure;
}
std::map<int,double>  WaveStaggered::getboundaryPressure() const {
		return _boundaryPressure;
}
std::map<int,double>  WaveStaggered::getboundaryVelocity() const {
		return _boundaryVelocity;
}

void WaveStaggered::setWallBoundIndex(int j ){
	_WallBoundFaceSet.push_back(j);
	_isWall = true;
}
void WaveStaggered::setSteggerBoundIndex(int j ){
	_SteggerBoundFaceSet.push_back(j);
}
void WaveStaggered::setInteriorIndex(int j ){
	_InteriorFaceSet.push_back(j);
}

bool  WaveStaggered::IsFaceBoundaryNotComputedInPeriodic(int j ) {
	std::map<int,int>::iterator it = _FacePeriodicMap.begin();
	while ( ( j !=it->second) && (it !=_FacePeriodicMap.end() ) )
		it++;
	return _indexFacePeriodicSet && (it !=  _FacePeriodicMap.end());
}
bool  WaveStaggered::IsFaceBoundaryComputedInPeriodic(int j ) {
	return _FacePeriodicMap.find(j) != _FacePeriodicMap.end() &&  (std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end()) ;
}


void WaveStaggered::setOrientation(int j,std::vector<double> vec_normal_sigma){
	for (int idim = 0; idim < _Ndim; ++idim)
		_vec_sigma[j].push_back(vec_normal_sigma[idim]);
}
//TODO what about checkorientation in periodic ?
double WaveStaggered::getOrientation(int l, Cell Cint) {	
	double checktouternormal;
	std::vector<double> vec(_Ndim, 0.0);
	std::map<int,int>::iterator it = _FacePeriodicMap.begin();
	while ( ( l !=it->second) && (it !=_FacePeriodicMap.end() ) )it++;
	for(int m=0; m<Cint.getNumberOfFaces(); m++){//we look for l the index of the face Fj for the cell Ctemp1
		if (l == Cint.getFacesId()[m]  ){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);
			checktouternormal = ( _mesh.getFace(l).x()  - Cint.x() ) *vec[0];
			if(_Ndim ==2) checktouternormal+= ( _mesh.getFace(l).y() - Cint.y()  ) *vec[1];
			if (checktouternormal <-1e-11){
				for (int idim = 0; idim < _Ndim; ++idim)
					vec[idim] = -Cint.getNormalVector(m,idim);

			}
		}
		else if ((_FacePeriodicMap.find(l) != _FacePeriodicMap.end()) && (_FacePeriodicMap.find(l)->second == Cint.getFacesId()[m])){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);
			checktouternormal = ( _mesh.getFace(_FacePeriodicMap.find(l)->second ).x()  - Cint.x() ) *vec[0];
			if(_Ndim ==2) checktouternormal+= ( _mesh.getFace(_FacePeriodicMap.find(l)->second ).y() - Cint.y()  ) *vec[1];
			if (checktouternormal <-1e-11){
				for (int idim = 0; idim < _Ndim; ++idim)
					vec[idim] = -Cint.getNormalVector(m,idim);
			}
		}
		else if (it !=_FacePeriodicMap.end() && ( it->first == Cint.getFacesId()[m])){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);
			checktouternormal = ( _mesh.getFace(it->first ).x()  - Cint.x() ) *vec[0];
			if(_Ndim ==2) checktouternormal+= ( _mesh.getFace(it->first ).y() - Cint.y()  ) *vec[1];
			if (checktouternormal <-1e-11){
				for (int idim = 0; idim < _Ndim; ++idim)
					vec[idim] = -Cint.getNormalVector(m,idim);

			}

		}
	}
	double dotprod = 0;
	for (int idim = 0; idim < _Ndim; ++idim)
		dotprod += vec[idim] * _vec_sigma.find(l)->second[idim]; 
	return dotprod ; 

}


void WaveStaggered::setPeriodicFaces( Mesh &M, const char &Direction, int ncells, double inf, double sup){ 
	for (int j=0;j<M.getNumberOfFaces() ; j++){
		Face my_face=M.getFace(j);
		double e;
		double tol = 1.0/(ncells *4);
		if (Direction == 'x')
			e=my_face.x();
		else if (Direction == 'y' && _Ndim ==2)
			e=my_face.y();
		if (my_face.getNumberOfCells() ==1 ){ 
			if( (_Ndim==2) &&  e>tol+ inf && e< (sup-tol) ){
				for (int iface=0;iface<M.getNumberOfFaces() ; iface++){
					Face face_i=M.getFace(iface);
					double ei;
					if (Direction == 'x')
						ei=face_i.x();
					else if (Direction == 'y')
						ei=face_i.y();
					if (face_i.getNumberOfCells() ==1 && iface !=j && ( abs(e-ei)<tol) && (_FacePeriodicMap.find(iface) == _FacePeriodicMap.end())){ 
						_FacePeriodicMap[j]=iface;
						setInteriorIndex(j);
					}
				}
			}
			else if (_Ndim == 1){
				for (int iface=0;iface<M.getNumberOfFaces() ; iface++){
					Face face_i=M.getFace(iface);
					if (face_i.getNumberOfCells() ==1 && iface !=j && (_FacePeriodicMap.find(iface) == _FacePeriodicMap.end())){ 
						_FacePeriodicMap[j]=iface;
						setInteriorIndex(j);
					}
				}	
			}	
		}
	}
	_indexFacePeriodicSet = true;
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
		_Velocity_at_Cells = Field("Velocity at cells results", CELLS, _mesh,3);
		_DivVelocity = Field("velocity divergence", CELLS, _mesh, 1);	 
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

	double * initialFieldVelocity = new double[_Nfaces];
	double * initialFieldPressure = new double[_Nmailles];
	
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

	_d = 1/( sqrt(2*_neibMaxNbCells) );
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

	// Création matrice Q tq U^n+1 - U^n = dt V^{-1} _A U^n pour schéma explicite
	MatCreate(PETSC_COMM_SELF, & _A); 
	MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _globalNbUnknowns, _globalNbUnknowns );
	MatSetFromOptions(_A);
	MatSetUp(_A);
	MatZeroEntries(_A);

	/********* Metrics Related ********/
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

	/********* Pressure Related ********/
	// matrice DIVERGENCE (|K|div(u))
	MatCreate(PETSC_COMM_SELF, & _Div); 
	MatSetSizes(_Div, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces );
	MatSetFromOptions(_Div);
	MatSetUp(_Div);
	MatZeroEntries(_Div);

	// matrix _LaplacianPressure (without boundary terms)
	MatCreate(PETSC_COMM_SELF, & _LaplacianPressure); 
	MatSetSizes(_LaplacianPressure, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles ); 
	MatSetFromOptions(_LaplacianPressure);
	MatSetUp(_LaplacianPressure);
	MatZeroEntries(_LaplacianPressure);

	// Vector BoundaryTerms for Pressure
	VecCreate(PETSC_COMM_SELF, & _BoundaryTerms); 
	VecSetSizes(_BoundaryTerms, PETSC_DECIDE, _globalNbUnknowns); 
	VecSetFromOptions(_BoundaryTerms);
	VecSetUp(_BoundaryTerms);
	VecZeroEntries(_BoundaryTerms);

	/********* Velocity Related ********/
	// matrix GRADIENT (we will impose to _Be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _DivTranspose); 
	MatSetSizes(_DivTranspose, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nmailles );
	MatSetFromOptions(_DivTranspose);
	MatSetUp(_DivTranspose);
	MatZeroEntries(_DivTranspose);

	// matrix GRAD DIV TILDE (we will impose to _Be 0 on faces so that u^n+1 = u^n at the _Boundary)
	MatCreate(PETSC_COMM_SELF, & _GradDivTilde); 
	MatSetSizes(_GradDivTilde, PETSC_DECIDE, PETSC_DECIDE, _Nfaces, _Nfaces );
	MatSetFromOptions(_GradDivTilde);
	MatSetUp(_GradDivTilde);
	MatZeroEntries(_GradDivTilde);

	AssembleMetricsMatrices();
	if(_system)
	{
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

	computeHodgeDecompositionWithBoundaries();

	save();//save initial data
}

//TOD0 integration formula is not good 
double WaveStaggered::MassLumping(const Cell &K, const int &idcell, const Face & Facej, const int &j){
	double masslumping_on_K =0;
	std::vector<Cell> Support_f, Support_j ;
	Support_f.push_back(K);
	Support_j.push_back(K);
	std::vector<Node> K_Nodes;
	for (int i=0; i < K.getNodesId().size(); i++)
		K_Nodes.push_back(_mesh.getNode(K.getNodesId()[i]) );

	for (int l =0; l <K.getNumberOfFaces(); l ++){
		Point Xl =  _mesh.getFace( K.getFacesId()[l] ).getBarryCenter();
		double weight = (_Ndim ==2) ? ( (K.getNumberOfFaces() == 4) ? 1/4.0 : 1/6.0 ) : 1/2.0;
		double Area = fabs( det( JacobianTransfor_K_X( xToxhat(K, Xl, K_Nodes), K_Nodes) ) ) * weight;
		std::vector<double> Psi_j_in_Xl = PhysicalBasisFunctionRaviartThomas(K, idcell, Support_j, Facej,j, Xl);

		for (int f=0; f <K.getNumberOfFaces(); f ++){
			Face Facef = _mesh.getFace( K.getFacesId()[f] );
			std::vector<double> Psi_f_in_Xl = PhysicalBasisFunctionRaviartThomas(K, idcell, Support_f, Facef, K.getFacesId()[f], Xl);
			for (int e=0; e<_Ndim; e++)
				masslumping_on_K +=  Area *  Psi_j_in_Xl[e] * Psi_f_in_Xl[e];	
		} 
	}		
	return masslumping_on_K;
}


void WaveStaggered::AssembleMetricsMatrices(){
	for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
 		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		std::vector< int > idFaces = Ctemp1.getFacesId();
		PetscScalar det, InvD_sigma, InvVol1, InvVol2, InvPerimeter1, InvPerimeter2;
		PetscInt IndexFace = _Nmailles + j;
		double D_sigma_K, D_sigma_L,HalfDiamondCell;
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;	

		D_sigma_K =  MassLumping(Ctemp1, idCells[0], Fj, j);
		InvPerimeter1 = (_Ndim ==2) ? 1.0/_perimeters(idCells[0]) : 1.0 ;
		InvVol1 = 1.0/Ctemp1.getMeasure();

		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end())
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) ? _mesh.getFace(_FacePeriodicMap.find(j)->second ) :  Fj;

			D_sigma_L = MassLumping(Ctemp2, idCells[1], Fj_physical, j);
			InvPerimeter2 = (_Ndim ==2) ? 1.0/_perimeters(idCells[1]) : 1.0 ;
			InvVol2 = 1.0/Ctemp2.getMeasure();

			InvD_sigma = 1.0/(D_sigma_K + D_sigma_L);
			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1 , INSERT_VALUES );
			MatSetValue(_InvVol, idCells[1],idCells[1], InvVol2, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES); 	
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvSurface,idCells[1],idCells[1], InvPerimeter2, INSERT_VALUES );
		}
		else if (IsWallBound || IsSteggerBound ) { 
			InvD_sigma = 1.0/D_sigma_K ; 	
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES);  
		}
	}
	MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(_InvSurface,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
}



void  WaveStaggered::computeHodgeDecompositionWithBoundaries(){
	assert(_nbTimeStep == 0);
    Mat B, Btranspose, HodgeLaplacian;
    Vec phi, b, nullVec, Bu, u_0, u_phi;
    KSP ksp;
    MatNullSpace nullspace;

	//Transfer u_0 into PETSC VEC
	VecCreateSeq(PETSC_COMM_SELF, _Nfaces, &u_0);
	VecZeroEntries(u_0);
	for(int i =0; i<_Nfaces; i++){
		VecSetValue(u_0, i, _Velocity(i), INSERT_VALUES);
	}
	
    //Hodge Laplacian is created from BB^t with B =-div 
    MatCreate(PETSC_COMM_WORLD, &HodgeLaplacian);
    MatSetSizes(HodgeLaplacian, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nmailles);
    MatSetFromOptions(HodgeLaplacian);
    MatSetUp(HodgeLaplacian);

	MatCreate(PETSC_COMM_WORLD, &B);
    MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, _Nmailles, _Nfaces);
    MatSetFromOptions(B);
    MatSetUp(B);
	MatZeroEntries(B); 

	for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
		if ( IsInterior ){	 //TODO periodic ?
			Cell Ctemp1 = _mesh.getCell(idCells[0]);
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			MatSetValue(B, idCells[0], j ,  -getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES ); 
			MatSetValue(B, idCells[1], j , - getOrientation(j,Ctemp2) * Fj.getMeasure(), ADD_VALUES );  
		}
	}
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
	MatTranspose(B, MAT_INITIAL_MATRIX, &Btranspose);
	MatMatMult(B,Btranspose, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &HodgeLaplacian);

    // Eliminating constants from the kernel
    VecCreateSeq(PETSC_COMM_SELF, _Nmailles, &nullVec);
    VecSet(nullVec, 1.0);  
	VecNormalize(nullVec, NULL);
   	MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &nullVec, &nullspace);
    MatSetNullSpace(HodgeLaplacian, nullspace);

    // Source term
		// -/int_{\partial \Omega} u_b\cdot n d\Gamma
    VecCreateSeq(PETSC_COMM_SELF, _Nmailles, &b);
	VecZeroEntries(b);
	VecDuplicate(b, &Bu);
    for (PetscInt i = 0; i < _Nmailles; ++i) {
		Cell K = _mesh.getCell(i);
		std::vector<int> FacesOfKi = K.getFacesId();
		for (int nei =0; nei <FacesOfKi.size(); nei ++){
			if (_mesh.getFace(FacesOfKi[nei]).getNumberOfCells() ==1 )
				VecSetValue(b, i, -_mesh.getFace(FacesOfKi[nei]).getMeasure() * getboundaryVelocity().find(FacesOfKi[nei])->second * getOrientation(FacesOfKi[nei], K) , ADD_VALUES);
		}  
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
		// /int_\Omega u_0 \cdot (-div)^*f dx
	MatMult(B, u_0, Bu);
	VecAXPY(b, 1, Bu);

    // Project source term on range of HodgeLaplacian
    MatNullSpaceRemove(nullspace, b);

    VecDuplicate(b, &phi);

    // --- Solver KSP ---
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, HodgeLaplacian, HodgeLaplacian);
    KSPSetType(ksp, KSPCG);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	PC pc;
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCJACOBI);
    KSPSetUp(ksp);

    KSPSolve(ksp, b, phi);

	// Construct (u_0)_psi := u_0  -  (Bphi + tr(u - u_b) )
		// First  u_phi := (Bphi + tr(u - u_b) )
	VecDuplicate(u_0, &u_phi);
	MatMultTranspose(B, phi, u_phi);
	for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
		if ( !IsInterior )	 //TODO periodic ?
			VecSetValue(u_phi, j , _Velocity(j) - getboundaryVelocity().find(j)->second ,INSERT_VALUES);
	}
		//Then (u_0)_psi := u_0  -  u_phi
	_Velocity_0_Psi = Field("(Velocity_0)_Psi", FACES, _mesh,3);
	_Velocity_0_Psi_at_Cells = Field("Velocity_Psi", CELLS, _mesh,3);
	double u_0_j, u_phi_j;
	for (int j=0; j<_Nfaces;j++){ 
		VecGetValues(u_0, 1, &j, &u_0_j);
		VecGetValues(u_phi, 1, &j, &u_phi_j);
		_Velocity_0_Psi(j) = u_0_j - u_phi_j;
	}
	// Saving the Field
	_Velocity_0_Psi.setTime(_time,_nbTimeStep);
	_Velocity_0_Psi.setInfoOnComponent(0,"(Velocity_0)_Psi. n_sigma_(m/s)");

	//check that it is divergence free
	for (int m=0; m <_Nmailles; m++){
		double div = 0;
		Cell K = _mesh.getCell(m);
		std::vector<int> idFaces = K.getFacesId();
		for (int f=0; f<K.getNumberOfFaces(); f++){
			Face Facef = _mesh.getFace(idFaces[f] );
			div += Facef.getMeasure() * getOrientation(idFaces[f], K) * _Velocity_0_Psi(idFaces[f]);	
		}
		if (fabs(div)>1e-10)
			cout << " div = "<< div <<endl;
		assert( fabs(div)<1e-10 );
	}
	//check that it is equal to u_b on d\Omega
	for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
		double delta;
		if ( !IsInterior ){	 //TODO periodic ?
			delta = fabs(getboundaryVelocity().find(j)->second - _Velocity_0_Psi(j));
		}
		if (fabs(delta)>1e-10)
			cout << " tr(ub - u_psi) = "<< delta <<endl;
		assert(fabs(delta)<1e-10);
		if (IsInterior){
			VecGetValues(u_phi, 1, &j, &u_phi_j);
			double u = _Velocity_0_Psi(j) + u_phi_j;
			assert( fabs(_Velocity(j) - u) <1e-10 );
		}
	}
	
    // Delete
    VecDestroy(&phi);
    VecDestroy(&b);
	VecDestroy(&Bu);
	VecDestroy(&u_0);
	VecDestroy(&u_phi);
    VecDestroy(&nullVec);
    MatDestroy(&HodgeLaplacian);
	MatDestroy(&B);
	MatDestroy(&Btranspose);
	//MatDestroy(&C);
    MatNullSpaceDestroy(&nullspace);
    KSPDestroy(&ksp);
}

double WaveStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme
	//The matrices are assembled only in the first time step since linear problem
	MatZeroEntries(_LaplacianPressure);
	MatZeroEntries(_Div); 
	MatZeroEntries(_DivTranspose); 
	VecZeroEntries(_b); 

	if ( _nbTimeStep == 0  ){   
		// Assembly of matrices 
		for (int j=0; j<_Nfaces;j++){
			Face Fj = _mesh.getFace(j);
			std::vector< int > idCells = Fj.getCellsId();
			Cell Ctemp1 = _mesh.getCell(idCells[0]);
			PetscInt IndexFace = _Nmailles + j;

			bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
			bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
			bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;

			if ( IsInterior ){	
				if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  )
					idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
				
				/******************* Pressure equation ***********************/
				MatSetValue(_A, idCells[0], IndexFace, -1/_rho * getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_A, idCells[1], IndexFace, -1/_rho * getOrientation(j,Ctemp2) * Fj.getMeasure(), ADD_VALUES );  
				MatSetValue(_A, idCells[0], idCells[0], -_d*_c *Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_A, idCells[0], idCells[1], _d*_c *Fj.getMeasure(), ADD_VALUES );  
				MatSetValue(_A, idCells[1], idCells[1], -_d*_c *Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_A, idCells[1], idCells[0], _d*_c *Fj.getMeasure(), ADD_VALUES );  
				
				/******************* Velocity equation ***********************/
				MatSetValue(_A, IndexFace, idCells[0], _kappa * getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES ); 
				MatSetValue(_A, IndexFace, idCells[1], _kappa *getOrientation(j,Ctemp2) * Fj.getMeasure(), ADD_VALUES );
				for (int nei =0; nei <2; nei ++){ // we are in the case where an interior face has two neighbours
					Cell K = _mesh.getCell(idCells[nei]); 
					std::vector<int> idFaces = K.getFacesId();
					//If  the face is periodic and K isn't the direct neighbour of Fj, recover the geometric informations of the associated periodic cell 
					Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[nei] ? _mesh.getFace(_FacePeriodicMap.find(j)->second ):  Fj;
					
					for (int f =0; f <K.getNumberOfFaces(); f ++){
						Face Facef = _mesh.getFace( idFaces[f] );
						std::vector< int> idCellsOfFacef =  Facef.getCellsId();
						std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
						PetscInt I = _Nmailles + idFaces[f];
						double gradiv = - _d * _c * Fj_physical.getMeasure() * getOrientation(j,K) *Facef.getMeasure()* getOrientation(idFaces[f], K) /( (_Ndim==2 )? _perimeters[idCells[nei]] : 1.0);
						MatSetValue(_A, IndexFace, I, gradiv, ADD_VALUES ); 
					}
				}	
			}
			else if (IsSteggerBound || IsWallBound ) { 
				//MatSetValue(_A, idCells[0], IndexFace, -1/_rho * getOrientation(j,Ctemp1) * Fj.getMeasure(), ADD_VALUES ); 
				
				if (IsSteggerBound){
					//********* pressure equation *************//
					MatSetValue(_A, idCells[0], IndexFace,  -1/_rho * getOrientation(j,Ctemp1) * Fj.getMeasure()/2.0, ADD_VALUES );
					VecSetValue(_BoundaryTerms, idCells[0], -1/_rho * getOrientation(j,Ctemp1) * Fj.getMeasure()/2.0 * getboundaryVelocity().find(j)->second, ADD_VALUES );
					MatSetValue(_A, idCells[0], idCells[0], -_d * _c * Fj.getMeasure(), ADD_VALUES );
					VecSetValue(_BoundaryTerms, idCells[0],  _d * _c * Fj.getMeasure() * getboundaryPressure().find(j)->second, ADD_VALUES );
					//************* Velocity equation *************//
					// pressure gradient
					MatSetValue(_A, IndexFace , idCells[0],  _kappa * getOrientation(j,Ctemp1) * Fj.getMeasure()/2.0, ADD_VALUES );
					VecSetValue(_BoundaryTerms, IndexFace , -_kappa * getOrientation(j,Ctemp1) * Fj.getMeasure()/2.0 * getboundaryPressure().find(j)->second, ADD_VALUES );
					// jump velocity
					MatSetValue(_A, IndexFace , IndexFace,  -_c/2.0  * Fj.getMeasure(), ADD_VALUES );
					VecSetValue(_BoundaryTerms, IndexFace , _c/2.0 * Fj.getMeasure()* getboundaryVelocity().find(j)->second, ADD_VALUES ); 
					
					Cell K = _mesh.getCell(idCells[0]); 
					std::vector<int> idFaces = K.getFacesId();
					//If  the face is periodic and K isn't the direct neighbour of Fj, recover the geometric informations of the associated periodic cell 
					Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) && Fj.getCellsId()[0] != idCells[0] ? _mesh.getFace(_FacePeriodicMap.find(j)->second ):  Fj;
					
					for (int f =0; f <K.getNumberOfFaces(); f ++){
						Face Facef = _mesh.getFace( idFaces[f] );
						std::vector< int> idCellsOfFacef =  Facef.getCellsId();
						std::vector< int > idNodesOfFacef = Facef.getNodesId(); 
						PetscInt I = _Nmailles + idFaces[f];
						double gradiv = - _d * _c * Fj_physical.getMeasure() * getOrientation(j,K) *Facef.getMeasure()* getOrientation(idFaces[f], K) /( (_Ndim==2 )? _perimeters[idCells[0]] : 1.0);
						MatSetValue(_A, IndexFace, I, gradiv, ADD_VALUES ); 
					}
				}
				if (IsWallBound){
					// Velocity equation
					// jump velocity
					MatSetValue(_A, IndexFace , IndexFace, -_c * Fj.getMeasure(), ADD_VALUES );
				}
			}	
		}
		MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
		VecAssemblyBegin(_BoundaryTerms);
		VecAssemblyEnd(_BoundaryTerms);
		Mat Prod;
		MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
		MatCopy(Prod,_A, SAME_NONZERO_PATTERN); 
		MatDestroy(& Prod);


		ComputeMinCellMaxPerim();

	}

	double dt;
	Vec Prod2;
	VecDuplicate(_BoundaryTerms, &Prod2);
	MatMult(_InvVol, _BoundaryTerms, Prod2);  
	if (_timeScheme == Explicit){
		MatMult(_A,_primitiveVars, _b); 
		dt =  _cfl * _minCell / (_maxPerim * _c ) ;
		VecScale(_b,  dt);
	}
	else if (_timeScheme == Implicit){
		VecCopy(_primitiveVars, _b); 
		dt =  10 * _minCell / (_maxPerim * _c ) ;
		if (_nbTimeStep ==0){
			MatScale(_A,  - dt);
			MatShift(_A,  1);
		}
		
	}
	VecAXPY(_b,     dt , Prod2);
	VecDestroy(& Prod2); 
		
	double PreviousTime = _Time.back();
	_Time.push_back(PreviousTime+ dt);
	ComputeEnergyAtTimeT();
	
	return dt ;
}



void WaveStaggered::ComputeMinCellMaxPerim(){
	if(_nbTimeStep == 0){
		Vec V, W;
		PetscScalar minInvSurf, maxInvVol;
		
		VecCreate(PETSC_COMM_SELF, & V);
		VecCreate(PETSC_COMM_SELF, & W);
		VecSetSizes(V, PETSC_DECIDE, _globalNbUnknowns);
		VecSetSizes(W, PETSC_DECIDE, _Nmailles);
		VecSetFromOptions(V);
		VecSetFromOptions(W);

		// Minimum size of mesh volumes
		MatGetDiagonal(_InvVol,V);
		VecMax(V, NULL, &maxInvVol);
		_minCell = 1.0/maxInvVol;
		//Maximum size of surfaces
		MatGetDiagonal(_InvSurface, W);
		VecMin(W, NULL, &minInvSurf);
		_maxPerim = 1.0/minInvSurf;

		VecDestroy(& V);
		VecDestroy(& W);
	} 
}

void WaveStaggered::computeNewtonVariation(){
	if(_timeScheme == Explicit)
		VecCopy(_b,_newtonVariation); //DELTA U = _b = delta t*  V^{-1}Au + delta t * V^{-1}_Boundterms
	else if (_timeScheme == Implicit)
		VecCopy(_primitiveVars, _newtonVariation); //DELTA U = U^n
}

bool WaveStaggered::iterateTimeStep(bool &converged){
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
	if(_timeScheme == Explicit){
		converged=true;
		VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk+relaxation*deltaV
	}
	else if (_timeScheme == Implicit){
		#if PETSC_VERSION_GREATER_3_5
        KSPSetOperators(_ksp, _A, _A);
#else
        KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif
        if(_conditionNumber)  KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
		KSPSolve(_ksp, _b, _primitiveVars);
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
            converged = false;
        }
        else {
            if( _MaxIterLinearSolver < _PetscIts)  _MaxIterLinearSolver = _PetscIts;
            VecAXPY(_newtonVariation,  -1, _primitiveVars );	//DELTA U = U^n - U^{n+1}
			VecScale(_newtonVariation, -1.0 );					//DELTA U = - (U^n - U^{n+1})       
            converged =  true ; //TODO (_erreur_rel <= _precision) ;
        }
	}	
	
	return converged;

}
void WaveStaggered::validateTimeStep(){
	//Calcul de la variation Un+1-Un
	_erreur_rel= 0;
	double x, dx;
	for(int j=0; j<_globalNbUnknowns; j++){
		VecGetValues(_newtonVariation, 1, &j, &dx);
		VecGetValues(_primitiveVars, 1, &j, &x);
		if (fabs(x)< _precision){
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

bool WaveStaggered::initTimeStep(double dt){
	_dt = dt;
	return _dt>0;//No need to call MatShift as the linear system matrix is filled at each Newton iteration (unlike linear problem)
}

void WaveStaggered::abortTimeStep(){
	_dt = 0;
}

void WaveStaggered::ComputeEnergyAtTimeT(){
	double E = 0;
	double pb = getboundaryPressure().begin()->second;	//Warning pb should be constant

	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		PetscInt I = _Nmailles + j;
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		PetscScalar InvD_sigma, InvCell1measure, InvCell2measure, pressure_in, pressure_out, velocity;

		/* double r =  sqrt(Fj.x()*Fj.x() + Fj.y()*Fj.y());
		double theta = atan(Fj.y()/Fj.x());
		double dotprodExact = pow(r1,2)/(pow(r1,2) - pow(r0,2))*( (1 - pow(r0,2)/pow(r,2) * cos(2*theta)) *_vec_sigma.find(j)->second[0] - pow(r0,2)/pow(r,2) * sin(2*theta) *_vec_sigma.find(j)->second[1]  ); 
		*/
		if (Fj.getNumberOfCells()==2  ){	// Fj is inside the domain or is a boundary periodic face (computed)
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			MatGetValues(_InvVol, 1, &I,1, &I, &InvD_sigma);
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell1measure );
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell2measure );
			VecGetValues(_primitiveVars, 1, &idCells[0], &pressure_in );
			VecGetValues(_primitiveVars, 1, &idCells[1], &pressure_out );
			VecGetValues(_primitiveVars, 1, &I, &velocity );

			double pressure_int=  1/(InvCell1measure*Ctemp1.getNumberOfFaces()) *pow( (pressure_in - pb) ,2) ;
			double pressure_ext=  1/(InvCell2measure*Ctemp2.getNumberOfFaces()) *pow( (pressure_out - pb) ,2) ;
			double velocity_part = 1/(InvD_sigma) *pow( (velocity - _Velocity_0_Psi(j)) ,2) ; 
			E += pressure_int + pressure_ext + velocity_part ;
						
		}
		else if (Fj.getNumberOfCells()==1 ) { //if boundary face and face index is different from periodic faces not computed 	
			MatGetValues(_InvVol,1,&I, 1, &I,&InvD_sigma);
			MatGetValues(_InvVol, 1, &idCells[0],1, &idCells[0], &InvCell1measure );
			VecGetValues(_primitiveVars, 1, &idCells[0], &pressure_in );
			VecGetValues(_primitiveVars, 1, &I, &velocity );

			double pressure_part_cellint=  1/(InvCell1measure*Ctemp1.getNumberOfFaces()) *pow( (pressure_in - pb) ,2) ;
			double velocity_part = 1/(InvD_sigma) *pow( (velocity - _Velocity_0_Psi(j)) ,2) ; 
			E += pressure_part_cellint + velocity_part ;
		}	
	}
	_Energy.push_back(E);
}


//################ Raviart-Thomas ##############//
std::vector<double> WaveStaggered::ReferenceBasisFunctionRaviartThomas(const int &i, const Point &Xhat, const std::vector<Node> &K_Nodes){
	std::vector<double> Psihat(_Ndim, 0.0);
	if (_Ndim ==1){
		if (i ==0)   Psihat[0] = Xhat.x()-1;
		else if (i ==1)   Psihat[0] = Xhat.x();
	}
	else if (_Ndim ==2){
		if ( K_Nodes.size() == 4){
			if (i ==0)      Psihat[0] = Xhat.x();
			else if (i ==1) Psihat[1] = Xhat.y();
			else if (i ==2)	Psihat[0]=  Xhat.x()-1;
			else if (i ==3) Psihat[1] = Xhat.y()-1;	
			double areaKleft = fabs( (K_Nodes[0].x()-K_Nodes[1].x())* (K_Nodes[2].y()-K_Nodes[1].y()) - (K_Nodes[0].y()-K_Nodes[1].y())* (K_Nodes[2].x()-K_Nodes[1].x())  )/2.0;
			double areaKright = fabs( (K_Nodes[2].x()-K_Nodes[3].x())* (K_Nodes[0].y()-K_Nodes[3].y()) - (K_Nodes[2].y()-K_Nodes[3].y())* (K_Nodes[0].x()-K_Nodes[3].x())  )/2.0;
			double distortedQuadsX = ((K_Nodes[0].x()-K_Nodes[1].x()) * (K_Nodes[3].y()-K_Nodes[2].y()) - (K_Nodes[0].y()-K_Nodes[1].y())* ( K_Nodes[3].x()-K_Nodes[2].x() ) )/(2.0 *(areaKleft+ areaKright) );
			double distortedQuadsY = ((K_Nodes[1].x()-K_Nodes[2].x()) * (K_Nodes[3].y()-K_Nodes[0].y()) - (K_Nodes[1].y()-K_Nodes[2].y())* ( K_Nodes[3].x()-K_Nodes[0].x() ) )/(2.0 *(areaKleft+ areaKright) );
			Psihat[0] += distortedQuadsX * Xhat.x() *(Xhat.x() - 1 );
			Psihat[1] += distortedQuadsY * Xhat.y() *(Xhat.y() - 1 );
		}
		if ( K_Nodes.size() == 3){
			if (i ==0){      
				Psihat[0] = Xhat.x();
				Psihat[1] = Xhat.y() - 1;
			}
			else if (i ==1){      
				Psihat[0] = Xhat.x();
				Psihat[1] = Xhat.y();
			}
			else if (i ==2){      
				Psihat[0] = Xhat.x()-1 ;
				Psihat[1] = Xhat.y();
			}
		}
	}
	return Psihat;
}

Point WaveStaggered::xToxhat(const Cell & K, const Point &X, const std::vector<Node> & K_Nodes){
	Point Xhat;
	Point Xb = K.getBarryCenter();
	if (_Ndim ==1){
		if (fabs(X.x() - Xb.x()) < 1e-11) // if X is the center of gravity
			Xhat = Point(0.5, 0, 0);
		else if (fabs(X.x() - K_Nodes[0].x()) < 1e-11) // if X is the Node X_1
			Xhat = Point(0, 0, 0);
		else if (fabs(X.x() - K_Nodes[1].x()) < 1e-11) // if X is the Node X_2
			Xhat = Point(1, 0, 0);
	}
	else if (_Ndim ==2){
		if (K_Nodes.size() == 4){
			if (fabs(X.x() - Xb.x()) < 1e-11 && fabs(X.y() - Xb.y()) < 1e-11) //if X is the center of gravity
				Xhat = Point(0.5, 0.5, 0);
			//Nodes related 
			else if (fabs(X.x() - K_Nodes[1].x()) < 1e-11 && fabs(X.y() - K_Nodes[1].y()) < 1e-11) // if X is the Node X_1
				Xhat = Point(0, 0, 0);
			else if (fabs(X.x() - K_Nodes[2].x()) < 1e-11 && fabs(X.y() - K_Nodes[2].y()) < 1e-11) // if X is the Node X_2
				Xhat = Point(0, 1, 0);
			else if (fabs(X.x() - K_Nodes[3].x()) < 1e-11 && fabs(X.y() - K_Nodes[3].y()) < 1e-11) // if X is the Node X_3
				Xhat = Point(1, 1, 0);
			else if (fabs(X.x() - K_Nodes[0].x()) < 1e-11 && fabs(X.y() - K_Nodes[0].y()) < 1e-11) // if X is the Node X_4
				Xhat = Point(1, 0, 0);
			//Faces related 
			else if (fabs(X.x() - (K_Nodes[1].x() + K_Nodes[2].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[1].y() + K_Nodes[2].y())/2.0) < 1e-11) // if X is the middle point of the first face
				Xhat = Point(0, 0.5, 0);
			else if (fabs(X.x() - (K_Nodes[2].x() + K_Nodes[3].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[2].y() + K_Nodes[3].y())/2.0) < 1e-11) // if X is the middle point of the second face
				Xhat = Point(0.5, 1, 0);
			else if (fabs(X.x() - (K_Nodes[3].x() + K_Nodes[0].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[3].y() + K_Nodes[0].y())/2.0) < 1e-11) // if X is the middle point of the third face
				Xhat = Point(1, 0.5, 0);
			else if (fabs(X.x() - (K_Nodes[0].x() + K_Nodes[1].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[0].y() + K_Nodes[1].y())/2.0) < 1e-11) // if X is the middle point of the fourth face
				Xhat = Point(0.5, 0, 0);
		}
		else if (K_Nodes.size() == 3){
			if (fabs(X.x() - Xb.x()) < 1e-11 && fabs(X.y() - Xb.y()) < 1e-11) // if X is the center of gravity
				Xhat = Point(1.0/3.0, 1.0/3.0, 0); 
			// Nodes related 
			else if (fabs(X.x() - K_Nodes[0].x()) < 1e-11 && fabs(X.y() - K_Nodes[0].y()) < 1e-11) // if X is the Node X_1
				Xhat = Point(0, 0, 0);
			else if (fabs(X.x() - K_Nodes[1].x()) < 1e-11 && fabs(X.y() - K_Nodes[1].y()) < 1e-11) // if X is the Node X_2
				Xhat = Point(0, 1, 0);
			else if (fabs(X.x() - K_Nodes[2].x()) < 1e-11 && fabs(X.y() - K_Nodes[2].y()) < 1e-11) // if X is the Node X_3
				Xhat = Point(1, 0, 0);
			// Faces related 
			else if (fabs(X.x() - (K_Nodes[0].x() + K_Nodes[1].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[0].y() + K_Nodes[1].y())/2.0) < 1e-11) // if X is the middle point of the first face
				Xhat = Point(0, 0.5, 0); 
			else if (fabs(X.x() - (K_Nodes[1].x() + K_Nodes[2].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[1].y() + K_Nodes[2].y())/2.0) < 1e-11) // if X is the middle point of the second face
				Xhat = Point(0.5, 0.5, 0);
			else if (fabs(X.x() - (K_Nodes[2].x() + K_Nodes[0].x())/2.0) < 1e-11 && fabs(X.y() - (K_Nodes[2].y() + K_Nodes[0].y())/2.0) < 1e-11) // if X is the middle point of the third face
				Xhat = Point(0.5, 0, 0);
		}
		
	}	
	return Xhat;
}


std::vector<double> WaveStaggered::JacobianTransfor_K_X(const Point &X, const std::vector<Node> &K_Nodes){
	std::vector<double> JacobianTransfor_K(_Ndim * _Ndim);
	if (_Ndim ==1)  JacobianTransfor_K[0] = (K_Nodes[1].x()-K_Nodes[0].x()) ; 
	else if (_Ndim ==2){ 
		if (K_Nodes.size()== 4){
			//first column : should only depend on Xhat_y
			JacobianTransfor_K[0] = ( (K_Nodes[0].x()  - K_Nodes[1].x())*(1 - X.y()) + (K_Nodes[3].x()  - K_Nodes[2].x())*X.y() );
			JacobianTransfor_K[2] = ( (K_Nodes[0].y()  - K_Nodes[1].y())*(1 - X.y()) + (K_Nodes[3].y()  - K_Nodes[2].y())*X.y() );
			//second column: should only depend on Xhat_x
			JacobianTransfor_K[1] = ( (K_Nodes[2].x()  - K_Nodes[1].x())*(1 - X.x()) + (K_Nodes[3].x()  - K_Nodes[0].x())*X.x()  );
			JacobianTransfor_K[3] = ( (K_Nodes[2].y()  - K_Nodes[1].y())*(1 - X.x()) + (K_Nodes[3].y()  - K_Nodes[0].y())*X.x()  );			
		}
		else if (K_Nodes.size()== 3){
			JacobianTransfor_K[0] = K_Nodes[2].x() - K_Nodes[0].x();
			JacobianTransfor_K[2] = K_Nodes[2].y() - K_Nodes[0].y();
			JacobianTransfor_K[1] = K_Nodes[1].x() - K_Nodes[0].x() ;
			JacobianTransfor_K[3] = K_Nodes[1].y() - K_Nodes[0].y();
		}
	}
	return JacobianTransfor_K;
}


bool WaveStaggered::FindlocalBasis(const int &m,const Face &Facej, const int &j, const  Cell& K, const std::vector<Node> &K_Nodes ){
	Point Xhat_j =  xToxhat(K,  Facej.getBarryCenter(), K_Nodes) ; 
	std::vector<double>  JacobianTransfor_K_Xhatf = JacobianTransfor_K_X(Xhat_j, K_Nodes );
	double J = fabs( det(JacobianTransfor_K_Xhatf) ) ;
	std::vector<double> PhysicalPsij(_Ndim, 0.0);

	for (int k =0; k < _Ndim ; k++){
		for (int l =0; l < _Ndim ; l++)
			PhysicalPsij[k] += JacobianTransfor_K_Xhatf[k* _Ndim +l] * ReferenceBasisFunctionRaviartThomas(m, Xhat_j, K_Nodes)[l]* Facej.getMeasure()/J ; 
	}
	double cdot =0;
	for (int e=0; e <_Ndim; e++)
		cdot += PhysicalPsij[e] * _vec_sigma.find(j)->second[e]*getOrientation(j,K) ;
	return (abs(cdot -1)< 1e-11); 
}


std::vector<double> WaveStaggered::PhysicalBasisFunctionRaviartThomas(Cell K, int idcell, std::vector<Cell> Support, Face Facej, int j, Point X){
	std::vector<double> xf;
	xf.push_back(X.x());
	if (_Ndim ==2) xf.push_back(X.y());
	std::vector<double> PhysicalPsif_in_X(_Ndim, 0.0); 

	if (_nbTimeStep ==0 && _BasisFunctionAlreadyComputed == false){ 
		std::vector<Node> K_Nodes;
		for (int i=0; i < K.getNodesId().size(); i++)
			K_Nodes.push_back(_mesh.getNode(K.getNodesId()[i]) );
		std::vector<double>  JacobianTransfor_K_Xhat = JacobianTransfor_K_X( xToxhat(K, X, K_Nodes), K_Nodes );
		double J = fabs( det(JacobianTransfor_K_Xhat) );

		bool K_is_in_Support = false;
		for (const auto &cell: Support){
			if (K.x() == cell.x() && K.y() == cell.y())
				K_is_in_Support =true;
		}
		
		for (int  m = 0; m < K.getNumberOfFaces(); m ++){
			if ( FindlocalBasis(m, Facej, j, K,  K_Nodes) == true && K_is_in_Support == true ){
				for (int k =0; k < _Ndim ; k++){	
					for (int l =0; l < _Ndim ; l++)
						PhysicalPsif_in_X[k] += getOrientation(j,K) * Facej.getMeasure()/fabs(J) * JacobianTransfor_K_Xhat[k*_Ndim + l] * ReferenceBasisFunctionRaviartThomas(m, xToxhat(K, X, K_Nodes),  K_Nodes )[l] ;
				}
			}
			
		}
			
		bool already_there =false;
		if (_PhysicalPsif.find(idcell) != _PhysicalPsif.end() ){
			if ( _PhysicalPsif.find(idcell)->second.find(j) != _PhysicalPsif.find(idcell)->second.end()){
				for (const auto & e : _PhysicalPsif.find(idcell)->second.find(j)->second ){
					if ( (_Ndim==1 && e.first[0] == X.x()) || (_Ndim ==2  && e.first[0] == X.x() && e.first[1] == X.y() ))
						already_there = true;
				}
			}
		}
		if (already_there == false)
			_PhysicalPsif[idcell][j].push_back( make_pair(xf, PhysicalPsif_in_X) );
	}
	else{ 
		for (const auto & e :  _PhysicalPsif.find(idcell)->second.find(j)->second ){
			if ( (_Ndim==1 && e.first[0] == X.x()) || (_Ndim ==2  && e.first[0] == X.x() && e.first[1] == X.y()) ){
				for (int k= 0; k < PhysicalPsif_in_X.size(); k++)
					PhysicalPsif_in_X[k]= e.second[k];
			}
		}
	}
	return PhysicalPsif_in_X;
}

double WaveStaggered::det(const std::vector<double> & mat){
	assert(mat.size() == _Ndim*_Ndim);
	double det;
	if (_Ndim ==1) det = mat[0];
	else if (_Ndim ==2) det = mat[0]*mat[3] - mat[1]*mat[2];
	return det; 

}








vector<string> WaveStaggered::getInputFieldsNames(){
	vector<string> result(1);
	result[0]="NOT DEFINED";
	return result;
}
void WaveStaggered::setInputField(const string& nameField, Field& inputField ){}

double WaveStaggered::getTimeStep(){
	return _dt;
}

void WaveStaggered::setExactVelocityFieldAtCells(const Field &atCells){

	_ExactVelocityInftyAtCells = atCells;

	_ExactVelocityInftyAtCells.setName("_ExactVelocityInftyAtCells");
	_time=_ExactVelocityInftyAtCells.getTime();
	_mesh=_ExactVelocityInftyAtCells.getMesh();
	_ExactVelocityInftyAtCells.setInfoOnComponent(0,"ExactVelocityInfty_x(m/s)");
	_ExactVelocityInftyAtCells.setInfoOnComponent(1,"ExactVelocityInfty_y(m/s)");
	string prim(_path+"/");
	prim+=_fileName;
	switch(_saveFormat)
	{
	case VTK :
		_ExactVelocityInftyAtCells.writeVTK(prim+"_ExactVelocityInftyAtCells");
		break;
	}
}


void WaveStaggered::InterpolateFromFacesToCells(std::vector<double> atFaces){ 
	assert( atFaces.size() == _mesh.getNumberOfFaces() ) ;
	Field atCells("ExactVelocity", CELLS, _mesh, 3);
	for (int l=0; l < _Nmailles ; l++){
		for (int k=0; k< 3; k++){
			atCells(l, k) =0;
		}
	}
	for (int i = 0 ; i < _Nfaces ; i++){
		Face Fj = _mesh.getFace(i);
		std::vector< int > idCells = Fj.getCellsId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		double orien1 = getOrientation(i,Ctemp1);
		
		std::vector<double> M1(_Ndim), M2(_Ndim);
		Point xK = Ctemp1.getBarryCenter();
		Point xsigma = Fj.getBarryCenter();

		M1[0] = orien1 * Fj.getMeasure()*(xsigma.x()- xK.x());
		if (_Ndim >1) M1[1] = orien1 * Fj.getMeasure()*(xsigma.y()- xK.y());

		if (Fj.getNumberOfCells() == 2){
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			xK = Ctemp2.getBarryCenter();
			double orien2 = getOrientation(i,Ctemp2);

			M2[0] = orien2 * Fj.getMeasure()*(xsigma.x()- xK.x());
			if (_Ndim >1) M2[1] = orien2 * Fj.getMeasure()*(xsigma.y()- xK.y());
		
			for (int k=0; k< _Ndim; k++){
				atCells(idCells[0], k) += orien1 *  atFaces[i] * M1[k]/Ctemp1.getMeasure(); 
				atCells(idCells[1], k) += orien2 * atFaces[i] * M2[k]/Ctemp2.getMeasure(); 
			}
		}
		else if  (Fj.getNumberOfCells() == 1){
			for (int k=0; k< _Ndim; k++)
				atCells(idCells[0], k) += atFaces[i] * M1[k]/Ctemp1.getMeasure(); 
		}
	}
	string prim(_path+"/");
	string primCells = prim +_fileName + atCells.getName();
	cout << primCells <<endl;

	switch(_saveFormat)
	{
	case VTK :
		atCells.writeVTK(primCells);
		break;
	}
}

double WaveStaggered::ErrorL2VelocityAtFaces(const std::vector<double> &ExactVelocity){
	double errorface =0;
	for (int j=0; j < _Nfaces; j++){
		PetscInt I = _Nmailles + j;
		double InvD_sigma;
		MatGetValues(_InvVol, 1, &I,1, &I, &InvD_sigma);
		double Dsigma = 1/InvD_sigma;
		errorface += Dsigma * (_Velocity(j) - ExactVelocity[j])*(_Velocity(j) - ExactVelocity[j]);
	}
	return errorface;
}

double WaveStaggered::ErrorInftyVelocityBoundary( std::map<int ,double> &BoundaryVelocity ){
	std::map<int, double>::iterator it = BoundaryVelocity.begin();
	double errorboundary =0;
	for (int j=0; j < _Nfaces; j++){
		bool is_j_Interior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
		if ( _mesh.getFace(j).getNumberOfCells() == 1){
			/* cout <<"vel = "<<	 _Velocity(j) << "  "<< BoundaryVelocity.find(j)->second <<endl;
			cout << "diff = "<< _Velocity(j) -BoundaryVelocity.find(j)->second<< endl;
			cout <<" diff^2 = "<< (_Velocity(j) -BoundaryVelocity.find(j)->second)*(_Velocity(j) -BoundaryVelocity.find(j)->second)<<endl; */
			double diff = _Velocity(j) -   BoundaryVelocity.find(j)->second  ;
			errorboundary += pow(diff, 2) ;
			//cout <<"error boundary = "<< errorboundary<<endl;
			}  
	}
	return errorboundary;
}


void WaveStaggered::terminate(){ 
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
	MatDestroy(& _GradDivTilde); 

	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);

	if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_primitiveVars_seq);
}

void WaveStaggered::RelativeEnergyBalanceEq(){
	double flux = 0;
	double residual =0;
	double pb = getboundaryPressure().begin()->second;
	for (int m=0; m<_Nmailles; m++){
		std::vector<int> idFaces = _mesh.getCell(m).getFacesId();
		double divtilde =0;
		for (int f=0; f <idFaces.size(); f++ ){
			if ( _mesh.getFace(idFaces[f]).getNumberOfCells() == 2){
				flux += (_Pressure(m) - pb )* _mesh.getFace(idFaces[f]).getMeasure() * getOrientation(idFaces[f], _mesh.getCell(m)) * _Velocity(idFaces[f]);
			}
			else
				flux += (_Pressure(m) - pb )* _mesh.getFace(idFaces[f]).getMeasure() * getOrientation(idFaces[f], _mesh.getCell(m)) *( _Velocity(idFaces[f]) + getboundaryVelocity().find(idFaces[f])->second)/2.0;
			divtilde += sqrt(_perimeters[m]) * _mesh.getFace(idFaces[f]).getMeasure() * getOrientation(idFaces[f], _mesh.getCell(m)) * (_Velocity(idFaces[f]) - _Velocity_0_Psi(idFaces[f]));
		}
		residual -= pow(divtilde, 2);
	}
	for (int j=0; j<_Nfaces;j++){
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		double grad = 0;
		double velocityjump =0;
		if (Fj.getNumberOfCells() == 2){
			for (int c=0; c< Fj.getNumberOfCells(); c++){
				flux -= (_Velocity(j) - _Velocity_0_Psi(j))*Fj.getMeasure() * getOrientation(j, _mesh.getCell(idCells[c])) * _Pressure(idCells[c]);
				grad +=  sqrt(Fj.getMeasure()) * getOrientation(j, _mesh.getCell(idCells[c])) *( _Pressure(idCells[c]) - pb);
			}
		}
		if (Fj.getNumberOfCells() == 1){
			flux += (_Velocity(j) - _Velocity_0_Psi(j))*Fj.getMeasure() * getOrientation(j, _mesh.getCell(idCells[0])) * (pb- _Pressure(idCells[0]))/2.0 ;
			grad += sqrt(Fj.getMeasure()) *( _Pressure(idCells[0]) -pb ) ;
			velocityjump = sqrt(Fj.getMeasure()) * (_Velocity(j) - _Velocity_0_Psi(j) );
		}
		residual -= pow(grad,2) + pow(velocityjump,2);
	}
	PetscReal norm;
	VecNorm(_newtonVariation, NORM_2, &norm);
	cout << "|| U^n+1 - U^n || = " <<norm <<endl;
	cout << "flux = "<< flux <<"        residual = "<< residual <<",         cfl*|| U^n+1 - U^n ||  ="<< _dt * _c *_maxPerim  / (2 * _minCell ) *_d *  norm<<  ",        sum = "<<residual + _dt * _c *_maxPerim  / (2 * _minCell ) *_d *  norm  <<endl;
}



void WaveStaggered::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/");///Results
	prim+=_fileName;

	RelativeEnergyBalanceEq();
	if (_nbTimeStep >2 ){//&& ( (_Energy.back() -_Energy[_Energy.size() - 2])/_dt) >1e-13
		
		cout <<" Relative Energy( "<< _time <<") = "<<_Energy.back() <<endl;
		cout <<"d_t E =  "<< (_Energy.back() -_Energy[_Energy.size() - 2])/_dt<<endl; 	//std::setprecision(15) << std::fixed<< 
	} 

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
		if (_nbTimeStep == 0) _Velocity_0_Psi_at_Cells.setTime(_time,_nbTimeStep);
		_DivVelocity.setTime(_time,_nbTimeStep);
		for (int l=0; l < _Nmailles ; l++){
			_DivVelocity(l) =0;
			for (int k=0; k< 3; k++){
				_Velocity_at_Cells(l, k) =0;
				if (_nbTimeStep == 0) _Velocity_0_Psi_at_Cells(l, k) =0;
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
			if (_mpi_size > 1)
				VecGetValues(_primitiveVars_seq,1,&I,&_Velocity(i));
			else
				VecGetValues(_primitiveVars,1,&I,&_Velocity(i));

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
					if (_nbTimeStep == 0){
						_Velocity_0_Psi_at_Cells(idCells[0], k) += _Velocity_0_Psi(i) * M1[k]/Ctemp1.getMeasure(); 
						_Velocity_0_Psi_at_Cells(idCells[1], k) += _Velocity_0_Psi(i) * M2[k]/Ctemp2.getMeasure(); 
					}
				}
				_DivVelocity( idCells[0]) += orien1 * Fj.getMeasure() * _Velocity(i)/Ctemp1.getMeasure();
				_DivVelocity( idCells[1]) += orien2 * Fj.getMeasure() * _Velocity(i)/Ctemp2.getMeasure();

			}
			else if  (Fj.getNumberOfCells() == 1){
				for (int k=0; k< _Ndim; k++){
					_Velocity_at_Cells(idCells[0], k) += _Velocity(i) * M1[k]/Ctemp1.getMeasure(); 
					if (_nbTimeStep == 0) _Velocity_0_Psi_at_Cells(idCells[0], k) += _Velocity_0_Psi(i) * M1[k]/Ctemp1.getMeasure();

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
		_Velocity_0_Psi_at_Cells.setInfoOnComponent(0,"x_(m/s)");
		_Velocity_0_Psi_at_Cells.setInfoOnComponent(1,"y_(m/s)");
	
		switch(_saveFormat)
		{
		case VTK :
			_Velocity_at_Cells.writeVTK(prim+"_VelocityAtCells");
			_Velocity_0_Psi_at_Cells.writeVTK(prim+"_Velocity_0_Psi_AtCells");
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
