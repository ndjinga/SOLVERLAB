/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"
#include <numeric>
#include <math.h>

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
	if (_Ndim == 3){	
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!WaveStaggered pas dispo en 3D, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered pas dispo en 3D, arret de calcul");				
	}
				
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
std::map<int,double>  WaveStaggered::getboundaryPressure() const {
		return _boundaryPressure;
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
	return (it !=  _FacePeriodicMap.end());
}
bool  WaveStaggered::IsFaceBoundaryComputedInPeriodic(int j ) {
	return  std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
}
std::map<int,int>  WaveStaggered::getFacePeriodicMap() const{
	return _FacePeriodicMap;
}

void WaveStaggered::setOrientation(int j,std::vector<double> vec_normal_sigma){
	for (int idim = 0; idim < _Ndim; ++idim)
		_vec_sigma[j].push_back(vec_normal_sigma[idim]);
}
double WaveStaggered::getOrientation(int l, Cell Cint) {
	double *vec =new double [_Ndim];	
	vec[0] =0;
	vec[1] =0;	
	std::map<int,int>::iterator it = _FacePeriodicMap.begin();
	while ( ( l !=it->second) && (it !=_FacePeriodicMap.end() ) )it++;
	for(int m=0; m<Cint.getNumberOfFaces(); m++){//we look for l the index of the face Fj for the cell Ctemp1
		if (l == Cint.getFacesId()[m]  ){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);
		}
		else if ((_FacePeriodicMap.find(l) != _FacePeriodicMap.end()) && (_FacePeriodicMap.find(l)->second == Cint.getFacesId()[m])){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);
		}
		else if (it !=_FacePeriodicMap.end() && ( it->first == Cint.getFacesId()[m])){
			for (int idim = 0; idim < _Ndim; ++idim)
				vec[idim] = Cint.getNormalVector(m,idim);

		}
	}
	double dotprod = 0;
	for (int idim = 0; idim < _Ndim; ++idim)
		dotprod += vec[idim] * _vec_sigma.find(l)->second[idim]; 
	delete[] vec;
	return dotprod;

}



void WaveStaggered::setExactVelocityFieldAtCells(const Field &atCells){

	_ExactVelocityInftyAtCells = atCells;

	_ExactVelocityInftyAtCells.setName("_ExactVelocityInftyAtCells");
	_time=_ExactVelocityInftyAtCells.getTime();
	_mesh=_ExactVelocityInftyAtCells.getMesh();
	_ExactVelocityInftyAtCells.setInfoOnComponent(0,"ExactVelocityInfty_x(m/s)");
	_ExactVelocityInftyAtCells.setInfoOnComponent(1,"ExactVelocityInfty_y(m/s)");
	string prim(_path+"/");///Results
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


void WaveStaggered::InterpolateFromFacesToCells(const Field &atFaces, Field &atCells){ 
	assert( atFaces.getTypeOfField() == FACES);
	assert( atCells.getTypeOfField() == CELLS);
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
		double fac;

		if (Ctemp1.getNumberOfFaces() == _Ndim*2)
			fac = 1;
		else if (Ctemp1.getNumberOfFaces() ==  _Ndim + 1)
			fac = -1;

		M1[0] = fac * Fj.getMeasure()*(xsigma.x()- xK.x());
		if (_Ndim >1)
			M1[1] = fac * Fj.getMeasure()*(xsigma.y()- xK.y());

		if (Fj.getNumberOfCells() == 2){
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			xK = Ctemp2.getBarryCenter();
			if (Ctemp2.getNumberOfFaces() == _Ndim*2)
				fac = 1;
			else if (Ctemp2.getNumberOfFaces() ==  _Ndim + 1)
				fac = -1;

			M2[0] = fac * Fj.getMeasure()*(xsigma.x()- xK.x());
			if (_Ndim >1)
				M2[1] = fac * Fj.getMeasure()*(xsigma.y()- xK.y());
		
			for (int k=0; k< _Ndim; k++){
				atCells(idCells[0], k) += atFaces(i) * M1[k]/Ctemp1.getMeasure(); 
				atCells(idCells[1], k) -= atFaces(i) * M2[k]/Ctemp2.getMeasure(); 
			}
		}
		else if  (Fj.getNumberOfCells() == 1){
			for (int k=0; k< _Ndim; k++)
				atCells(idCells[0], k) += atFaces(i) * M1[k]/Ctemp1.getMeasure(); 
		}
	}
	string prim(_path+"/");///Results
	string primCells = prim +_fileName + atCells.getName();
	string primFaces = prim + _fileName  +atFaces.getName();
	cout << primFaces <<endl;
	cout << primCells <<endl;

	switch(_saveFormat)
	{
	case VTK :
		atCells.writeVTK(primCells);
		atFaces.writeVTK(primFaces);
		break;
	}
}

std::vector<double> WaveStaggered::ErrorL2VelocityInfty(const Field &ExactVelocityInftyAtFaces, const Field &ExactVelocityInftyAtCells ){
	double errorface =0;
	double errorcell =0;
	std::vector<double> Error(2);
	for (int j=0; j < _Nfaces; j++){
		PetscInt I = _Nmailles + j;
		double InvD_sigma;
		MatGetValues(_InvVol, 1, &I,1, &I, &InvD_sigma);
		double Dsigma = 1/InvD_sigma;
		errorface += Dsigma * (_Velocity(j) - ExactVelocityInftyAtFaces(j))*(_Velocity(j) - ExactVelocityInftyAtFaces(j));
	}
	for (int j=0; j < _Nmailles; j++){
		double InvK;
		MatGetValues(_InvVol, 1, &j,1, &j, &InvK);
		double K = 1/InvK;
		for (int k =0; k < _Ndim; k++)
			errorcell += K * (_Velocity_at_Cells(j,k) - ExactVelocityInftyAtCells(j,k))*(_Velocity_at_Cells(j,k) - ExactVelocityInftyAtCells(j,k));
	}
	Error[0] = sqrt(errorface);
	Error[1] = sqrt(errorcell);
	return Error;
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
	save();//save initial data
}


double WaveStaggered::MassLumping(const Cell &K, const int &idcell, const Face & Facej, const int &j){
	double masslumping_on_K =0;
	std::vector<Cell> Support_f, Support_j ;
	Support_f.push_back(K);
	Support_j.push_back(K);

	for (int l =0; l <K.getNumberOfFaces(); l ++){
		Node vertex1, vertex2;
		vertex1 = _mesh.getNode( _mesh.getFace( K.getFacesId()[l] ).getNodesId()[0] );
		Point Xl =  _mesh.getFace( K.getFacesId()[l] ).getBarryCenter();
		if (_Ndim ==2 )  vertex2 = _mesh.getNode( _mesh.getFace( K.getFacesId()[l] ).getNodesId()[1] );
		double Area = (_Ndim==2 )? fabs((K.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (K.y() - vertex1.y() )* (vertex2.x() - vertex1.x() ) )/2.0 : abs(K.x() - vertex1.x()) ;
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

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;	
		
		double mlump = MassLumping(Ctemp1, idCells[0], Fj, j);
		InvPerimeter1 = (_Ndim ==2) ? 1.0/_perimeters(idCells[0]) : 1.0 ;
		InvVol1 = 1.0/Ctemp1.getMeasure();
		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end())
				idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) ? _mesh.getFace(_FacePeriodicMap.find(j)->second ) :  Fj;
			InvPerimeter2 = (_Ndim ==2) ? 1.0/_perimeters(idCells[1]) : 1.0 ;
			mlump += MassLumping(Ctemp2,idCells[1], Fj_physical, j);

			InvD_sigma = 1.0/mlump;
			PetscScalar InvVol2 = 1.0/Ctemp2.getMeasure();
			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1 , INSERT_VALUES );
			MatSetValue(_InvVol, idCells[1],idCells[1], InvVol2, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES); 	
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvSurface,idCells[1],idCells[1], InvPerimeter2, INSERT_VALUES );
		}
		else if (IsWallBound || IsSteggerBound ) { 
			InvD_sigma = 1.0/mlump; 	
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES);  
		}
	}
	//TODo à supprimer : à permis de comprendre d'où venait le probleme de métrique sur masslumping
	/* for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
 		Cell Ctemp1 = _mesh.getCell(idCells[0]);

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;	
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;			
	
		if (  (Fj.x() < -1e-10 ) ){
			Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) ? _mesh.getFace(_FacePeriodicMap.find(j)->second ) :  Fj;
			if (IsInterior){	
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
				double mlump = MassLumping(Ctemp2, idCells[1], Fj, j)+ MassLumping(Ctemp1, idCells[0], Fj, j);

				for (int k=0; k <_Nfaces ; k++){
					Face Fk = _mesh.getFace(k);
					std::vector< int > idCellsFk = Fk.getCellsId();
 					Cell C1 = _mesh.getCell(idCellsFk[0]);
					if ( fabs( fabs(Fk.x()) -fabs(Fj.x())) <1e-10 && fabs( Fk.y()-Fj.y()) < 1e-10 && Fk.x() > 1e-10 && std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),k ) != _InteriorFaceSet.end() && Fk.x() > 1e-10){
						Cell C2 = _mesh.getCell(idCellsFk[1]);
						if (fabs( mlump -(MassLumping(C1, idCellsFk[0], Fk, k)+ MassLumping(C2, idCellsFk[1], Fk, k))) >1e-10){
							cout << "mlump = " <<mlump <<" and mlump wrong = "<< MassLumping(C1, idCellsFk[0], Fk, k)+ MassLumping(C2, idCellsFk[1], Fk, k)<< " k = " <<k <<" Fk = ("<< Fk.x() <<", "<< Fk.y()<< ") theta ="<< atan(Fk.y()/Fk.x())<<endl;
							cout <<"k = "<< k <<" theta = "<< atan(Fk.y()/Fk.x())<<endl; 
						}
					}
				}
			}
			if (IsWallBound || IsSteggerBound){
				double mlump = MassLumping(Ctemp1, idCells[0], Fj, j);
				for (int k=0; k <_Nfaces ; k++){
					Face Fk = _mesh.getFace(k);
					std::vector< int > idCellsFk = Fk.getCellsId();
 					Cell C1 = _mesh.getCell(idCellsFk[0]);
					if ( fabs( fabs(Fk.x()) -fabs(Fj.x())) <1e-10 && fabs( Fk.y()-Fj.y()) < 1e-10 && Fk.x() > 1e-10 ){
						if (fabs( mlump -(MassLumping(C1, idCellsFk[0], Fk, k)) >1e-10) ){
							cout << "mlump = " <<mlump <<" and mlump wrong = "<< MassLumping(C1, idCellsFk[0], Fk, k)<< " k = " <<k <<" Fk = ("<< Fk.x() <<", "<< Fk.y()<< ") theta ="<< atan(Fk.y()/Fk.x())<<endl;
							cout <<"k = "<< k <<" theta = "<< atan(Fk.y()/Fk.x())<<endl; 
						}
					}
				}
			}
		}
	} */
	
	MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(_InvSurface,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
}


/* void WaveStaggered::AssembleMetricsMatrices(){
	MatZeroEntries(_InvVol); 
	MatZeroEntries(_InvSurface); 
	for (int j=0; j<_Nfaces;j++){ 
		Face Fj = _mesh.getFace(j);
		std::vector< int > idCells = Fj.getCellsId();
		std::vector< int > NodesFj =  Fj.getNodesId();
		Cell Ctemp1 = _mesh.getCell(idCells[0]);
		std::vector< int > idFaces = Ctemp1.getFacesId();
		PetscScalar det, InvD_sigma, InvVol1, InvVol2, InvPerimeter1, InvPerimeter2;
		PetscInt IndexFace = _Nmailles + j;

		bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;    
		bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
		bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;    
		
		InvPerimeter1 = (_Ndim ==2) ? 1.0/_perimeters(idCells[0]) : 1.0 ;
		InvVol1 = 1.0/Ctemp1.getMeasure();
		if (IsInterior){
			if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end())
					idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
			Cell Ctemp2 = _mesh.getCell(idCells[1]);
			Face Fj_physical =  ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end() ) ? _mesh.getFace(_FacePeriodicMap.find(j)->second ) :  Fj;
			InvPerimeter2 = (_Ndim ==2) ? 1.0/(_perimeters(idCells[1]) ) : 1.0;
			InvVol2 = 1.0/Ctemp2.getMeasure();
			if (_Ndim == 1) det = Ctemp2.x() - Ctemp1.x();
			if (_Ndim ==2){
				std::vector<int> nodes =  Fj.getNodesId();
				Node vertex = _mesh.getNode( nodes[0] );
				// determinant of the vectors forming the diamond cell around the face sigma
				det = (Ctemp1.x() - vertex.x() )* (Ctemp2.y() - vertex.y() ) - (Ctemp1.y() - vertex.y() )* (Ctemp2.x() - vertex.x() );
			}
			InvD_sigma = 1.0/PetscAbsReal(det);

			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1 , INSERT_VALUES );
			MatSetValue(_InvVol, idCells[1],idCells[1], InvVol2, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES);         
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvSurface,idCells[1],idCells[1], InvPerimeter2, INSERT_VALUES );
		}
		else if (IsWallBound || IsSteggerBound ) { 
			if (_Ndim == 1)
					det = InvD_sigma = 2.0/Ctemp1.getMeasure() ;
			if (_Ndim == 2){
				std::vector< int > nodes =  Fj.getNodesId();
				Node vertex1 = _mesh.getNode( nodes[0] );
				Node vertex2 = _mesh.getNode( nodes[1] );
				det = (Ctemp1.x() - vertex1.x() )* (vertex2.y() - vertex1.y() ) - (Ctemp1.y() - vertex1.y() )* (vertex2.x() - vertex1.x() );
				// determinant of the vectors forming the interior half diamond cell around the face sigma
			}
			InvD_sigma = 1.0/PetscAbsReal(det);          
			MatSetValue(_InvSurface,idCells[0],idCells[0], InvPerimeter1, INSERT_VALUES );
			MatSetValue(_InvVol, idCells[0],idCells[0], InvVol1, INSERT_VALUES );
			MatSetValue(_InvVol, IndexFace, IndexFace,  InvD_sigma, INSERT_VALUES);

		}
		MatAssemblyBegin(_InvVol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvVol, MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(_InvSurface,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_InvSurface, MAT_FINAL_ASSEMBLY);
	}
} */



double WaveStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme
	//The matrices are assembled only in the first time step since linear problem
	MatZeroEntries(_LaplacianPressure);
	MatZeroEntries(_Div); 
	MatZeroEntries(_DivTranspose); 

	if (_timeScheme == Explicit ){ 
		if ( _nbTimeStep == 0  ){   
			// Assembly of matrices 
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				std::vector< int > idCells = Fj.getCellsId();
				Cell Ctemp1 = _mesh.getCell(idCells[0]);

				bool IsInterior = std::find(_InteriorFaceSet.begin(), _InteriorFaceSet.end(),j ) != _InteriorFaceSet.end() ;
				bool IsWallBound = std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(),j ) != _WallBoundFaceSet.end() ;
				bool IsSteggerBound = std::find(_SteggerBoundFaceSet.begin(), _SteggerBoundFaceSet.end(),j ) != _SteggerBoundFaceSet.end() ;

				// Metrics
				double orien = getOrientation(j,Ctemp1);
				PetscScalar orientedFaceArea = orien * Fj.getMeasure();
				PetscScalar orientedMinusFaceArea = -orientedFaceArea;
				PetscScalar FaceArea = Fj.getMeasure();
				PetscScalar MinusFaceArea = -FaceArea;
				PetscScalar det, InvPerimeter1, InvPerimeter2, InvD_sigma, InvVol1, InvVol2;
				PetscInt IndexFace = _Nmailles + j;
				InvVol1 = 1.0/(Ctemp1.getMeasure()*Ctemp1.getNumberOfFaces());
				
				if ( IsInterior ){	
					if ( _FacePeriodicMap.find(j) != _FacePeriodicMap.end()  )
						idCells.push_back( _mesh.getFace(_FacePeriodicMap.find(j)->second).getCellsId()[0]  );
					Cell Ctemp2 = _mesh.getCell(idCells[1]);
					
					/******************* Pressure equation ***********************/
					MatSetValue(_Div, idCells[0], j, orientedFaceArea, ADD_VALUES ); 
					MatSetValue(_Div, idCells[1], j, orientedMinusFaceArea, ADD_VALUES );  
					MatSetValue(_LaplacianPressure, idCells[0], idCells[0], MinusFaceArea, ADD_VALUES ); 
					MatSetValue(_LaplacianPressure, idCells[0], idCells[1], FaceArea, ADD_VALUES );  
					MatSetValue(_LaplacianPressure, idCells[1], idCells[1], MinusFaceArea, ADD_VALUES ); 
					MatSetValue(_LaplacianPressure, idCells[1], idCells[0], FaceArea, ADD_VALUES );  
					
					/******************* Velocity equation ***********************/
					MatSetValue(_DivTranspose, j, idCells[0], orientedFaceArea, ADD_VALUES ); 
					MatSetValue(_DivTranspose, j, idCells[1], orientedMinusFaceArea, ADD_VALUES );				
				
				}
				else if (IsSteggerBound || IsWallBound ) { // && (periodicFaceNotComputed == false) if boundary face and face index is different from periodic faces not computed 		

					/***************** Pressure equation related matrices ******************/
					MatSetValue(_Div, idCells[0], j, orientedFaceArea, ADD_VALUES ); 
					MatSetValue(_LaplacianPressure, idCells[0], idCells[0], MinusFaceArea, ADD_VALUES );

					//Is the face a wall boundarycondition face
					if (IsWallBound ){
						VecGetValues(_primitiveVars,1,&idCells[0],&_pInt);
						_pExt =  Fj.getMeasure()*_pInt; //_pExt = pin so (grad p)_j = 0
					}
					else if (IsSteggerBound){ //Imposed boundaryconditions
						std::map<int,double> boundaryPressure = getboundaryPressure(); 
						std::map<int,double>::iterator it = boundaryPressure.find(j);
						_pExt = Fj.getMeasure()*boundaryPressure[it->first]; 
					}
					VecSetValues(_BoundaryTerms, 1,&idCells[0], &_pExt, INSERT_VALUES );
				}	
			}
			
			MatAssemblyBegin(_Div,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_Div, MAT_FINAL_ASSEMBLY);
			MatAssemblyBegin(_DivTranspose, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_DivTranspose, MAT_FINAL_ASSEMBLY);

			MatAssemblyBegin(_LaplacianPressure,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(_LaplacianPressure, MAT_FINAL_ASSEMBLY);
			VecAssemblyBegin(_BoundaryTerms);
			VecAssemblyEnd(_BoundaryTerms);
			VecScale(_BoundaryTerms, _d * _c);

			// _A = (dc _LaplacianPressure  ;  -1/rho B         )
			//      (kappa B^t     ;  dc -B^t(1/|dK|) B ) 
			MatScale(_DivTranspose, -1.0);
			if (_nbTimeStep == 0) 
				MatMatMatMult(_DivTranspose,_InvSurface, _Div , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &_GradDivTilde); 
			MatScale(_LaplacianPressure, _d*_c );
			MatScale(_Div, -1.0/_rho);
			MatScale(_DivTranspose, -1.0*_kappa);
			MatScale(_GradDivTilde, _d*_c); 
			Mat G[4], Prod;
			G[0] = _LaplacianPressure;
			G[1] = _Div;
			G[2] = _DivTranspose;
			G[3] = _GradDivTilde; 
			
			// WARNING MatCreateNest() memory problem when time evolution of matrices
			assert( _nbTimeStep == 0);
			MatCreateNest(PETSC_COMM_WORLD,2, NULL, 2, NULL , G, &_A); 
			MatConvert(_A, MATAIJ, MAT_INPLACE_MATRIX, & _A);
			MatMatMult(_InvVol, _A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & Prod); 
			MatCopy(Prod,_A, SAME_NONZERO_PATTERN); 
			MatDestroy(& Prod);

			ComputeMinCellMaxPerim();

			
		}
		if (_isWall && _nbTimeStep >0 ){	
			for (int j=0; j<_Nfaces;j++){
				Face Fj = _mesh.getFace(j);
				if (Fj.getNumberOfCells()==1) { //if boundary face 
					//Is the face a wall boundarycondition face
					if (std::find(_WallBoundFaceSet.begin(), _WallBoundFaceSet.end(), j)!=_WallBoundFaceSet.end()){
						std::vector< int > idCells = Fj.getCellsId();
						VecGetValues(_primitiveVars,1,&idCells[0],&_pInt);
						_pExt = _d * _c * Fj.getMeasure()*_pInt; //_pExt = pin so (grad p)_j = 0
						VecSetValues(_BoundaryTerms, 1,&idCells[0], &_pExt, INSERT_VALUES );
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
		VecDestroy(& Prod2); 	

	}
	
	ComputeEnergyAtTimeT();
	double dt = _cfl * _minCell / (_maxPerim * _c ) ;
	double PreviousTime = _Time.back();
	_Time.push_back(PreviousTime+ dt);
	return dt ;
}



void WaveStaggered::ComputeMinCellMaxPerim(){
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
		}
		if ( K_Nodes.size() == 3){
			if (i ==0){      
				Psihat[0] = Xhat.x();
				Psihat[1] = Xhat.y() - 1;
			}
			else if (i ==1){      
				Psihat[0] = Xhat.x();// * 2  * sqrt(2.0);
				Psihat[1] = Xhat.y();// * 2 * sqrt(2.0);
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
	double J = Jacobian(JacobianTransfor_K_Xhatf);
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

	if (_nbTimeStep ==0){ 
		std::vector<Node> K_Nodes;
		for (int i=0; i < K.getNodesId().size(); i++)
			K_Nodes.push_back(_mesh.getNode(K.getNodesId()[i]) );
		std::vector<double>  JacobianTransfor_K_Xhat = JacobianTransfor_K_X( xToxhat(K, X, K_Nodes), K_Nodes );
		double J = Jacobian(JacobianTransfor_K_Xhat);

		bool K_is_in_Support = false;
		for (const auto &cell: Support){
			if (K.getBarryCenter().x() == cell.getBarryCenter().x() && K.getBarryCenter().y() == cell.getBarryCenter().y())
				K_is_in_Support =true;
		}
		
		for (int  m = 0; m < K.getNumberOfFaces(); m ++){
			if ( FindlocalBasis(m, Facej, j, K,  K_Nodes) == true && K_is_in_Support == true ){
				for (int k =0; k < _Ndim ; k++){	
					for (int l =0; l < _Ndim ; l++)
						PhysicalPsif_in_X[k] += getOrientation(j,K) * Facej.getMeasure()/J * JacobianTransfor_K_Xhat[k*_Ndim + l] * ReferenceBasisFunctionRaviartThomas(m, xToxhat(K, X, K_Nodes),  K_Nodes )[l] ;
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

//TODO or fabs( ) ?
double WaveStaggered::Jacobian(const std::vector<double> & mat){
	assert(mat.size() == _Ndim*_Ndim);
	double jacobian;
	if (_Ndim ==1) jacobian = mat[0];
	else if (_Ndim ==2) jacobian = mat[0]*mat[3] - mat[1]*mat[2];
	return jacobian; 

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
	//converged=convergence des iterations
	if(_timeScheme == Explicit)
		converged=true;
	VecAXPY(_primitiveVars, 1, _newtonVariation);//Vk+1=Vk+relaxation*deltaV
	
	return converged;

}

void WaveStaggered::computeNewtonVariation(){
	if(_timeScheme == Explicit){	
		VecCopy(_b,_newtonVariation);
		VecScale(_newtonVariation, _dt);
		if(_verbose && (_nbTimeStep-1)%_freqSave ==0){
			cout<<"Vecteur _newtonVariation =_b*dt"<<endl;
			VecView(_newtonVariation,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
	}
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

void WaveStaggered::setPeriodicFaces( Mesh &M, const char &Direction, int ncells){ 
	for (int j=0;j<M.getNumberOfFaces() ; j++){
		Face my_face=M.getFace(j);
		double e;
		double tol = 1.0/(ncells *4);
		if (Direction == 'x')
			e=my_face.x();
		else if (Direction == 'y' && _Ndim ==2)
			e=my_face.y();
		if (my_face.getNumberOfCells() ==1 ){ 
			if( (_Ndim==2) &&  e>tol && e< (1.0-tol) ){
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

vector<string> WaveStaggered::getInputFieldsNames(){
	vector<string> result(1);
	result[0]="NOT DEFINED";
	return result;
}
void WaveStaggered::setInputField(const string& nameField, Field& inputField ){}

double WaveStaggered::getTimeStep(){
	return _dt;
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

void WaveStaggered::save(){
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


			//Multiply by -1 on triangles (why ?)
			M1[0] =  orien1 * Fj.getMeasure()*(xsigma.x()- xK.x()) * ((_Ndim ==2) ? ( ( Ctemp1.getNumberOfFaces() == 3 ) ? -1 : 1.0  ): 1.0);
			if (_Ndim >1) M1[1] = orien1 * Fj.getMeasure()*(xsigma.y()- xK.y()) * ((_Ndim ==2) ? ( ( Ctemp1.getNumberOfFaces() == 3 ) ? -1 : 1.0  ): 1.0);

			if (Fj.getNumberOfCells() == 2){
				Cell Ctemp2 = _mesh.getCell(idCells[1]);
				double orien2 = getOrientation(i,Ctemp2);
				Point xK = Ctemp2.getBarryCenter();

				M2[0] = orien2 * Fj.getMeasure()*(xsigma.x()- xK.x()) *  ((_Ndim ==2) ? ( ( Ctemp2.getNumberOfFaces() == 3 ) ? -1 : 1.0  ): 1.0);
				if (_Ndim >1) M2[1] = orien2 * Fj.getMeasure()*(xsigma.y()- xK.y()) * ((_Ndim ==2) ? ( ( Ctemp2.getNumberOfFaces() == 3 ) ? -1 : 1.0  ): 1.0);
			
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
			_Velocity_at_Cells.writeVTK(prim+"_Velocity at cells");
			_Velocity.writeVTK(prim+"_Velocity");
			_DivVelocity.writeVTK(prim+"Divergence Velocity");
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
