#include "LinearElasticityModel.hxx"
#include "SparseMatrixPetsc.hxx"
#include "Node.hxx"
#include "math.h"
#include <algorithm> 
#include <fstream>
#include <sstream>

using namespace std;

LinearElasticityModel::LinearElasticityModel(int dim, bool FECalculation,  double rho, double lambda, double mu,MPI_Comm comm){
	/* Initialisation of PETSC */
	//check if PETSC is already initialised
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
	{//check if MPI is already initialised
		int mpiInitialized;
		MPI_Initialized(&mpiInitialized);
		if(mpiInitialized)
			PETSC_COMM_WORLD = comm;
		PetscInitialize(NULL,NULL,0,0);//Note this is ok if MPI has been been initialised independently from PETSC

#if CMAKE_BUILD_TYPE==DEBUG
		int argc = 2;
		char **argv = new char*[argc];
		argv[0] = (char*)"LinearElasticityModel";
		argv[1] = (char*)"-on_error_attach_debugger";
		PetscInitialize(&argc, &argv, 0, 0);//Note this is ok if MPI has been been initialised independently from PETSC
#else
		PetscInitialize(NULL,NULL,0,0);//Note this is ok if MPI has been been initialised independently from PETSC
#endif
	}
	MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&_size);
	PetscPrintf(PETSC_COMM_WORLD,"Simulation on %d processors\n",_size);//Prints to standard out, only from the first processor in the communicator. Calls from other processes are ignored. 
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Processor [%d] ready for action\n",_rank);//Prints synchronized output from several processors. Output of the first processor is followed by that of the second, etc. 
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    if(lambda < 0.)
    {
        std::cout<<"First Lamé coefficient="<<lambda<<endl;
        throw CdmathException("Error : First Lamé coefficient lambda cannot  be negative");
    }
    if(2*mu+dim*lambda < 0.)
    {
        std::cout<<"First Lamé coefficient="<<lambda<<", second Lamé coefficient="<<mu<<", 2*mu+dim*lambda= "<<2*mu+dim*lambda<<endl;
        throw CdmathException("Error : 2*mu+dim*lambda cannot  be negative");
    }
    if(dim<=0)
    {
        std::cout<<"space dimension="<<dim<<endl;
        throw CdmathException("Error : parameter dim cannot  be negative");
    }

    _FECalculation=FECalculation;
    _onlyNeumannBC=false;    
    
	_Ndim=dim;
	_nVar=dim;
	_initializedMemory=false;

    //Mesh data
    _neibMaxNbCells=0;    
    _meshSet=false;
    _neibMaxNbNodes=0;    

    //Boundary conditions
    _boundaryNodeIds=std::vector< int >(0);
    _dirichletNodeIds=std::vector< int >(0);
    _NboundaryNodes=0;
    _NdirichletNodes=0;
    _NunknownNodes=0;
    _dirichletValuesSet=false;   
    _neumannValuesSet=false;   
    
    //Linear solver data
	_precision=1.e-6;
	_precision_Newton=_precision;
	_MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
	_maxPetscIts=50;
	int _PetscIts=0;//the number of iterations of the linear solver
	_ksptype = (char*)&KSPGMRES;
	_pctype = (char*)&PCLU;
	_conditionNumber=false;
	_erreur_rel= 0;

    //parameters for monitoring simulation
	_verbose = false;
	_system = false;
	_runLogFile=new ofstream;

	//result save parameters
	_fileName = "LinearElasticityProblem";
	char result[ PATH_MAX ];//extracting current directory
	getcwd(result, PATH_MAX );
	_path=string( result );
	_saveFormat=VTK;
    _computationCompletedSuccessfully=false;
    
    //heat transfer parameters
	_lambda= lambda;
	_mu    = mu;
	_rho   = rho;
	_densityFieldSet=false;
}

void LinearElasticityModel::initialize()
{
	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation

	if(!_meshSet)
		throw CdmathException("LinearElasticityModel::initialize() set mesh first");
	else
    {
		cout<<"!!!! Initialisation of the computation of the elastic deformation of a solid using ";
        *_runLogFile<<"!!!!! Initialisation of the computation of the elastic deformation of a solid using ";
        if(!_FECalculation)
        {
            cout<< "Finite volumes method"<<endl<<endl;
            *_runLogFile<< "Finite volumes method"<<endl<<endl;
        }
        else
        {
            cout<< "Finite elements method"<<endl<<endl;
            *_runLogFile<< "Finite elements method"<<endl<<endl;
        }
    }
	/**************** Field creation *********************/

	if(!_densityFieldSet){
        if(_FECalculation){
            _densityField=Field("Density",NODES,_mesh,1);
            for(int i =0; i<_Nnodes; i++)
                _densityField(i) = _rho;
        }
        else{
            _densityField=Field("Density",CELLS,_mesh,1);
            for(int i =0; i<_Nmailles; i++)
                _densityField(i) = _rho;
        }
        _densityFieldSet=true;
    }

    /* Détection des noeuds frontière avec une condition limite de Dirichlet */
    if(_FECalculation)
    {
        if(_NboundaryNodes==_Nnodes)
            cout<<"!!!!! Warning : all nodes are boundary nodes !!!!!"<<endl<<endl;

        for(int i=0; i<_NboundaryNodes; i++)
        {
            std::map<int,double>::iterator it=_dirichletBoundaryValues.find(_boundaryNodeIds[i]);
            if( it != _dirichletBoundaryValues.end() )
                _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
            else if( _mesh.getNode(_boundaryNodeIds[i]).getGroupNames().size()==0 )
            {
                cout<<"!!! No boundary value set for boundary node" << _boundaryNodeIds[i]<< endl;
                *_runLogFile<< "!!! No boundary value set for boundary node" << _boundaryNodeIds[i]<<endl;
                _runLogFile->close();
                throw CdmathException("Missing boundary value");
            }
            else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType==NoneBCLinearElasticity)
            {
                cout<<"!!! No boundary condition set for boundary node " << _boundaryNodeIds[i]<< endl;
                *_runLogFile<< "!!!No boundary condition set for boundary node " << _boundaryNodeIds[i]<<endl;
                _runLogFile->close();
                throw CdmathException("Missing boundary condition");
            }
            else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType==DirichletLinearElasticity)
                _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
            else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType!=NeumannLinearElasticity)
            {
                cout<<"!!! Wrong boundary condition "<< _limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType<< " set for boundary node " << _boundaryNodeIds[i]<< endl;
                cout<<"!!! Accepted boundary conditions are Dirichlet "<< DirichletLinearElasticity <<" and Neumann "<< NeumannLinearElasticity << endl;
                *_runLogFile<< "Wrong boundary condition "<< _limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType<< " set for boundary node " << _boundaryNodeIds[i]<<endl;
                *_runLogFile<< "Accepted boundary conditions are Dirichlet "<< DirichletLinearElasticity <<" and Neumann "<< NeumannLinearElasticity <<endl;
                _runLogFile->close();
                throw CdmathException("Wrong boundary condition");
            }
        }	
        _NdirichletNodes=_dirichletNodeIds.size();
        _NunknownNodes=_Nnodes - _NdirichletNodes;
        cout<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
		*_runLogFile<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
    }

	//creation de la matrice
    if(!_FECalculation)
        MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles*_nVar, _Nmailles*_nVar, (1+_neibMaxNbCells), NULL, &_A);
    else
        MatCreateSeqAIJ(PETSC_COMM_SELF, _NunknownNodes*_nVar, _NunknownNodes*_nVar, (1+_neibMaxNbNodes), NULL, &_A);

	VecCreate(PETSC_COMM_SELF, &_displacements);

	VecDuplicate(_displacements, &_b);//RHS of the linear system

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

    //Checking whether all boundary conditions are Neumann boundary condition
    if(!_neumannValuesSet)//Boundary conditions set via LimitField structure
    {
        map<string, LimitFieldLinearElasticity>::iterator it = _limitField.begin();
        while(it != _limitField.end() and (it->second).bcType == NeumannLinearElasticity)
            it++;
        _onlyNeumannBC = (it == _limitField.end() && _limitField.size()>0);//what if _limitField.size()==0 ???
    }
    else
        if(_FECalculation)
            _onlyNeumannBC = _neumannBoundaryValues.size()==_NboundaryNodes;
        else
            _onlyNeumannBC = _neumannBoundaryValues.size()==_mesh.getBoundaryFaceIds().size();

    //If only Neumann BC, then matrix is singular and solution should be sought in space of mean zero vectors
    if(_onlyNeumannBC)
    {
        std::cout<<"## Warning all boundary conditions are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
        std::cout<<"## As a consequence we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;
        *_runLogFile<<"## Warning all boundary condition are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
        *_runLogFile<<"## As a consequence we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;

		//Check that the matrix is symmetric
		PetscBool isSymetric;
		MatIsSymmetric(_A,_precision,&isSymetric);
		if(!isSymetric)
			{
				cout<<"Singular matrix is not symmetric, tolerance= "<< _precision<<endl;
				throw CdmathException("Singular matrix should be symmetric with kernel composed of constant vectors");
			}
		MatNullSpace nullsp;
		MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
		MatSetNullSpace(_A, nullsp);
		MatSetTransposeNullSpace(_A, nullsp);
		MatNullSpaceDestroy(&nullsp);
		//PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);
		//PCFactorSetShiftAmount(_pc,1e-10);
    }

	_initializedMemory=true;

}

double LinearElasticityModel::computeStiffnessMatrix(bool & stop)
{
    double result;
    
    if(_FECalculation)
        result=computeStiffnessMatrixFE(stop);
    else
        result=computeStiffnessMatrixFV(stop);

    if(_verbose or _system)
        MatView(_A,PETSC_VIEWER_STDOUT_SELF);

    return  result;
}

double LinearElasticityModel::computeStiffnessMatrixFE(bool & stop){
	Cell Cj;
	string nameOfGroup;
	double dn, coeff;
	MatZeroEntries(_A);
	VecZeroEntries(_b);
    
    Matrix M(_Ndim+1,_Ndim+1);//cell geometry matrix
    std::vector< Vector > GradShapeFuncs(_Ndim+1);//shape functions of cell nodes
    std::vector< int > nodeIds(_Ndim+1);//cell node Ids
    std::vector< Node > nodes(_Ndim+1);//cell nodes
    int i_int, j_int; //index of nodes j and k considered as unknown nodes
    bool dirichletCell_treated;
    
    std::vector< vector< double > > values(_Ndim+1,vector< double >(_Ndim+1,0));//values of shape functions on cell node
    for (int idim=0; idim<_Ndim+1;idim++)
        values[idim][idim]=1;

    /* parameters for boundary treatment */
    vector< Vector > valuesBorder(_Ndim+1);
    Vector GradShapeFuncBorder(_Ndim+1);
    
	for (int j=0; j<_Nmailles;j++)
    {
		Cj = _mesh.getCell(j);

        for (int idim=0; idim<_Ndim+1;idim++){
            nodeIds[idim]=Cj.getNodeId(idim);
            nodes[idim]=_mesh.getNode(nodeIds[idim]);
            for (int jdim=0; jdim<_Ndim;jdim++)
                M(idim,jdim)=nodes[idim].getPoint()[jdim];
            M(idim,_Ndim)=1;
        }
        for (int idim=0; idim<_Ndim+1;idim++)
            GradShapeFuncs[idim]=DiffusionEquation::gradientNodal(M,values[idim])/DiffusionEquation::fact(_Ndim);

        /* Loop on the edges of the cell */
        for (int idim=0; idim<_Ndim+1;idim++)
        {
            if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[idim])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[idim])
            {//First node of the edge is not Dirichlet node
                i_int=DiffusionEquation::unknownNodeIndex(nodeIds[idim], _dirichletNodeIds);//assumes Dirichlet boundary node numbering is strictly increasing
                dirichletCell_treated=false;
                for (int jdim=0; jdim<_Ndim+1;jdim++)
                {
                    if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[jdim])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[jdim])
                    {//Second node of the edge is not Dirichlet node
                        j_int= DiffusionEquation::unknownNodeIndex(nodeIds[jdim], _dirichletNodeIds);//assumes Dirichlet boundary node numbering is strictly increasing
                        MatSetValue(_A,i_int,j_int,_conductivity*(_DiffusionTensor*GradShapeFuncs[idim])*GradShapeFuncs[jdim]/Cj.getMeasure(), ADD_VALUES);
                    }
                    else if (!dirichletCell_treated)
                    {//Second node of the edge is a Dirichlet node
                        dirichletCell_treated=true;
                        for (int kdim=0; kdim<_Ndim+1;kdim++)
                        {
							std::map<int,double>::iterator it=_dirichletBoundaryValues.find(nodeIds[kdim]);
							if( it != _dirichletBoundaryValues.end() )
                            {
                                if( _dirichletValuesSet )//New way of storing BC
                                    valuesBorder[kdim]=_dirichletBoundaryValues[it->second];
                                else    //old way of storing BC
                                    valuesBorder[kdim]=_limitField[_mesh.getNode(nodeIds[kdim]).getGroupName()].displacement;
                            }
                            else
                                valuesBorder[kdim]=Vector(_Ndim);                            
                        }
                        GradShapeFuncBorder=DiffusionEquation::gradientNodal(M,valuesBorder)/DiffusionEquation::fact(_Ndim);
                        double coeff =-_conductivity*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
                        VecSetValue(_b,i_int,coeff, ADD_VALUES);                        
                    }
                }
            }
        }            
	}
    
    //Calcul de la contribution de la condition limite de Neumann au second membre
    if( _NdirichletNodes !=_NboundaryNodes)
    {
        vector< int > boundaryFaces = _mesh.getBoundaryFaceIds();
        int NboundaryFaces=boundaryFaces.size();
        for(int i = 0; i< NboundaryFaces ; i++)//On parcourt les faces du bord
        {
            Face Fi = _mesh.getFace(boundaryFaces[i]);
            for(int j = 0 ; j<_Ndim ; j++)//On parcourt les noeuds de la face
            {
                if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),Fi.getNodeId(j))==_dirichletNodeIds.end())//node j is an Neumann BC node (not a Dirichlet BC node)
                {
                    j_int=DiffusionEquation::unknownNodeIndex(Fi.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
                    if( _neumannValuesSet )
                        coeff =Fi.getMeasure()/(_Ndim+1)*_neumannBoundaryValues[Fi.getNodeId(j)];
                    else    
                        coeff =Fi.getMeasure()/(_Ndim+1)*_limitField[_mesh.getNode(Fi.getNodeId(j)).getGroupName()].normalForce;
                    VecSetValue(_b, j_int, coeff, ADD_VALUES);
                }
            }
        }
    }

    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);

    stop=false ;

	return INFINITY;
}

double LinearElasticityModel::computeStiffnessMatrixFV(bool & stop){
	long nbFaces = _mesh.getNumberOfFaces();
	Face Fj;
	Cell Cell1,Cell2;
	string nameOfGroup;
	double inv_dxi, inv_dxj;
	double barycenterDistance,coeff;
	Vector normale(_Ndim);
	double dn;
	PetscInt idm, idn;
	std::vector< int > idCells;
	MatZeroEntries(_A);
	VecZeroEntries(_b);

	for (int j=0; j<nbFaces;j++){
		Fj = _mesh.getFace(j);

		// compute the normal vector corresponding to face j : from idCells[0] to idCells[1]
		idCells = Fj.getCellsId();
		Cell1 = _mesh.getCell(idCells[0]);
		idm = idCells[0];
        for(int l=0; l<Cell1.getNumberOfFaces(); l++){
            if (j == Cell1.getFacesId()[l]){
                for (int idim = 0; idim < _Ndim; ++idim)
                    normale[idim] = Cell1.getNormalVector(l,idim);
                break;
            }
        }

		//Compute velocity at the face Fj
		dn=_conductivity*(_DiffusionTensor*normale)*normale;

		// compute 1/dxi = volume of Ci/area of Fj
        inv_dxi = Fj.getMeasure()/Cell1.getMeasure();

		// If Fj is on the boundary
		if (Fj.getNumberOfCells()==1) {
			if(_verbose )
			{
				cout << "face numero " << j << " cellule frontiere " << idCells[0] << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << normale[p] << ",";
				cout << ") "<<endl;
			}

            std::map<int,double>::iterator it=_dirichletBoundaryValues.find(j);
            if( it != _dirichletBoundaryValues.end() )
            {
                barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
                MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                                     , ADD_VALUES);
                VecSetValue(_b,idm,    dn*inv_dxi/barycenterDistance*it->second, ADD_VALUES);
            }
            else
            {
                nameOfGroup = Fj.getGroupName();
    
                if (_limitField[nameOfGroup].bcType==NeumannLinearElasticity){
                    if( _neumannValuesSet )
                        coeff =Fi.getMeasure()*_neumannBoundaryValues[j];
                    else    
                        coeff =Fi.getMeasure()*_limitField[nameOfGroup].normalForce;
                    VecSetValue(_b,idm,    coeff, ADD_VALUES);
                }
                else if(_limitField[nameOfGroup].bcType==DirichletLinearElasticity){
                    barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
                    MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                           , ADD_VALUES);
                    VecSetValue(_b,idm,    dn*inv_dxi/barycenterDistance*_limitField[nameOfGroup].displacement, ADD_VALUES);
                }
                else {
                    stop=true ;
                    cout<<"!!!!!!!!!!!!!!! Error LinearElasticityModel::computeStiffnessMatrixFV !!!!!!!!!!"<<endl;
                    cout<<"!!!!!! Boundary condition not accepted for boundary named !!!!!!!!!!"<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
                    cout<<"Accepted boundary conditions are Neumann "<<NeumannLinearElasticity<< " and Dirichlet "<<DirichletLinearElasticity<<endl;
                    *_runLogFile<<"!!!!!! Boundary condition not accepted for boundary named !!!!!!!!!!"<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
                    _runLogFile->close();
                    throw CdmathException("Boundary condition not accepted");
                }
            }
			// if Fj is inside the domain
		} else 	if (Fj.getNumberOfCells()==2 ){
			if(_verbose )
			{
				cout << "face numero " << j << " cellule gauche " << idCells[0] << " cellule droite " << idCells[1];
				cout << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << normale[p] << ",";
				cout << ") "<<endl;
			}
			Cell2 = _mesh.getCell(idCells[1]);
			idn = idCells[1];
			if (_Ndim > 1)
				inv_dxj = Fj.getMeasure()/Cell2.getMeasure();
			else
				inv_dxj = 1/Cell2.getMeasure();
			
			barycenterDistance=Cell1.getBarryCenter().distance(Cell2.getBarryCenter());

			MatSetValue(_A,idm,idm, dn*inv_dxi/barycenterDistance, ADD_VALUES);
			MatSetValue(_A,idm,idn,-dn*inv_dxi/barycenterDistance, ADD_VALUES);
			MatSetValue(_A,idn,idn, dn*inv_dxj/barycenterDistance, ADD_VALUES);
			MatSetValue(_A,idn,idm,-dn*inv_dxj/barycenterDistance, ADD_VALUES);
		}
		else
        {
            *_runLogFile<<"LinearElasticityModel::computeStiffnessMatrixFV(): incompatible number of cells around a face"<<endl;
			throw CdmathException("LinearElasticityModel::computeStiffnessMatrixFV(): incompatible number of cells around a face");
        }
	}

	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);
    
    stop=false ;
	
    return INFINITY;
}

double LinearElasticityModel::computeRHS(bool & stop)//Contribution of the PDE RHS to the linear systemm RHS (boundary conditions do contribute to the system RHS via the function computeStiffnessMatrix)
{
	VecAssemblyBegin(_b);

    if(!_FECalculation)
        for (int i=0; i<_Nmailles;i++)
			for (int j=0; j<_Ndim;j++)
				VecSetValue(_b,i*nVar+j,_gravity(j)*_densityField(i),ADD_VALUES);
    else
        {
            Cell Ci;
            std::vector< int > nodesId;
            for (int i=0; i<_Nmailles;i++)
            {
                Ci=_mesh.getCell(i);
                nodesId=Ci.getNodesId();
                for (int j=0; j<nodesId.size();j++)
                    if(!_mesh.isBorderNode(nodesId[j]))
						for (int k=0; k<_Ndim; k++)
							VecSetValue(_b,DiffusionEquation::unknownNodeIndex(nodesId[j], _dirichletNodeIds)*nVar+k, _gravity(k)*_densityField(j)*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
            }
        }
    
	VecAssemblyEnd(_b);

    if(_verbose or _system)
        VecView(_b,PETSC_VIEWER_STDOUT_SELF);

    stop=false ;
	return INFINITY;
}

bool LinearElasticityModel::solveLinearSystem()
{
	bool resolutionOK;

    //Only implicit scheme considered
    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);

#if PETSC_VERSION_GREATER_3_5
    KSPSetOperators(_ksp, _A, _A);
#else
    KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

    if(_conditionNumber)
        KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
    KSPSolve(_ksp, _b, _displacements);

	KSPConvergedReason reason;
	KSPGetConvergedReason(_ksp,&reason);
    KSPGetIterationNumber(_ksp, &_PetscIts);
    double residu;
    KSPGetResidualNorm(_ksp,&residu);
	if (reason!=2 and reason!=3)
    {
        cout<<"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc."<<endl;
        cout<<"!!!!!!!!!!!!! Itérations maximales "<<_maxPetscIts<<" atteintes, résidu="<<residu<<", précision demandée= "<<_precision<<endl;
        cout<<"Solver used "<<  _ksptype<<", preconditioner "<<_pctype<<", Final number of iteration= "<<_PetscIts<<endl;
		*_runLogFile<<"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc."<<endl;
        *_runLogFile<<"!!!!!!!!!!!!! Itérations maximales "<<_maxPetscIts<<" atteintes, résidu="<<residu<<", précision demandée= "<<_precision<<endl;
        *_runLogFile<<"Solver used "<<  _ksptype<<", preconditioner "<<_pctype<<", Final number of iteration= "<<_PetscIts<<endl;
		_runLogFile->close();
		
		resolutionOK = false;
    }
    else{
        if( _MaxIterLinearSolver < _PetscIts)
            _MaxIterLinearSolver = _PetscIts;
        cout<<"## Système linéaire résolu en "<<_PetscIts<<" itérations par le solveur "<<  _ksptype<<" et le preconditioneur "<<_pctype<<", précision demandée= "<<_precision<<endl<<endl;
		*_runLogFile<<"## Système linéaire résolu en "<<_PetscIts<<" itérations par le solveur "<<  _ksptype<<" et le preconditioneur "<<_pctype<<", précision demandée= "<<_precision<<endl<<endl;
        
        resolutionOK = true;
    }

	return resolutionOK;
}

void LinearElasticityModel::setMesh(const Mesh &M)
{
	if(_Ndim != M.getSpaceDimension() or _Ndim!=M.getMeshDimension())//for the moment we must have space dim=mesh dim
	{
        cout<< "Problem : dimension defined is "<<_Ndim<< " but mesh dimension= "<<M.getMeshDimension()<<", and space dimension is "<<M.getSpaceDimension()<<endl;
		*_runLogFile<< "Problem : dim = "<<_Ndim<< " but mesh dim= "<<M.getMeshDimension()<<", mesh space dim= "<<M.getSpaceDimension()<<endl;
		*_runLogFile<<"LinearElasticityModel::setMesh: mesh has incorrect dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("LinearElasticityModel::setMesh: mesh has incorrect  dimension");
	}

	_mesh=M;
	_Nmailles = _mesh.getNumberOfCells();
	_Nnodes =   _mesh.getNumberOfNodes();
    
    cout<<"Mesh has "<< _Nmailles << " cells and " << _Nnodes << " nodes"<<endl<<endl;;
	*_runLogFile<<"Mesh has "<< _Nmailles << " cells and " << _Nnodes << " nodes"<<endl<<endl;
    
	// find  maximum nb of neibourghs
    if(!_FECalculation)
    {
    	_VV=Field ("Displacements", CELLS, _mesh, _Ndim);
        _neibMaxNbCells=_mesh.getMaxNbNeighbours(CELLS);
    }
    else
    {
        if(_Ndim==1 )//The 1D cdmath mesh is necessarily made of segments
			cout<<"1D Finite element method on segments"<<endl;
        else if(_Ndim==2)
        {
			if( _mesh.isTriangular() )//Mesh dim=2
				cout<<"2D Finite element method on triangles"<<endl;
			else if (_mesh.getMeshDimension()==1)//Mesh dim=1
				cout<<"1D Finite element method on a 2D network : space dimension is "<<_Ndim<< ", mesh dimension is "<<_mesh.getMeshDimension()<<endl;			
			else
			{
				cout<<"Error Finite element with Space dimension "<<_Ndim<< ", and mesh dimension  "<<_mesh.getMeshDimension()<< ", mesh should be either triangular either 1D network"<<endl;
				*_runLogFile<<"LinearElasticityModel::setMesh: mesh has incorrect dimension"<<endl;
				_runLogFile->close();
				throw CdmathException("LinearElasticityModel::setMesh: mesh has incorrect cell types");
			}
        }
        else if(_Ndim==3)
        {
			if( _mesh.isTetrahedral() )//Mesh dim=3
				cout<<"3D Finite element method on tetrahedra"<<endl;
			else if (_mesh.getMeshDimension()==2 and _mesh.isTriangular())//Mesh dim=2
				cout<<"2D Finite element method on a 3D surface : space dimension is "<<_Ndim<< ", mesh dimension is "<<_mesh.getMeshDimension()<<endl;			
			else if (_mesh.getMeshDimension()==1)//Mesh dim=1
				cout<<"1D Finite element method on a 3D network : space dimension is "<<_Ndim<< ", mesh dimension is "<<_mesh.getMeshDimension()<<endl;			
			else
			{
				cout<<"Error Finite element with Space dimension "<<_Ndim<< ", and mesh dimension  "<<_mesh.getMeshDimension()<< ", mesh should be either tetrahedral, either a triangularised surface or 1D network"<<endl;
				*_runLogFile<<"LinearElasticityModel::setMesh: mesh has incorrect dimension"<<endl;
				_runLogFile->close();
				throw CdmathException("LinearElasticityModel::setMesh: mesh has incorrect cell types");
			}
        }

		_VV=Field ("Temperature", NODES, _mesh, _Ndim);

        _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
        _boundaryNodeIds = _mesh.getBoundaryNodeIds();
        _NboundaryNodes=_boundaryNodeIds.size();
    }

	_meshSet=true;
}

void LinearElasticityModel::setLinearSolver(linearSolver kspType, preconditioner pcType)
{
	//_maxPetscIts=maxIterationsPetsc;
	// set linear solver algorithm
	if (kspType==GMRES)
		_ksptype = (char*)&KSPGMRES;
	else if (kspType==CG)
		_ksptype = (char*)&KSPCG;
	else if (kspType==BCGS)
		_ksptype = (char*)&KSPBCGS;
	else {
		cout << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
		*_runLogFile << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'GMRES', 'CG' or 'BCGS' algorithm is acceptable !!!");
	}
	// set preconditioner
	if (pcType == NONE)
		_pctype = (char*)&PCNONE;
	else if (pcType ==LU)
		_pctype = (char*)&PCLU;
	else if (pcType == ILU)
		_pctype = (char*)&PCILU;
	else if (pcType ==CHOLESKY)
		_pctype = (char*)&PCCHOLESKY;
	else if (pcType == ICC)
		_pctype = (char*)&PCICC;
	else {
		cout << "!!! Error : only 'NONE', 'LU', 'ILU', 'CHOLESKY' or 'ICC' preconditioners are acceptable !!!" << endl;
		*_runLogFile << "!!! Error : only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
		_runLogFile->close();
		throw CdmathException("!!! Error : only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" );
	}
}

bool LinearElasticityModel::solveStationaryProblem()
{
	if(!_initializedMemory)
	{
		*_runLogFile<< "ProblemCoreFlows::run() call initialize() first"<< _fileName<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::run() call initialize() first");
	}
	bool stop=false; // Does the Problem want to stop (error) ?

	cout<< "!!! Running test case "<< _fileName << " using ";
	*_runLogFile<< "!!! Running test case "<< _fileName<< " using ";

    if(!_FECalculation)
    {
        cout<< "Finite volumes method"<<endl<<endl;
		*_runLogFile<< "Finite volumes method"<<endl<<endl;
	}
    else
	{
        cout<< "Finite elements method"<<endl<<endl;
		*_runLogFile<< "Finite elements method"<< endl<<endl;
	}

    computeStiffnessMatrix( stop);
    if (stop){
        cout << "Error : failed computing diffusion matrix, stopping calculation"<< endl;
        *_runLogFile << "Error : failed computing diffusion matrix, stopping calculation"<< endl;
 		_runLogFile->close();
       throw CdmathException("Failed computing diffusion matrix");
    }
    computeRHS(stop);
    if (stop){
        cout << "Error : failed computing right hand side, stopping calculation"<< endl;
        *_runLogFile << "Error : failed computing right hand side, stopping calculation"<< endl;
        throw CdmathException("Failed computing right hand side");
    }
    stop = !solveLinearSystem();
    if (stop){
        cout << "Error : failed solving linear system, stopping calculation"<< endl;
        *_runLogFile << "Error : failed linear system, stopping calculation"<< endl;
		_runLogFile->close();
        throw CdmathException("Failed solving linear system");
    }
    
    _computationCompletedSuccessfully=true;
    save();
    
    *_runLogFile<< "!!!!!! Computation successful"<< endl;
	_runLogFile->close();

	return !stop;
}

void LinearElasticityModel::save(){
    cout<< "Saving numerical results"<<endl<<endl;
    *_runLogFile<< "Saving numerical results"<< endl<<endl;

	string resultFile(_path+"/LinearElasticityModel");//Results

	resultFile+="_";
	resultFile+=_fileName;

	// create mesh and component info
    string suppress ="rm -rf "+resultFile+"_*";
    system(suppress.c_str());//Nettoyage des précédents calculs identiques
    
    if(_verbose or _system)
        VecView(_displacements,PETSC_VIEWER_STDOUT_SELF);

    double uk; 
    if(!_FECalculation)
        for(int i=0; i<_Nmailles; i++)
            {
				for(int j=0; j<_nVar; j++)
				{
					int k=i*_nVar+j;
					VecGetValues(_displacements, 1, &k, &uk);
					_VV(i,j)=uk;
				}
            }
    else
    {
        int globalIndex;
        for(int i=0; i<_NunknownNodes; i++)
        {
			globalIndex = DiffusionEquation::globalNodeIndex(i, _dirichletNodeIds);
			for(int j=0; j<_nVar; j++)
			{
				int k=i*_nVar+j;
				VecGetValues(_displacements, 1, &k, &uk);
				_VV(globalIndex,j)=uk;
			}
        }

        Node Ni;
        string nameOfGroup;
        for(int i=0; i<_NdirichletNodes; i++)
        {
            Ni=_mesh.getNode(_dirichletNodeIds[i]);
            nameOfGroup = Ni.getGroupName();
			for(int j=0; j<_nVar; j++)
				_VV(_dirichletNodeIds[i])=_limitField[nameOfGroup].displacement[i];
        }
    }

    switch(_saveFormat)
    {
        case VTK :
            _VV.writeVTK(resultFile);
            break;
        case MED :
            _VV.writeMED(resultFile);
            break;
        case CSV :
            _VV.writeCSV(resultFile);
            break;
    }
}
Field 
LinearElasticityModel::getOutputDisplacementField()
{
    if(!_computationCompletedSuccessfully)
        throw("Computation not performed yet or failed. No displacement field available");
    else
        return _VV;
}

void LinearElasticityModel::terminate()
{
	VecDestroy(&_displacements);
	VecDestroy(&_b);
	MatDestroy(&_A);
}
void 
LinearElasticityModel::setDirichletValues(map< int, double> dirichletBoundaryValues)
{
    _dirichletValuesSet=true;
    _dirichletBoundaryValues=dirichletBoundaryValues;
}
void 
LinearElasticityModel::setNeumannValues(map< int, double> neumannBoundaryValues)
{
    _neumannValuesSet=true;
    _neumannBoundaryValues=neumannBoundaryValues;
}

double 
LinearElasticityModel::getConditionNumber(bool isSingular, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getConditionNumber( isSingular, tol);
}
