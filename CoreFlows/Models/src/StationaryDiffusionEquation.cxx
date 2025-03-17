#include "StationaryDiffusionEquation.hxx"
#include "DiffusionEquation.hxx"
#include "Node.hxx"
#include "SparseMatrixPetsc.hxx"
#include "math.h"
#include <algorithm> 
#include <fstream>
#include <sstream>

using namespace std;

StationaryDiffusionEquation::StationaryDiffusionEquation(int dim, bool FECalculation, double lambda,MPI_Comm comm){
    /* Initialisation of PETSC */
    //check if PETSC is already initialised
    PetscBool petscInitialized;
    PetscInitialized(&petscInitialized);
    if(!petscInitialized)
    {//check if MPI is already initialised
        int mpiInitialized;
        MPI_Initialized(&mpiInitialized);
        if(mpiInitialized)
            PETSC_COMM_WORLD = comm;// run PETSc on ONLY a subset of MPI_COMM_WORLD

#if CMAKE_BUILD_TYPE==DEBUG
        int argc = 2;
        char **argv = new char*[argc];
        argv[0] = (char*)"StationaryDiffusionEquation";
        argv[1] = (char*)"-on_error_attach_debugger";
        PetscInitialize(&argc, &argv, 0, 0);//Note this is ok if MPI has been been initialised independently from PETSC
#else
        PetscInitialize(NULL,NULL,0,0);//Note this is ok if MPI has been been initialised independently from PETSC
#endif
    }
    MPI_Comm_rank(PETSC_COMM_WORLD,&_mpi_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_mpi_size);
    PetscPrintf(PETSC_COMM_WORLD,"\n Simulation on %d processors\n",_mpi_size);//Prints to standard out, only from the first processor in the communicator. Calls from other processes are ignored. 
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Processor [%d] ready for action\n",_mpi_rank);//Prints synchronized output from several processors. Output of the first processor is followed by that of the second, etc. 
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    if(lambda < 0.)
    {
        std::cout<<"Conductivity="<<lambda<<endl;
        throw CdmathException("!!!!!!!!Error : conductivity parameter lambda cannot  be negative");
    }
    if(dim<=0)
    {
        std::cout<<"Space dimension="<<dim<<endl;
        throw CdmathException("!!!!!!!!Error : parameter dim cannot  be negative");
    }

    _FECalculation=FECalculation;
    _onlyNeumannBC=false;    
    
    _Ndim=dim;
    _nVar=1;//scalar prolem
    _dt_src=0;
    _diffusionMatrixSet=false;
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
    
    //Linear solver data
    _precision=1.e-6;
    _precision_Newton=_precision;
    _MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
    _maxPetscIts=50;
    _maxNewtonIts=50;
    _NEWTON_its=0;
    _PetscIts=0;//the number of iterations of the linear solver
    _ksptype = (char*)&KSPGMRES;
    _pctype = (char*)&PCILU;

    _conditionNumber=false;
    _erreur_rel= 0;

    //parameters for monitoring simulation
    _verbose = false;
    _system = false;
    _runLogFile=new ofstream;

    //result save parameters
    _fileName = "StationaryDiffusionProblem";
    char result[ PATH_MAX ];//extracting current directory
    getcwd(result, PATH_MAX );
    _path=string( result );
    _saveFormat=VTK;
    _computationCompletedSuccessfully=false;
    
    //heat transfer parameters
    _conductivity=lambda;
    _fluidTemperatureFieldSet=false;
    _fluidTemperature=0;
    _heatPowerFieldSet=false;
    _heatTransfertCoeff=0;
    _heatSource=0;

    /* Default diffusion tensor is diagonal */
    _DiffusionTensor=Matrix(_Ndim);
    for(int idim=0;idim<_Ndim;idim++)
        _DiffusionTensor(idim,idim)=_conductivity;

    PetscPrintf(PETSC_COMM_WORLD,"\n Stationary diffusion problem with conductivity %.2f", lambda);
    if(FECalculation)
        PetscPrintf(PETSC_COMM_WORLD," and finite elements method\n\n");
    else
        PetscPrintf(PETSC_COMM_WORLD," and finite volumes method\n\n");
}

void StationaryDiffusionEquation::initialize()
{
    if(_mpi_rank==0)
    {
        /**************** Affichages par le proc 0 (on évite affichages redondants) *********************/
    
        _runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation
    
        if(!_meshSet)
            throw CdmathException("!!!!!!!!StationaryDiffusionEquation::initialize() set mesh first");
        else
        {
            cout<<"\n Initialisation of the computation of the temperature diffusion in a solid using ";
            *_runLogFile<<"\n Initialisation of the computation of the temperature diffusion in a solid using ";
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
        
        /**************** Field creation par le proc 0 (le seul à connaitre le maillage) *********************/
    
        if(!_heatPowerFieldSet){
            _heatPowerField=Field("Heat power",_VV.getTypeOfField(),_mesh,1);
            for(int i =0; i<_VV.getNumberOfElements(); i++)
                _heatPowerField(i) = _heatSource;
            _heatPowerFieldSet=true;
        }
        if(!_fluidTemperatureFieldSet){
            _fluidTemperatureField=Field("Fluid temperature",_VV.getTypeOfField(),_mesh,1);
            for(int i =0; i<_VV.getNumberOfElements(); i++)
                _fluidTemperatureField(i) = _fluidTemperature;
            _fluidTemperatureFieldSet=true;
        }
    
        /******** Détection par le proc 0 (le seul à connaitre le maillage) des noeuds frontière avec une condition limite de Dirichlet ********/
        if(_FECalculation)
        {
            if(_NboundaryNodes==_Nnodes)
                cout<<"!!!!! Warning : all nodes are boundary nodes !!!!!"<<endl<<endl;
    
            for(int i=0; i<_NboundaryNodes; i++)
            {
                std::map<int,double>::iterator it_dirichlet=_dirichletBoundaryValues.find(_boundaryNodeIds[i]);
                std::map<int,double>::iterator it_neumann=_neumannBoundaryValues.find(_boundaryNodeIds[i]);
                if( it_dirichlet != _dirichletBoundaryValues.end() || it_neumann!=_neumannBoundaryValues.end())
                {
                    if( it_dirichlet != _dirichletBoundaryValues.end() )
                        _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
                }
                else if( _mesh.getNode(_boundaryNodeIds[i]).getGroupNames().size()==0 )
                {
                    cout<<"!!! No boundary group set for boundary node" << _boundaryNodeIds[i]<< endl;
                    *_runLogFile<< "!!! No boundary group set for boundary node" << _boundaryNodeIds[i]<<endl;
                    _runLogFile->close();
                    throw CdmathException("Missing boundary group");
                }
                else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType==NoneBCStationaryDiffusion)
                {
                    cout<<"!!! No boundary condition set for boundary node " << _boundaryNodeIds[i]<< endl;
                    cout<<"!!! Accepted boundary conditions are DirichletStationaryDiffusion "<< DirichletStationaryDiffusion <<" and NeumannStationaryDiffusion "<< NeumannStationaryDiffusion << endl;
                    *_runLogFile<< "!!! No boundary condition set for boundary node " << _boundaryNodeIds[i]<<endl;
                    *_runLogFile<< "!!! Accepted boundary conditions are DirichletStationaryDiffusion "<< DirichletStationaryDiffusion <<" and NeumannStationaryDiffusion "<< NeumannStationaryDiffusion <<endl;
                    _runLogFile->close();
                    throw CdmathException("Missing boundary condition");
                }
                else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType==DirichletStationaryDiffusion)
                    _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
                else if(_limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType!=NeumannStationaryDiffusion)
                {
                    cout<<"!!! Wrong boundary condition "<< _limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType<< " set for boundary node " << _boundaryNodeIds[i]<< endl;
                    cout<<"!!! Accepted boundary conditions are DirichletStationaryDiffusion "<< DirichletStationaryDiffusion <<" and NeumannStationaryDiffusion "<< NeumannStationaryDiffusion << endl;
                    *_runLogFile<< "!!! Wrong boundary condition "<< _limitField[_mesh.getNode(_boundaryNodeIds[i]).getGroupName()].bcType<< " set for boundary node " << _boundaryNodeIds[i]<<endl;
                    *_runLogFile<< "!!! Accepted boundary conditions are DirichletStationaryDiffusion "<< DirichletStationaryDiffusion <<" and NeumannStationaryDiffusion "<< NeumannStationaryDiffusion <<endl;
                    _runLogFile->close();
                    throw CdmathException("Wrong boundary condition");
                }
            }    
            _NdirichletNodes=_dirichletNodeIds.size();
            _NunknownNodes=_Nnodes - _NdirichletNodes;
            cout<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
            *_runLogFile<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
        }
    }
    
    if(!_FECalculation)
        _globalNbUnknowns = _Nmailles*_nVar;
    else{
        MPI_Bcast(&_NunknownNodes, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        _globalNbUnknowns = _NunknownNodes*_nVar;
    }
    
    /* Parallel vectors creations (all procs) */
    VecCreate(PETSC_COMM_WORLD, &_Tk);//main unknown
    VecSetSizes(_Tk,PETSC_DECIDE,_globalNbUnknowns);
    VecSetFromOptions(_Tk);
    VecGetLocalSize(_Tk, &_localNbUnknowns);
    
    VecDuplicate(_Tk, &_Tkm1);
    VecDuplicate(_Tk, &_deltaT);
    VecDuplicate(_Tk, &_b);//RHS of the linear system: _b=Tn/dt + _b0 + puisance volumique + couplage thermique avec le fluide

    /* Parallel matrix creation (all procs) */
       MatCreateAIJ(PETSC_COMM_WORLD, _localNbUnknowns, _localNbUnknowns, _globalNbUnknowns, _globalNbUnknowns, _d_nnz, NULL, _o_nnz, NULL, &_A);
    
    /* Local sequential vector creation (all procs) */
    if(_mpi_size>1 && _mpi_rank == 0)
        VecCreateSeq(PETSC_COMM_SELF,_globalNbUnknowns,&_Tk_seq);//For saving results on proc 0
    if(_mpi_size==0)
        _Tk_seq=_Tk;
    VecScatterCreateToZero(_Tk,&_scat,&_Tk_seq);

    /* PETSc Linear solver (all procs) */
    KSPCreate(PETSC_COMM_WORLD, &_ksp);
    KSPSetType(_ksp, _ksptype);
    KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
    KSPGetPC(_ksp, &_pc);
    //PETSc preconditioner
    if(_mpi_size==1 )
        PCSetType(_pc, _pctype);
    else
    {
        PCSetType(_pc, PCBJACOBI);//Global preconditioner is block jacobi
        if(_pctype != (char*)&PCILU)//Default pc type is ilu
        {
            PetscOptionsSetValue(NULL,"-sub_pc_type ",_pctype);
            PetscOptionsSetValue(NULL,"-sub_ksp_type ","preonly");
            //If the above setvalue does not work, try the following
            /*
            KSPSetUp(_ksp);//to set the block Jacobi data structures (including creation of an internal KSP context for each block)
            KSP * subKSP;
            PC subpc;
            int nlocal;//nb local blocs (should equal 1)
            PCBJacobiGetSubKSP(_pc,&nlocal,NULL,&subKSP);
            if(nlocal==1)
            {
                KSPSetType(subKSP[0], KSPPREONLY);//local block solver is same as global
                KSPGetPC(subKSP[0],&subpc);
                PCSetType(subpc,_pctype);
            }
            else
                throw CdmathException("PC Block Jacobi, more than one block in this processor!!");
            */ 
        }
    }

    //Checking whether at least one boundary conditions is imposed
    if( _limitField.size()==0 && _neumannBoundaryValues.size()==0 && _dirichletBoundaryValues.size()==0 && _dirichletBoundaryField.getMEDCouplingField()==NULL && _neumannBoundaryField.getMEDCouplingField()==NULL )
        throw CdmathException("No boundary condition imposed. Cannot initialize simulation.");
        
    //Checking on all procs whether all boundary conditions are Neumann boundary condition ->singular system
    if(_limitField.size()!=0)//Boundary conditions set via LimitField structure
    {
        map<string, LimitFieldStationaryDiffusion>::iterator it = _limitField.begin();
        while(it != _limitField.end() and (it->second).bcType == NeumannStationaryDiffusion)
            it++;
        _onlyNeumannBC = (it == _limitField.end());
    }
    else//Boundary conditions set via boundary values
        _onlyNeumannBC = _dirichletBoundaryValues.size()==0;

    //If only Neumann BC, then matrix is singular and solution should be sought in space of mean zero vectors
    if(_onlyNeumannBC)
    {
        if(_mpi_rank==0)//Avoid redundant printing
        {
            std::cout<<"### Warning : all boundary conditions are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
            std::cout<<"### Check the compatibility condition between the right hand side and the boundary data. For homogeneous Neumann BCs, the right hand side must have integral equal to zero."<<std::endl;
            std::cout<<"### The system matrix being singular, we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;
            *_runLogFile<<"### Warning : all boundary condition are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
            *_runLogFile<<"### The system matrix being singular, we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;
            *_runLogFile<<"### Check the compatibility condition between the right hand side and the boundary data. For homogeneous Neumann BCs, the right hand side must have integral equal to zero."<<std::endl;
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

double StationaryDiffusionEquation::computeTimeStep(bool & stop){
    if(!_diffusionMatrixSet)//The diffusion matrix is computed once and for all time steps
        computeDiffusionMatrix(stop);

    _dt_src=computeRHS(stop);
    stop=false;
    return _dt_src;
}

double StationaryDiffusionEquation::computeDiffusionMatrix(bool & stop)
{
    double result;
    
    MatZeroEntries(_A);
    VecZeroEntries(_b);

    if(_FECalculation)
        result=computeDiffusionMatrixFE(stop);
    else
        result=computeDiffusionMatrixFV(stop);

    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);

    //Contribution from the solid/fluid heat exchange with assumption of constant heat transfer coefficient
    //update value here if variable  heat transfer coefficient
    if(_heatTransfertCoeff>_precision)
        MatShift(_A,_heatTransfertCoeff);//Contribution from the liquit/solid heat transfer
        
    if(_verbose or _system)
        MatView(_A,PETSC_VIEWER_STDOUT_WORLD);

    return  result;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFE(bool & stop){
    if(_mpi_rank == 0)
        {
        Cell Cj;
        string nameOfGroup;
        double coeff;//Diffusion coefficients between nodes i and j (to be inserted in stiffness matrix), or boundary coefficient (to be inserted in RHS)
        
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
        vector< double > valuesBorder(_Ndim+1);
        Vector GradShapeFuncBorder(_Ndim+1);
        std::map<int,double>::iterator it;
        
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
                        {//Second node of the edge is not Dirichlet node -> contribution to the stiffness matrix
                            j_int= DiffusionEquation::unknownNodeIndex(nodeIds[jdim], _dirichletNodeIds);//assumes Dirichlet boundary node numbering is strictly increasing
                            coeff = (_DiffusionTensor*GradShapeFuncs[idim])*GradShapeFuncs[jdim]/Cj.getMeasure();
#if CMAKE_BUILD_TYPE==DEBUG
                            if( idim != jdim && coeff>0)//non acute triangle/tetrahedron -> violation of mawimum principle
                                PetscPrintf(PETSC_COMM_WORLD,"\n !!! Warning : non acute triangle/tetrahedron cell %d, nodes %d and %d, coeff=%.2f, possible violation of the maximum principle \n",j, nodeIds[idim], nodeIds[jdim], coeff);
#endif
                            MatSetValue(_A,i_int,j_int, coeff, ADD_VALUES);
                        }
                        else if (!dirichletCell_treated)
                        {//Second node of the edge is a Dirichlet node -> contribution to the right hand side
                            dirichletCell_treated=true;
                            for (int kdim=0; kdim<_Ndim+1;kdim++)
                            {
                                if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[kdim])!=_dirichletNodeIds.end())//node kdim is a Dirichlet BC node
                                {
                                    it=_dirichletBoundaryValues.find(nodeIds[kdim]);
                                    if( it != _dirichletBoundaryValues.end() )//Une valeur limite est associée au noeud//BC set via setDirichletValues
                                        valuesBorder[kdim]=it->second;
                                    else//Une valeur limite est associée au groupe frontière//BC set via setDirichletBoundaryCondition    
                                        valuesBorder[kdim]=_limitField[_mesh.getNode(nodeIds[kdim]).getGroupName()].T;
                                }
                                else
                                    valuesBorder[kdim]=0;                      
                            }
                            GradShapeFuncBorder=DiffusionEquation::gradientNodal(M,valuesBorder)/DiffusionEquation::fact(_Ndim);
                            coeff =-1.*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
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
                    if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),Fi.getNodeId(j))==_dirichletNodeIds.end())//node j is a Neumann BC node (not a Dirichlet BC node)
                    {
                        j_int=DiffusionEquation::unknownNodeIndex(Fi.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
                        if( _neumannBoundaryValues.size()!=0 )//Une valeur limite est associée à chaque noeud frontière
                            coeff =Fi.getMeasure()/_Ndim*_neumannBoundaryValues[Fi.getNodeId(j)];
                        else//Une valeur limite est associée à chaque groupe frontière
                            coeff =Fi.getMeasure()/_Ndim*_limitField[_mesh.getNode(Fi.getNodeId(j)).getGroupName()].normalFlux;
                        VecSetValue(_b, j_int, coeff, ADD_VALUES);
                    }
                }
            }
        }
    }

    if(_onlyNeumannBC)    //Check that the matrix is symmetric
    {
        PetscBool isSymetric;
        MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
        MatIsSymmetric(_A,_precision,&isSymetric);
        if(!isSymetric)
        {
            PetscPrintf(PETSC_COMM_WORLD,"\n Singular matrix is not symmetric, tolerance= %.2f \n",_precision);
            throw CdmathException("Singular matrix should be symmetric with kernel composed of constant vectors");
        }
    }

    _diffusionMatrixSet=true;
    stop=false ;

    return INFINITY;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFV(bool & stop){
    if(_mpi_rank == 0)
    {
        long nbFaces = _mesh.getNumberOfFaces();
        Face Fj;
        Cell Cell1,Cell2;
        string nameOfGroup;
        double inv_dxi, inv_dxj;
        double barycenterDistance;
        Vector normale(_Ndim);
        double dn;
        PetscInt idm, idn;
        std::vector< int > idCells;
    
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
            dn=(_DiffusionTensor*normale)*normale;
    
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
                if( it != _dirichletBoundaryValues.end() )//Une valeur limite est associée à la face frontière
                {
                    barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
                    MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                                     , ADD_VALUES);
                    VecSetValue(_b,idm,    dn*inv_dxi/barycenterDistance*it->second, ADD_VALUES);
                }
                else//Une valeur limite est associée au groupe frontière
                {
                    nameOfGroup = Fj.getGroupName();
        
                    if (_limitField[nameOfGroup].bcType==NeumannStationaryDiffusion){
                        VecSetValue(_b,idm,    -dn*inv_dxi*_limitField[nameOfGroup].normalFlux, ADD_VALUES);
                    }
                    else if(_limitField[nameOfGroup].bcType==DirichletStationaryDiffusion){
                        barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
                        MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                           , ADD_VALUES);
                        VecSetValue(_b,idm,    dn*inv_dxi/barycenterDistance*_limitField[nameOfGroup].T, ADD_VALUES);
                    }
                    else {
                        stop=true ;
                        cout<<"!!!!!!!!!!!!!!! Error StationaryDiffusionEquation::computeDiffusionMatrixFV !!!!!!!!!!"<<endl;
                        cout<<"!!!!!! No boundary condition set for boundary named "<<nameOfGroup<< "!!!!!!!!!! _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
                        cout<<"Accepted boundary conditions are NeumannStationaryDiffusion "<<NeumannStationaryDiffusion<< " and DirichletStationaryDiffusion "<<DirichletStationaryDiffusion<<endl;
                        *_runLogFile<<"!!!!!! Boundary condition not set for boundary named "<<nameOfGroup<< "!!!!!!!!!! _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
                        _runLogFile->close();
                        throw CdmathException("Boundary condition not set");
                    }
                }
                // if Fj is inside the domain
            } else     if (Fj.getNumberOfCells()==2 ){
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
                *_runLogFile<<"StationaryDiffusionEquation::computeDiffusionMatrixFV(): incompatible number of cells around a face"<<endl;
                throw CdmathException("StationaryDiffusionEquation::computeDiffusionMatrixFV(): incompatible number of cells around a face");
            }
        }
    }
    
    if(_onlyNeumannBC)    //Check that the matrix is symmetric
    {
        PetscBool isSymetric;
        MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
        MatIsSymmetric(_A,_precision,&isSymetric);
        if(!isSymetric)
        {
            cout<<"Singular matrix is not symmetric, tolerance= "<< _precision<<endl;
            throw CdmathException("Singular matrix should be symmetric with kernel composed of constant vectors");
        }
    }

    _diffusionMatrixSet=true;
    stop=false ;
    
    return INFINITY;
}

double StationaryDiffusionEquation::computeRHS(bool & stop)//Contribution of the PDE RHS to the linear systemm RHS (boundary conditions do contribute to the system RHS via the function computeDiffusionMatrix
{
    if(_mpi_rank == 0)
    {
        if(!_FECalculation)
            for (int i=0; i<_Nmailles;i++)
                VecSetValue(_b,i, _heatTransfertCoeff*_fluidTemperatureField(i) + _heatPowerField(i) ,ADD_VALUES);
        else
        {
            double coeff;// Coefficient to be inserted in RHS
            std::vector< int > nodesId;
            //Contribution de la température fluide au second membre
            if( _heatTransfertCoeff>0 )
            {
                Cell Ci;
                for (int i=0; i<_Nmailles;i++)
                {
                    Ci=_mesh.getCell(i);
                    nodesId=Ci.getNodesId();
                    for (int j=0; j<nodesId.size();j++)
                        if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodesId[j])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[idim])
                        {
                            coeff = _heatTransfertCoeff*_fluidTemperatureField(nodesId[j]);
                            VecSetValue(_b,DiffusionEquation::unknownNodeIndex(nodesId[j], _dirichletNodeIds), coeff*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                        }
                }
            }

            //Contribution du chauffage au second membre
            if( _heatPowerField.getTypeOfField() == NODES )//first order quadrature
            {
                Cell Ci;
                for (int i=0; i<_Nmailles;i++)
                {
                    Ci=_mesh.getCell(i);
                    nodesId=Ci.getNodesId();
                    for (int j=0; j<nodesId.size();j++)
                        if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodesId[j])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[idim])
                        {
                            coeff = _heatPowerField(nodesId[j]);
                            VecSetValue(_b,DiffusionEquation::unknownNodeIndex(nodesId[j], _dirichletNodeIds), coeff*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                        }
                }
            }
            else if( _heatPowerField.getTypeOfField() == FACES )//second order quadrature
            {
                Face Fi;
                Cell Ci1,Ci2;
                int i1, i2;
                double valueFace;
                for (int i=0; i<_mesh.getNumberOfFaces();i++)
                {
                    Fi=_mesh.getFace(i);
                    nodesId = Fi.getNodesId();
                    valueFace = _heatPowerField(i);//_heatTransfertCoeff*_fluidTemperatureField(nodesId[inode]) 

                    i1=Fi.getCellId(0);
                    Ci1=_mesh.getCell(i1);
                    coeff = Ci1.getMeasure();
                    if( ! Fi.isBorder() )
                    {
                        i2=Fi.getCellId(1);
                        Ci2=_mesh.getCell(i2);
                        coeff += Ci2.getMeasure();
                    }
                    coeff *= valueFace/(_Ndim*(_Ndim+1));//Ci1.getNumberOfFaces()=_Ndim+1, Fi.getNumberOfNodes()=_Ndim
                    for(int inode=0; inode<Fi.getNumberOfNodes(); inode++)
                        if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodesId[inode])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[idim]) 
                            VecSetValue(_b,DiffusionEquation::unknownNodeIndex(nodesId[inode], _dirichletNodeIds), coeff,ADD_VALUES);                    
                }
            }
            else
                throw CdmathException("StationaryDiffusionEquation::computeRHS: field heatPowerField should be on NODES or FACES");            
        }
    }
    VecAssemblyBegin(_b);
    VecAssemblyEnd(_b);

    if(_verbose or _system)
        VecView(_b,PETSC_VIEWER_STDOUT_WORLD);

    stop=false ;
    return INFINITY;
}

bool StationaryDiffusionEquation::iterateNewtonStep(bool &converged)
{
    bool stop;

    //Only implicit scheme considered

#if PETSC_VERSION_GREATER_3_5
    KSPSetOperators(_ksp, _A, _A);
#else
    KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

    if(_conditionNumber)
        KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);

    if(_system)
    {
        cout << "Matrice du système linéaire" << endl;
        MatView(_A,PETSC_VIEWER_STDOUT_SELF);
        cout << endl;
        cout << "Second membre du système linéaire" << endl;
        VecView(_b, PETSC_VIEWER_STDOUT_SELF);
        cout << endl;
    }

    KSPSolve(_ksp, _b, _Tk);

    KSPConvergedReason reason;
    KSPGetConvergedReason(_ksp,&reason);
    KSPGetIterationNumber(_ksp, &_PetscIts);
    double residu;
    KSPGetResidualNorm(_ksp,&residu);
    if (reason!=2 and reason!=3)
    {
            PetscPrintf(PETSC_COMM_WORLD,"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc.\n");
            PetscPrintf(PETSC_COMM_WORLD,"!!!!!!!!!!!!! Itérations maximales %d atteintes, résidu = %1.2e, précision demandée= %1.2e.\n",_maxPetscIts,residu,_precision);
            PetscPrintf(PETSC_COMM_WORLD,"Solver used %s, preconditioner %s, Final number of iteration = %d.\n",_ksptype,_pctype,_PetscIts);
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile<<"!!!!!!!!!!!!! Erreur système linéaire : pas de convergence de Petsc."<<endl;
            *_runLogFile<<"!!!!!!!!!!!!! Itérations maximales "<<_maxPetscIts<<" atteintes, résidu="<<residu<<", précision demandée= "<<_precision<<endl;
            *_runLogFile<<"Solver used "<<  _ksptype<<", preconditioner "<<_pctype<<", Final number of iteration= "<<_PetscIts<<endl;
            _runLogFile->close();
        }
        converged = false;
        stop = true;
    }
    else{
        if( _MaxIterLinearSolver < _PetscIts)
            _MaxIterLinearSolver = _PetscIts;
        PetscPrintf(PETSC_COMM_WORLD,"## Système linéaire résolu en %d itérations par le solveur %s et le preconditioneur %s, précision demandée = %1.2e, résidu final  = %1.2e\n",_PetscIts,_ksptype,_pctype,_precision, residu);
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<<"## Système linéaire résolu en "<<_PetscIts<<" itérations par le solveur "<<  _ksptype<<" et le preconditioneur "<<_pctype<<", précision demandée= "<<_precision<<", résidu final = "<< residu<<endl<<endl;

        VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
        VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1

        VecAssemblyBegin(_Tk);
        VecAssemblyEnd(  _Tk);
        VecAssemblyBegin(_deltaT);
        VecAssemblyEnd(  _deltaT);

        if(_verbose)
            PetscPrintf(PETSC_COMM_WORLD,"Début calcul de l'erreur maximale\n");

        VecNorm(_deltaT,NORM_INFINITY,&_erreur_rel);

        if(_verbose)
            PetscPrintf(PETSC_COMM_WORLD,"Fin calcul de la variation relative, erreur maximale : %1.2e\n", _erreur_rel );

        stop=false;
        converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
    }

    VecCopy(_Tk, _Tkm1);

    return stop;
}

void StationaryDiffusionEquation::setMesh(const Mesh &M)
{
    if(_mpi_rank==0)
    {
        if(_Ndim != M.getSpaceDimension() or _Ndim!=M.getMeshDimension())//for the moment we must have space dim=mesh dim
        {
            cout<< "Problem : dimension defined is "<<_Ndim<< " but mesh dimension= "<<M.getMeshDimension()<<", and space dimension is "<<M.getSpaceDimension()<<endl;
            *_runLogFile<< "Problem : dim = "<<_Ndim<< " but mesh dim= "<<M.getMeshDimension()<<", mesh space dim= "<<M.getSpaceDimension()<<endl;
            *_runLogFile<<"StationaryDiffusionEquation::setMesh: mesh has incorrect dimension"<<endl;
            _runLogFile->close();
            throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect  dimension");
        }
    
        _mesh=M;
        _Nmailles = _mesh.getNumberOfCells();
        _Nnodes =   _mesh.getNumberOfNodes();
        
        cout<<"Mesh has "<< _Nmailles << " cells and " << _Nnodes << " nodes"<<endl<<endl;;
        *_runLogFile<<"Mesh has "<< _Nmailles << " cells and " << _Nnodes << " nodes"<<endl<<endl;
        
        // find  maximum nb of neibourghs
        if(!_FECalculation)
        {
            _VV=Field ("Temperature", CELLS, _mesh, 1);
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
                    *_runLogFile<<"StationaryDiffusionEquation::setMesh: mesh has incorrect dimension"<<endl;
                    _runLogFile->close();
                    throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect cell types");
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
                    *_runLogFile<<"StationaryDiffusionEquation::setMesh: mesh has incorrect dimension"<<endl;
                    _runLogFile->close();
                    throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect cell types");
                }
            }
    
            _VV=Field ("Temperature", NODES, _mesh, 1);
    
            _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
            _boundaryNodeIds = _mesh.getBoundaryNodeIds();
            _NboundaryNodes=_boundaryNodeIds.size();
        }
    }
    
    /* MPI distribution parameters */
    int nbVoisinsMax;//Mettre en attribut ?
    if(!_FECalculation){
        MPI_Bcast(&_Nmailles      , 1, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(&_neibMaxNbCells, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        nbVoisinsMax = _neibMaxNbCells;
    }
    else{
        MPI_Bcast(&_Nnodes        , 1, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(&_neibMaxNbNodes, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        nbVoisinsMax = _neibMaxNbNodes;
    }
    _d_nnz = (nbVoisinsMax+1)*_nVar;
    _o_nnz =  nbVoisinsMax   *_nVar;

    _meshSet=true;
}

void StationaryDiffusionEquation::setLinearSolver(linearSolver kspType, preconditioner pcType)
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
        if(_mpi_rank==0)//Avoid redundant printing
        {
            cout << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
            *_runLogFile << "!!! Error : only 'GMRES', 'CG' or 'BCGS' is acceptable as a linear solver !!!" << endl;
            _runLogFile->close();
        }
        throw CdmathException("!!! Error : only 'GMRES', 'CG' or 'BCGS' algorithm is acceptable !!!");
    }
    // set preconditioner
    if (pcType == NOPC)
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
        if(_mpi_rank==0)//Avoid redundant printing
        {
            cout << "!!! Error : only 'NOPC', 'LU', 'ILU', 'CHOLESKY' or 'ICC' preconditioners are acceptable !!!" << endl;
            *_runLogFile << "!!! Error : only 'NOPC' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
            _runLogFile->close();
        }
        throw CdmathException("!!! Error : only 'NOPC' or 'LU' or 'ILU' preconditioners are acceptable !!!" );
    }
}

bool StationaryDiffusionEquation::solveStationaryProblem()
{
    if(!_initializedMemory)
    {
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile<< "ProblemCoreFlows::run() call initialize() first"<< _fileName<<endl;
            _runLogFile->close();
        }
        throw CdmathException("ProblemCoreFlows::run() call initialize() first");
    }
    bool stop=false; // Does the Problem want to stop (error) ?
    bool converged=false; // has the newton scheme converged (end) ?

    PetscPrintf(PETSC_COMM_WORLD,"!!! Running test case %s using ",_fileName.c_str());
    if(_mpi_rank==0)//Avoid redundant printing
        *_runLogFile<< "!!! Running test case "<< _fileName<< " using ";

    if(!_FECalculation)
    {
        PetscPrintf(PETSC_COMM_WORLD,"Finite volumes method\n\n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<< "Finite volumes method"<<endl<<endl;
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD,"Finite elements method\n\n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<< "Finite elements method"<< endl<<endl;
    }

    computeDiffusionMatrix( stop);//For the moment the conductivity does not depend on the temperature (linear LHS)
    if (stop){
        PetscPrintf(PETSC_COMM_WORLD,"Error : failed computing diffusion matrix, stopping calculation\n");
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile << "Error : failed computing diffusion matrix, stopping calculation"<< endl;
             _runLogFile->close();
        }
       throw CdmathException("Failed computing diffusion matrix");
    }
    computeRHS(stop);//For the moment the heat power does not depend on the unknown temperature (linear RHS)
    if (stop){
        PetscPrintf(PETSC_COMM_WORLD,"Error : failed computing right hand side, stopping calculation\n");
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile << "Error : failed computing right hand side, stopping calculation"<< endl;
            _runLogFile->close();
        }
        throw CdmathException("Failed computing right hand side");
    }
    stop = iterateNewtonStep(converged);
    if (stop){
        PetscPrintf(PETSC_COMM_WORLD,"Error : failed solving linear system, stopping calculation\n");
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile << "Error : failed linear system, stopping calculation"<< endl;
            _runLogFile->close();
        }
        throw CdmathException("Failed solving linear system");
    }
    
    _computationCompletedSuccessfully=true;
    save();

    // Newton iteration loop for non linear problems
    /*
    while(!stop and !converged and _NEWTON_its<_maxNewtonIts)
    {
        computeDiffusionMatrix( stop);//case when the conductivity depends on the temperature (nonlinear LHS)
        computeRHS(stop);//case the heat power depends on the unknown temperature (nonlinear RHS)
        stop = iterateNewtonStep(converged);
        _NEWTON_its++;
    }
    if (stop){
        cout << "Error : failed solving Newton iteration "<<_NEWTON_its<<", stopping calculation"<< endl;
        *_runLogFile << "Error : failed solving Newton iteration "<<_NEWTON_its<<", stopping calculation"<< endl;
        throw CdmathException("Failed solving a Newton iteration");
    }
    else if(_NEWTON_its==_maxNewtonIts){
        cout << "Error : no convergence of Newton scheme. Maximum Newton iterations "<<_maxNewtonIts<<" reached, stopping calculation"<< endl;
        *_runLogFile << "Error : no convergence of Newton scheme. Maximum Newton iterations "<<_maxNewtonIts<<" reached, stopping calculation"<< endl;
        throw CdmathException("No convergence of Newton scheme");
    }
    else{
        cout << "Convergence of Newton scheme at iteration "<<_NEWTON_its<<", end of calculation"<< endl;
        *_runLogFile << "Convergence of Newton scheme at iteration "<<_NEWTON_its<<", end of calculation"<< endl;
        save();
    }
    */
    
    if(_mpi_rank==0)//Avoid redundant printing
    {
        *_runLogFile<< "!!!!!! Computation successful !!!!!!"<< endl;
        _runLogFile->close();
    }

    return !stop;
}

void StationaryDiffusionEquation::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results\n\n");
    if(_mpi_rank==0)//Avoid redundant printing
        *_runLogFile<< "Saving numerical results"<< endl<<endl;

    string resultFile(_path+"/StationaryDiffusionEquation");//Results

    resultFile+="_";
    resultFile+=_fileName;

    // create mesh and component info
    string suppress ="rm -rf "+resultFile+"_*";
    system(suppress.c_str());//Nettoyage des précédents calculs identiques
    
    if(_mpi_size>1){
        VecScatterBegin(_scat,_Tk,_Tk_seq,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(  _scat,_Tk,_Tk_seq,INSERT_VALUES,SCATTER_FORWARD);
    }

    if(_verbose or _system)
        VecView(_Tk,PETSC_VIEWER_STDOUT_WORLD);

    if(_mpi_rank==0){
        double Ti; 
        if(!_FECalculation)
            for(int i=0; i<_Nmailles; i++)
                {
                    if(_mpi_size>1)
                        VecGetValues(_Tk_seq, 1, &i, &Ti);
                    else
                        VecGetValues(_Tk    , 1, &i, &Ti);
                    _VV(i)=Ti;
                }
        else
        {
        /* Noeuds inconnus */
            int globalIndex;
            for(int i=0; i<_NunknownNodes; i++)
            {
                if(_mpi_size>1)
                    VecGetValues(_Tk_seq, 1, &i, &Ti);
                else
                    VecGetValues(_Tk    , 1, &i, &Ti);
                globalIndex = DiffusionEquation::globalNodeIndex(i, _dirichletNodeIds);
                _VV(globalIndex)=Ti;
            }
    
            /* Noeuds connus (condition limite de Dirichlet */
            Node Ni;
            string nameOfGroup;
            std::map<int,double>::iterator it;
            for(int i=0; i<_NdirichletNodes; i++)
            {
                it=_dirichletBoundaryValues.find(_dirichletNodeIds[i]);
                if( it != _dirichletBoundaryValues.end() )//Une valeur limite est associée au noeud
                    _VV(_dirichletNodeIds[i])=it->second;
                else//Une valeur limite est associée au groupe frontière    
                {
                    Ni=_mesh.getNode(_dirichletNodeIds[i]);
                    nameOfGroup = Ni.getGroupName();
                    _VV(_dirichletNodeIds[i])=_limitField[nameOfGroup].T;
                }
            }
        }
    
        _VV.setInfoOnComponent(0,"Temperature_(K)");
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
}

void StationaryDiffusionEquation::terminate()
{
    VecDestroy(&_Tk);
    VecDestroy(&_Tkm1);
    VecDestroy(&_deltaT);
    VecDestroy(&_b);
    MatDestroy(&_A);
    if(_mpi_size>1 && _mpi_rank == 0)
        VecDestroy(&_Tk_seq);

    PetscBool petscInitialized;
    PetscInitialized(&petscInitialized);
    if(petscInitialized)
        PetscFinalize();

    delete _runLogFile;
}
void 
StationaryDiffusionEquation::setDirichletValues(map< int, double> dirichletBoundaryValues)
{
    _dirichletBoundaryValues=dirichletBoundaryValues;
}

void 
StationaryDiffusionEquation::setNeumannValues(map< int, double> neumannBoundaryValues)
{
    _neumannBoundaryValues=neumannBoundaryValues;
}

double 
StationaryDiffusionEquation::getConditionNumber(bool isSingular, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getConditionNumber( isSingular, tol);
}
std::vector< double > 
StationaryDiffusionEquation::getEigenvalues(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  
  if(_FECalculation)//We need to scale the FE matrix, otherwise the eigenvalues go to zero as the mesh is refined
  {
      Vector nodal_volumes(_NunknownNodes);
      int j_int;
      for(int i = 0; i< _Nmailles ; i++)//On parcourt les cellules du maillage
      {
        Cell Ci = _mesh.getCell(i);
        for(int j = 0 ; j<_Ndim+1 ; j++)//On parcourt les noeuds de la cellule
        {
            if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),Ci.getNodeId(j))==_dirichletNodeIds.end())//node j is an unknown node (not a Dirichlet node)
            {
                j_int=DiffusionEquation::unknownNodeIndex(Ci.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
                nodal_volumes[j_int]+=Ci.getMeasure()/(_Ndim+1);
            }
        }
       }
      for( j_int = 0; j_int< _NunknownNodes ; j_int++)
        nodal_volumes[j_int]=1/nodal_volumes[j_int];
      A.leftDiagonalScale(nodal_volumes);
  }

  return A.getEigenvalues( nev, which, tol);
}
std::vector< Vector > 
StationaryDiffusionEquation::getEigenvectors(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  return A.getEigenvectors( nev, which, tol);
}
Field 
StationaryDiffusionEquation::getEigenvectorsField(int nev, EPSWhich which, double tol) const
{
  SparseMatrixPetsc A = SparseMatrixPetsc(_A);
  MEDCoupling::DataArrayDouble * d = A.getEigenvectorsDataArrayDouble( nev, which, tol);
  Field my_eigenfield;
  
    my_eigenfield = Field("Eigenvectors field", _VV.getTypeOfField(), _mesh, nev);

  my_eigenfield.setFieldByDataArrayDouble(d);
  
  return my_eigenfield;
}

Field&
StationaryDiffusionEquation::getOutputTemperatureField()
{
    if(!_computationCompletedSuccessfully)
        throw("Computation not performed yet or failed. No temperature field available");
    else
        return _VV;
}

Field& 
StationaryDiffusionEquation::getRodTemperatureField()
{
   return getOutputTemperatureField();
}

vector<string> 
StationaryDiffusionEquation::getInputFieldsNames()
{
    vector<string> result(2);
    
    result[0]="FluidTemperature";
    result[1]="HeatPower";
    
    return result;
}
vector<string> 
StationaryDiffusionEquation::getOutputFieldsNames()
{
    vector<string> result(1);
    
    result[0]="RodTemperature";
    
    return result;
}

Field& 
StationaryDiffusionEquation::getOutputField(const string& nameField )
{
    if(nameField=="RodTemperature" || nameField=="RODTEMPERATURE" || nameField=="TEMPERATURECOMBUSTIBLE" || nameField=="TemperatureCombustible" )
        return getRodTemperatureField();
    else
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Error : Field name %s is not an input field name, call getOutputFieldsNames first\n", nameField.c_str());
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile<< "Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first"<< endl;
            _runLogFile->close();
        }
        throw CdmathException("DiffusionEquation::getOutputField error : Unknown Field name");
    }
}

void
StationaryDiffusionEquation::setInputField(const string& nameField, Field& inputField )
{
    if(!_meshSet)
        throw CdmathException("!!!!!!!! StationaryDiffusionEquation::setInputField set the mesh first");

    if(nameField=="FluidTemperature" || nameField=="FLUIDTEMPERATURE" || nameField=="TemperatureFluide" || nameField=="TEMPERATUREFLUIDE")
        return setFluidTemperatureField( inputField) ;
    else if(nameField=="HeatPower" || nameField=="HEATPOWER" || nameField=="PuissanceThermique" || nameField=="PUISSANCETHERMIQUE" )
        return setHeatPowerField( inputField );
    else
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Error : Field name %s is not an input field name, call getInputFieldsNames first\n", nameField.c_str());
        if(_mpi_rank==0)//Avoid redundant printing
        {
            *_runLogFile<< "Error : Field name "<< nameField << " is not an input field name, call getInputFieldsNames first"<< endl;
            _runLogFile->close();
        }
        throw CdmathException("StationaryDiffusionEquation::setInputField error : Unknown Field name");
    }
}

void 
StationaryDiffusionEquation::setFluidTemperatureField(Field coupledTemperatureField){
    if(!_meshSet)
        throw CdmathException("!!!!!!!! StationaryDiffusionEquation::setFluidTemperatureField set initial field first");

    coupledTemperatureField.getMesh().checkFastEquivalWith(_mesh);
    _fluidTemperatureField=coupledTemperatureField;
    _fluidTemperatureFieldSet=true;
};

void 
StationaryDiffusionEquation::setHeatPowerField(Field heatPower){
    if(!_meshSet)
        throw CdmathException("!!!!!!!! StationaryDiffusionEquation::setHeatPowerField set initial field first");

    heatPower.getMesh().checkFastEquivalWith(_mesh);
    _heatPowerField=heatPower;
    _heatPowerFieldSet=true;
}

void 
StationaryDiffusionEquation::setHeatPowerField(string fileName, string fieldName, EntityType field_support_type, int iteration, int order, int meshLevel){
    if(!_meshSet)
        throw CdmathException("!!!!!!!! StationaryDiffusionEquation::setHeatPowerField set initial field first");

    _heatPowerField=Field(fileName, field_support_type,fieldName, iteration, order, meshLevel);
    _heatPowerField.getMesh().checkFastEquivalWith(_mesh);
    _heatPowerFieldSet=true;
}

void 
StationaryDiffusionEquation::setDirichletBoundaryCondition(string groupName, string fileName, string fieldName, int timeStepNumber, int order, int meshLevel, EntityType field_support_type){
    if(_FECalculation && field_support_type != NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<< "Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }
    else if(!_FECalculation && field_support_type == NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<<"Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }

    setDirichletBoundaryCondition( groupName, Field(fileName, field_support_type, fieldName, timeStepNumber, order, meshLevel));
}

void StationaryDiffusionEquation::setDirichletBoundaryCondition(string groupName, Field bc_field){
    if(_FECalculation && bc_field.getTypeOfField() != NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<< "Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }
    else if(!_FECalculation && bc_field.getTypeOfField() == NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<<"Warning : StationaryDiffusionEquation::setDirichletBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }
    
    _dirichletBoundaryField = bc_field;
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> dirichletBoundaryMesh = bc_field.getMesh().getMEDCouplingMesh();

    //* Check that the boundary field is based on the correct boundary mesh */
    int compType=2;//This is the weakest comparison policy for medcoupling meshes. It can be used by users not sensitive to cell orientation
    MEDCoupling::DataArrayIdType * arr;//DataArrayIdType to contain the correspondence between cells of the two meshes
    MEDCoupling::MEDCouplingUMesh* dirichletBoundaryUMesh = dynamic_cast<MEDCoupling::MEDCouplingUMesh*> ( dirichletBoundaryMesh.retn());
    
    if( !_mesh.getBoundaryMEDCouplingMesh()->areCellsIncludedIn(dirichletBoundaryUMesh, compType, arr) )
        throw CdmathException(" !!!!! StationaryDiffusionEquation::setDirichletBoundaryCondition : The boundary field is not based on the correct boundary mesh. Use mesh::getBoundaryMesh");

    int nBoundaryCells = dirichletBoundaryUMesh->getNumberOfCells();//Boundary cells that support the boundary field (subpart of the total boundary)
    std::map<int,double>::iterator it;
    long int iCell_global;
    if(!_FECalculation)//Finite volume simulation
        for(int i=0; i<nBoundaryCells ; i++)
        {
            arr->getTuple(i,&iCell_global);
            it=_dirichletBoundaryValues.find(iCell_global);
            if( it == _dirichletBoundaryValues.end() )//Aucune valeur limite n'est associée au noeud
                it->second = bc_field[i];
        }
    long int inode_global, length;
    if(_FECalculation)
    {
        const MEDCoupling::DataArrayIdType *globalNodeId  = dirichletBoundaryUMesh->computeFetchedNodeIds();
        for(int i=0; i<bc_field.getNumberOfElements(); i++)
        {
            globalNodeId->getTuple(i,&inode_global);
            it=_dirichletBoundaryValues.find(inode_global);
            if( it == _dirichletBoundaryValues.end() )//Aucune valeur limite n'est associée au noeud
                it->second = bc_field[i];  
        }
        globalNodeId->decrRef();
    }
    arr->decrRef();
};

void 
StationaryDiffusionEquation::setNeumannBoundaryCondition(string groupName, string fileName, string fieldName, int timeStepNumber, int order, int meshLevel, EntityType field_support_type){
    if(_FECalculation && field_support_type != NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setNeumannBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<<"Warning : StationaryDiffusionEquation::setNeumannBoundaryCondition : finite element simulation should have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }
    else if(!_FECalculation && field_support_type == NODES)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n Warning : StationaryDiffusionEquation::setNeumannBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type \n");
        if(_mpi_rank==0)//Avoid redundant printing
            *_runLogFile<<"Warning : StationaryDiffusionEquation::setNeumannBoundaryCondition : finite volume simulation should not have boundary field on nodes!!! Change parameter field_support_type"<< endl;
    }

    setNeumannBoundaryCondition( groupName, Field(fileName, field_support_type, fieldName, timeStepNumber, order, meshLevel) );    
}

void StationaryDiffusionEquation::setNeumannBoundaryCondition(string groupName, Field bc_field){
    _neumannBoundaryField = bc_field;
       MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> neumannBoundaryMesh = _neumannBoundaryField.getMesh().getMEDCouplingMesh();

    //* Check that the boundary field is based on the correct boundary mesh */
    int compType=2;//This is the weakest comparison policy for medcoupling meshes. It can be used by users not sensitive to cell orientation
    MEDCoupling::DataArrayIdType * arr;//DataArrayIdType to contain the correspondence between cells of the two meshes
    MEDCoupling::MEDCouplingUMesh* neumannBoundaryUMesh = dynamic_cast<MEDCoupling::MEDCouplingUMesh*> ( neumannBoundaryMesh.retn());

    if( !_mesh.getBoundaryMEDCouplingMesh()->areCellsIncludedIn(neumannBoundaryUMesh, compType, arr) )
        throw CdmathException(" !!!!! StationaryDiffusionEquation::setNeumannBoundaryCondition : The boundary field is not based on the correct boundary mesh. Use mesh::getBoundaryMesh");
};
