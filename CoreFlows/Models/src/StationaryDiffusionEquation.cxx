#include "StationaryDiffusionEquation.hxx"
#include "Node.hxx"
#include "SparseMatrixPetsc.hxx"
#include "math.h"
#include <algorithm> 
#include <fstream>
#include <sstream>

using namespace std;

int StationaryDiffusionEquation::fact(int n)
{
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}
int StationaryDiffusionEquation::unknownNodeIndex(int globalIndex, std::vector< int > dirichletNodes) const
{//assumes Dirichlet node numbering is strictly increasing
    int j=0;//indice de parcours des noeuds frontière avec CL Dirichlet
    int boundarySize=dirichletNodes.size();
    while(j<boundarySize and dirichletNodes[j]<globalIndex)
        j++;
    if(j==boundarySize)
        return globalIndex-boundarySize;
    else if (dirichletNodes[j]>globalIndex)
        return globalIndex-j;
    else
        throw CdmathException("StationaryDiffusionEquation::unknownNodeIndex : Error : node is a Dirichlet boundary node");
}

int StationaryDiffusionEquation::globalNodeIndex(int unknownNodeIndex, std::vector< int > dirichletNodes) const
{//assumes Dirichlet boundary node numbering is strictly increasing
    int boundarySize=dirichletNodes.size();
    /* trivial case where all boundary nodes are Neumann BC */
    if(boundarySize==0)
        return unknownNodeIndex;
        
    double unknownNodeMax=-1;//max unknown node number in the interval between jth and (j+1)th Dirichlet boundary nodes
    int j=0;//indice de parcours des noeuds frontière
    //On cherche l'intervale [j,j+1] qui contient le noeud de numéro interieur unknownNodeIndex
    while(j+1<boundarySize and unknownNodeMax<unknownNodeIndex)
    {
        unknownNodeMax += dirichletNodes[j+1]-dirichletNodes[j]-1;
        j++;
    }    

    if(j+1==boundarySize)
        return unknownNodeIndex+boundarySize;
    else //unknownNodeMax>=unknownNodeIndex, hence our node global number is between dirichletNodes[j-1] and dirichletNodes[j]
        return unknownNodeIndex - unknownNodeMax + dirichletNodes[j]-1;
}

StationaryDiffusionEquation::StationaryDiffusionEquation(int dim, bool FECalculation, double lambda){
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
		PetscInitialize(NULL,NULL,0,0);

    if(lambda < 0.)
    {
        std::cout<<"Conductivity="<<lambda<<endl;
        throw CdmathException("Error : conductivity parameter lambda cannot  be negative");
    }
    if(dim<=0)
    {
        std::cout<<"Space dimension="<<dim<<endl;
        throw CdmathException("Error : parameter dim cannot  be negative");
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
    _dirichletValuesSet=false;   
    _neumannValuesSet=false;   
    
    //Linear solver data
	_precision=1.e-6;
	_precision_Newton=_precision;
	_MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
	_maxPetscIts=50;
	_maxNewtonIts=50;
	_NEWTON_its=0;
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
}

void StationaryDiffusionEquation::initialize()
{
	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation

	if(!_meshSet)
		throw CdmathException("StationaryDiffusionEquation::initialize() set mesh first");
	else
    {
		cout<<"!!!! Initialisation of the computation of the temperature diffusion in a solid using ";
        *_runLogFile<<"!!!!! Initialisation of the computation of the temperature diffusion in a solid using ";
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
    
	_DiffusionTensor=Matrix(_Ndim);
	for(int idim=0;idim<_Ndim;idim++)
		_DiffusionTensor(idim,idim)=1;
	/**************** Field creation *********************/

	if(!_heatPowerFieldSet){
        if(_FECalculation){
            _heatPowerField=Field("Heat power",NODES,_mesh,1);
            for(int i =0; i<_Nnodes; i++)
                _heatPowerField(i) = _heatSource;
        }
        else{
            _heatPowerField=Field("Heat power",CELLS,_mesh,1);
            for(int i =0; i<_Nmailles; i++)
                _heatPowerField(i) = _heatSource;
        }
        _heatPowerFieldSet=true;
    }
	if(!_fluidTemperatureFieldSet){
        if(_FECalculation){
            _fluidTemperatureField=Field("Fluid temperature",NODES,_mesh,1);
            for(int i =0; i<_Nnodes; i++)
                _fluidTemperatureField(i) = _fluidTemperature;
        }
        else{
            _fluidTemperatureField=Field("Fluid temperature",CELLS,_mesh,1);
            for(int i =0; i<_Nmailles; i++)
                _fluidTemperatureField(i) = _fluidTemperature;
        }
        _fluidTemperatureFieldSet=true;
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

	//creation de la matrice
    if(!_FECalculation)
        MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles, _Nmailles, (1+_neibMaxNbCells), PETSC_NULL, &_A);
    else
        MatCreateSeqAIJ(PETSC_COMM_SELF, _NunknownNodes, _NunknownNodes, (1+_neibMaxNbNodes), PETSC_NULL, &_A);

	VecCreate(PETSC_COMM_SELF, &_Tk);

    if(!_FECalculation)
        VecSetSizes(_Tk,PETSC_DECIDE,_Nmailles);
    else
        VecSetSizes(_Tk,PETSC_DECIDE,_NunknownNodes);

	VecSetFromOptions(_Tk);
	VecDuplicate(_Tk, &_Tkm1);
	VecDuplicate(_Tk, &_deltaT);
	VecDuplicate(_Tk, &_b);//RHS of the linear system

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

    //Checking whether all boundary conditions are Neumann boundary condition
    //if(_FECalculation) _onlyNeumannBC = _NdirichletNodes==0;
    if(!_neumannValuesSet)//Boundary conditions set via LimitField structure
    {
        map<string, LimitFieldStationaryDiffusion>::iterator it = _limitField.begin();
        while(it != _limitField.end() and (it->second).bcType == NeumannStationaryDiffusion)
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
        std::cout<<"### Warning : all boundary conditions are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
        std::cout<<"### Check the compatibility condition between the right hand side and the boundary data. For homogeneous Neumann BCs, the right hand side must have integral equal to zero."<<std::endl;
        std::cout<<"### The system matrix being singular, we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;
        *_runLogFile<<"### Warning : all boundary condition are Neumann. System matrix is not invertible since constant vectors are in the kernel."<<std::endl;
        *_runLogFile<<"### The system matrix being singular, we seek a zero sum solution, and exact (LU and CHOLESKY) and incomplete factorisations (ILU and ICC) may fail."<<std::endl<<endl;
        *_runLogFile<<"### Check the compatibility condition between the right hand side and the boundary data. For homogeneous Neumann BCs, the right hand side must have integral equal to zero."<<std::endl;

		//Check that the matrix is symmetric
		PetscBool isSymetric;
		MatIsSymmetric(_A,_precision,&isSymetric);
		if(!isSymetric)
			{
				cout<<"Singular matrix is not symmetric, tolerance= "<< _precision<<endl;
				throw CdmathException("Singular matrix should be symmetric with kernel composed of constant vectors");
			}
		MatNullSpace nullsp;
		MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
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

Vector StationaryDiffusionEquation::gradientNodal(Matrix M, vector< double > values){
    vector< Matrix > matrices(_Ndim);
    
    for (int idim=0; idim<_Ndim;idim++){
        matrices[idim]=M.deepCopy();
        for (int jdim=0; jdim<_Ndim+1;jdim++)
			matrices[idim](jdim,idim) = values[jdim] ;
    }

	Vector result(_Ndim);
    for (int idim=0; idim<_Ndim;idim++)
        result[idim] = matrices[idim].determinant();

	return result;    
}

double StationaryDiffusionEquation::computeDiffusionMatrix(bool & stop)
{
    double result;
    
    if(_FECalculation)
        result=computeDiffusionMatrixFE(stop);
    else
        result=computeDiffusionMatrixFV(stop);

    //Contribution from the solid/fluid heat exchange with assumption of constant heat transfer coefficient
    //update value here if variable  heat transfer coefficient
    if(_heatTransfertCoeff>_precision)
        MatShift(_A,_heatTransfertCoeff);//Contribution from the liquit/solid heat transfer
        
    if(_verbose or _system)
        MatView(_A,PETSC_VIEWER_STDOUT_SELF);

    return  result;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFE(bool & stop){
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
    vector< double > valuesBorder(_Ndim+1);
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
            GradShapeFuncs[idim]=gradientNodal(M,values[idim])/fact(_Ndim);

        /* Loop on the edges of the cell */
        for (int idim=0; idim<_Ndim+1;idim++)
        {
            if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[idim])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[idim])
            {//First node of the edge is not Dirichlet node
                i_int=unknownNodeIndex(nodeIds[idim], _dirichletNodeIds);//assumes Dirichlet boundary node numbering is strictly increasing
                dirichletCell_treated=false;
                for (int jdim=0; jdim<_Ndim+1;jdim++)
                {
                    if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[jdim])==_dirichletNodeIds.end())//!_mesh.isBorderNode(nodeIds[jdim])
                    {//Second node of the edge is not Dirichlet node
                        j_int= unknownNodeIndex(nodeIds[jdim], _dirichletNodeIds);//assumes Dirichlet boundary node numbering is strictly increasing
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
                                if( _dirichletValuesSet )
                                    valuesBorder[kdim]=_dirichletBoundaryValues[it->second];
                                else    
                                    valuesBorder[kdim]=_limitField[_mesh.getNode(nodeIds[kdim]).getGroupName()].T;
                            }
                            else
                                valuesBorder[kdim]=0;                            
                        }
                        GradShapeFuncBorder=gradientNodal(M,valuesBorder)/fact(_Ndim);
                        coeff =-_conductivity*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
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
            Face Fi = _mesh.getFace(i);
            for(int j = 0 ; j<_Ndim ; j++)//On parcourt les noeuds de la face
            {
                if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),Fi.getNodeId(j))==_dirichletNodeIds.end())//node j is an Neumann BC node (not a Dirichlet BC node)
                {
                    j_int=unknownNodeIndex(Fi.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
                    if( _neumannValuesSet )
                        coeff =Fi.getMeasure()/_Ndim*_neumannBoundaryValues[Fi.getNodeId(j)];
                    else    
                        coeff =Fi.getMeasure()/_Ndim*_limitField[_mesh.getNode(Fi.getNodeId(j)).getGroupName()].normalFlux;
                    VecSetValue(_b, j_int, coeff, ADD_VALUES);
                }
            }
        }
    }
    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);

	_diffusionMatrixSet=true;
    stop=false ;

	return INFINITY;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFV(bool & stop){
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
                    cout<<"!!!!!! Boundary condition not accepted for boundary named !!!!!!!!!!"<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
                    cout<<"Accepted boundary conditions are NeumannStationaryDiffusion "<<NeumannStationaryDiffusion<< " and DirichletStationaryDiffusion "<<DirichletStationaryDiffusion<<endl;
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
            *_runLogFile<<"StationaryDiffusionEquation::computeDiffusionMatrixFV(): incompatible number of cells around a face"<<endl;
			throw CdmathException("StationaryDiffusionEquation::computeDiffusionMatrixFV(): incompatible number of cells around a face");
        }
	}

	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);
    
	_diffusionMatrixSet=true;
    stop=false ;
	
    return INFINITY;
}

double StationaryDiffusionEquation::computeRHS(bool & stop)//Contribution of the PDE RHS to the linear systemm RHS (boundary conditions do contribute to the system RHS via the function computeDiffusionMatrix)
{
	VecAssemblyBegin(_b);

    if(!_FECalculation)
        for (int i=0; i<_Nmailles;i++)
            VecSetValue(_b,i,_heatTransfertCoeff*_fluidTemperatureField(i) + _heatPowerField(i),ADD_VALUES);
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
                    {
                        double coeff = _heatTransfertCoeff*_fluidTemperatureField(nodesId[j]) + _heatPowerField(nodesId[j]);
                        VecSetValue(_b,unknownNodeIndex(nodesId[j], _dirichletNodeIds), coeff*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                    }
            }
        }
    
	VecAssemblyEnd(_b);

    if(_verbose or _system)
        VecView(_b,PETSC_VIEWER_STDOUT_SELF);

    stop=false ;
	return INFINITY;
}

bool StationaryDiffusionEquation::iterateNewtonStep(bool &converged)
{
	bool stop;

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
    KSPSolve(_ksp, _b, _Tk);

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
        converged = false;
        stop = true;
    }
    else{
        if( _MaxIterLinearSolver < _PetscIts)
            _MaxIterLinearSolver = _PetscIts;
        cout<<"## Système linéaire résolu en "<<_PetscIts<<" itérations par le solveur "<<  _ksptype<<" et le preconditioneur "<<_pctype<<", précision demandée= "<<_precision<<endl<<endl;
		*_runLogFile<<"## Système linéaire résolu en "<<_PetscIts<<" itérations par le solveur "<<  _ksptype<<" et le preconditioneur "<<_pctype<<", précision demandée= "<<_precision<<endl<<endl;
        VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
        VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1

        _erreur_rel= 0;
        double Ti, dTi;

        VecAssemblyBegin(_Tk);
        VecAssemblyEnd(  _Tk);
        VecAssemblyBegin(_deltaT);
        VecAssemblyEnd(  _deltaT);

        if(!_FECalculation)
            for(int i=0; i<_Nmailles; i++)
            {
                VecGetValues(_deltaT, 1, &i, &dTi);
                VecGetValues(_Tk, 1, &i, &Ti);
                if(_erreur_rel < fabs(dTi/Ti))
                    _erreur_rel = fabs(dTi/Ti);
            }
        else
            for(int i=0; i<_NunknownNodes; i++)
            {
                VecGetValues(_deltaT, 1, &i, &dTi);
                VecGetValues(_Tk, 1, &i, &Ti);
                if(_erreur_rel < fabs(dTi/Ti))
                    _erreur_rel = fabs(dTi/Ti);
            }
        stop=false;
        converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
    }

	VecCopy(_Tk, _Tkm1);

	return stop;
}

void StationaryDiffusionEquation::setMesh(const Mesh &M)
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

bool StationaryDiffusionEquation::solveStationaryProblem()
{
	if(!_initializedMemory)
	{
		*_runLogFile<< "ProblemCoreFlows::run() call initialize() first"<< _fileName<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::run() call initialize() first");
	}
	bool stop=false; // Does the Problem want to stop (error) ?
	bool converged=false; // has the newton scheme converged (end) ?

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

    computeDiffusionMatrix( stop);//For the moment the conductivity does not depend on the temperature (linear LHS)
    if (stop){
        cout << "Error : failed computing diffusion matrix, stopping calculation"<< endl;
        *_runLogFile << "Error : failed computing diffusion matrix, stopping calculation"<< endl;
 		_runLogFile->close();
       throw CdmathException("Failed computing diffusion matrix");
    }
    computeRHS(stop);//For the moment the heat power does not depend on the unknown temperature (linear RHS)
    if (stop){
        cout << "Error : failed computing right hand side, stopping calculation"<< endl;
        *_runLogFile << "Error : failed computing right hand side, stopping calculation"<< endl;
        throw CdmathException("Failed computing right hand side");
    }
    stop = iterateNewtonStep(converged);
    if (stop){
        cout << "Error : failed solving linear system, stopping calculation"<< endl;
        *_runLogFile << "Error : failed linear system, stopping calculation"<< endl;
		_runLogFile->close();
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
    
    *_runLogFile<< "!!!!!! Computation successful !!!!!!"<< endl;
	_runLogFile->close();

	return !stop;
}

void StationaryDiffusionEquation::save(){
    cout<< "Saving numerical results"<<endl<<endl;
    *_runLogFile<< "Saving numerical results"<< endl<<endl;

	string resultFile(_path+"/StationaryDiffusionEquation");//Results

	resultFile+="_";
	resultFile+=_fileName;

	// create mesh and component info
    string suppress ="rm -rf "+resultFile+"_*";
    system(suppress.c_str());//Nettoyage des précédents calculs identiques
    
    if(_verbose or _system)
        VecView(_Tk,PETSC_VIEWER_STDOUT_SELF);

    double Ti; 
    if(!_FECalculation)
        for(int i=0; i<_Nmailles; i++)
            {
                VecGetValues(_Tk, 1, &i, &Ti);
                _VV(i)=Ti;
            }
    else
    {
        int globalIndex;
        for(int i=0; i<_NunknownNodes; i++)
        {
            VecGetValues(_Tk, 1, &i, &Ti);
            globalIndex = globalNodeIndex(i, _dirichletNodeIds);
            _VV(globalIndex)=Ti;
        }

        Node Ni;
        string nameOfGroup;
        for(int i=0; i<_NdirichletNodes; i++)
        {
            Ni=_mesh.getNode(_dirichletNodeIds[i]);
            nameOfGroup = Ni.getGroupName();
            _VV(_dirichletNodeIds[i])=_limitField[nameOfGroup].T;
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
Field 
StationaryDiffusionEquation::getOutputTemperatureField()
{
    if(!_computationCompletedSuccessfully)
        throw("Computation not performed yet or failed. No temperature field available");
    else
        return _VV;
}

void StationaryDiffusionEquation::terminate()
{
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b);
	MatDestroy(&_A);
}
void 
StationaryDiffusionEquation::setDirichletValues(map< int, double> dirichletBoundaryValues)
{
    _dirichletValuesSet=true;
    _dirichletBoundaryValues=dirichletBoundaryValues;
}

void 
StationaryDiffusionEquation::setNeumannValues(map< int, double> neumannBoundaryValues)
{
    _neumannValuesSet=true;
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
                j_int=unknownNodeIndex(Ci.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
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
  
  if(_FECalculation)
    my_eigenfield = Field("Eigenvectors field", NODES, _mesh, nev);
  else
    my_eigenfield = Field("Eigenvectors field", CELLS, _mesh, nev);

  my_eigenfield.setFieldByDataArrayDouble(d);
  
  return my_eigenfield;
}
