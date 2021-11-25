#include "DiffusionEquation.hxx"
#include "math.h"
#include <algorithm> 
#include <fstream>
#include <sstream>

using namespace std;

int DiffusionEquation::fact(int n)
{
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}
int DiffusionEquation::unknownNodeIndex(int globalIndex, std::vector< int > dirichletNodes)
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
        throw CdmathException("DiffusionEquation::unknownNodeIndex : Error : node is a Dirichlet boundary node");
}

int DiffusionEquation::globalNodeIndex(int unknownNodeIndex, std::vector< int > dirichletNodes)
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

Vector DiffusionEquation::gradientNodal(Matrix M, vector< double > values)
{
	if(! M.isSquare() )
		throw CdmathException("DiffusionEquation::gradientNodal Matrix M should be square !!!");
		
	int Ndim = M.getNumberOfRows()-1;
    vector< Matrix > matrices(Ndim);
    
    for (int idim=0; idim<Ndim;idim++){
        matrices[idim]=M.deepCopy();
        for (int jdim=0; jdim<Ndim+1;jdim++)
			matrices[idim](jdim,idim) = values[jdim] ;
    }

	Vector result(Ndim);
    for (int idim=0; idim<Ndim;idim++)
        result[idim] = matrices[idim].determinant();

	return result;    
}

DiffusionEquation::DiffusionEquation(int dim, bool FECalculation,double rho,double cp, double lambda, MPI_Comm comm ):ProblemCoreFlows(comm)
{
    /* Control input value are acceptable */
    if(rho<_precision or cp<_precision)
    {
        PetscPrintf(PETSC_COMM_WORLD,"rho = %.2f, cp = %.2f, precision = %.2e\n",rho,cp,_precision);
        throw CdmathException("Error : parameters rho and cp should be strictly positive");
    }
    if(lambda < 0.)
    {
        PetscPrintf(PETSC_COMM_WORLD,"Conductivity = %.2f\n",lambda);
        throw CdmathException("Error : conductivity parameter lambda cannot  be negative");
    }
    if(dim<=0)
    {
        PetscPrintf(PETSC_COMM_WORLD,"Space dimension = %.2f\n",dim);
        throw CdmathException("Error : parameter dim cannot  be negative");
    }

    PetscPrintf(PETSC_COMM_WORLD,"\n Diffusion problem with density %.2e, specific heat %.2e, conductivity %.2e", rho,cp,lambda);
    if(FECalculation)
        PetscPrintf(PETSC_COMM_WORLD," and finite elements method\n\n");
    else
        PetscPrintf(PETSC_COMM_WORLD," and finite volumes method\n\n");
    
    _FECalculation=FECalculation;
    
    /* Finite element data */
    _boundaryNodeIds=std::vector< int >(0);
    _dirichletNodeIds=std::vector< int >(0);
    _NboundaryNodes=0;
    _NdirichletNodes=0;
    _NunknownNodes=0;
    _dirichletValuesSet=false;   
    _neumannValuesSet=false;   

    /* Physical parameters */
	_conductivity=lambda;
	_cp=cp;
	_rho=rho;
	_diffusivity=_conductivity/(_rho*_cp);
	_fluidTemperatureFieldSet=false;
	_fluidTemperature=0;
    
    /* Numerical parameters */
	_Ndim=dim;
	_nVar=1;
	_dt_diffusion=0;
	_dt_src=0;
	_diffusionMatrixSet=false;

	_fileName = "SolverlabDiffusionProblem";

    /* Default diffusion tensor is diagonal */
   	_DiffusionTensor=Matrix(_Ndim);
	for(int idim=0;idim<_Ndim;idim++)
		_DiffusionTensor(idim,idim)=_diffusivity;
}

void DiffusionEquation::initialize()
{
	if(_mpi_rank==0)
	{
		_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation
	
		if(_Ndim != _mesh.getSpaceDimension() or _Ndim!=_mesh.getMeshDimension())//for the moment we must have space dim=mesh dim
		{
	        PetscPrintf(PETSC_COMM_SELF,"Problem : dimension defined is %d but mesh dimension= %d, and space dimension is %d",_Ndim,_mesh.getMeshDimension(),_mesh.getSpaceDimension());
			*_runLogFile<< "Problem : dim = "<<_Ndim<< " but mesh dim= "<<_mesh.getMeshDimension()<<", mesh space dim= "<<_mesh.getSpaceDimension()<<endl;
			*_runLogFile<<"DiffusionEquation::initialize: mesh has incorrect dimension"<<endl;
			_runLogFile->close();
			throw CdmathException("!!!!!!!!DiffusionEquation::initialize: mesh has incorrect  dimension");
		}
	
		if(!_initialDataSet)
			throw CdmathException("!!!!!!!!DiffusionEquation::initialize() set initial data first");
		else
	        {
	            PetscPrintf(PETSC_COMM_SELF,"\n Initialising the diffusion of a solid temperature using ");
	            *_runLogFile<<"Initialising the diffusion of a solid temperature using ";
	            if(!_FECalculation)
	            {
	                PetscPrintf(PETSC_COMM_SELF,"Finite volumes method\n\n");
	                *_runLogFile<< "Finite volumes method"<<endl<<endl;
	            }
	            else
	            {
	                PetscPrintf(PETSC_COMM_SELF,"Finite elements method\n\n");
	                *_runLogFile<< "Finite elements method"<<endl<<endl;
	            }
	        }

		/**************** Field creation *********************/
	
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
	
	    /* Détection des noeuds frontière avec une condition limite de Dirichlet */
	    if(_FECalculation)
	    {
	        if(_Ndim==1 )//The 1D cdmath mesh is necessarily made of segments
				PetscPrintf(PETSC_COMM_SELF,"1D Finite element method on segments\n");
	        else if(_Ndim==2)
	        {
				if( _mesh.isTriangular() )//Mesh dim=2
					PetscPrintf(PETSC_COMM_SELF,"2D Finite element method on triangles\n");
				else if (_mesh.getMeshDimension()==1)//Mesh dim=1
					PetscPrintf(PETSC_COMM_SELF,"1D Finite element method on a 2D network : space dimension is %d, mesh dimension is %d\n",_Ndim,_mesh.getMeshDimension());			
				else
				{
					PetscPrintf(PETSC_COMM_SELF,"Error Finite element with space dimension %, and mesh dimension  %d, mesh should be either triangular either 1D network\n",_Ndim,_mesh.getMeshDimension());
					*_runLogFile<<"DiffusionEquation::initialize mesh has incorrect dimension"<<endl;
					_runLogFile->close();
					throw CdmathException("DiffusionEquation::initialize: mesh has incorrect cell types");
				}
	        }
	        else if(_Ndim==3)
	        {
				if( _mesh.isTetrahedral() )//Mesh dim=3
					PetscPrintf(PETSC_COMM_SELF,"3D Finite element method on tetrahedra\n");
				else if (_mesh.getMeshDimension()==2 and _mesh.isTriangular())//Mesh dim=2
					PetscPrintf(PETSC_COMM_SELF,"2D Finite element method on a 3D surface : space dimension is %d, mesh dimension is %d\n",_Ndim,_mesh.getMeshDimension());			
				else if (_mesh.getMeshDimension()==1)//Mesh dim=1
					PetscPrintf(PETSC_COMM_SELF,"1D Finite element method on a 3D network : space dimension is %d, mesh dimension is %d\n",_Ndim,_mesh.getMeshDimension());			
				else
				{
					PetscPrintf(PETSC_COMM_SELF,"Error Finite element with space dimension %d, and mesh dimension  %d, mesh should be either tetrahedral, either a triangularised surface or 1D network",_Ndim,_mesh.getMeshDimension());
					*_runLogFile<<"DiffusionEquation::initialize mesh has incorrect dimension"<<endl;
					_runLogFile->close();
					throw CdmathException("DiffusionEquation::initialize: mesh has incorrect cell types");
				}
	        }
				
	        _boundaryNodeIds = _mesh.getBoundaryNodeIds();
	        _NboundaryNodes=_boundaryNodeIds.size();
	
	        if(_NboundaryNodes==_Nnodes)
	            PetscPrintf(PETSC_COMM_SELF,"!!!!! Warning : all nodes are boundary nodes !!!!!");
	
	        for(int i=0; i<_NboundaryNodes; i++)
	            if(_limitField[(_mesh.getNode(_boundaryNodeIds[i])).getGroupName()].bcType==DirichletDiffusion)
	                _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
	        _NdirichletNodes=_dirichletNodeIds.size();
	        _NunknownNodes=_Nnodes - _NdirichletNodes;
	        PetscPrintf(PETSC_COMM_SELF,"Number of unknown nodes %d, Number of boundary nodes %d, Number of Dirichlet boundary nodes %d\n\n", _NunknownNodes,_NboundaryNodes, _NdirichletNodes);
	        *_runLogFile<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
	    }
	}

    if(!_FECalculation)
		_globalNbUnknowns = _Nmailles*_nVar;
    else{
    	MPI_Bcast(&_NunknownNodes, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		_globalNbUnknowns = _NunknownNodes*_nVar;
	}

	/* Vectors creations */
	VecCreate(PETSC_COMM_WORLD, &_Tk);//main unknown
    VecSetSizes(_Tk,PETSC_DECIDE,_globalNbUnknowns);
	VecSetFromOptions(_Tk);
	VecGetLocalSize(_Tk, &_localNbUnknowns);
	
	VecDuplicate(_Tk, &_Tn);
	VecDuplicate(_Tk, &_Tkm1);
	VecDuplicate(_Tk, &_deltaT);
	VecDuplicate(_Tk, &_b);//RHS of the linear system: _b=Tn/dt + _b0 + puisance volumique + couplage thermique avec le fluide
	VecDuplicate(_Tk, &_b0);//part of the RHS that comes from the boundary conditions. Computed only once at the first time step

	if(_mpi_rank == 0)//Process 0 reads and distributes initial data
		if(_FECalculation)
			for(int i = 0; i<_NunknownNodes; i++)
	        {
				int globalIndex = globalNodeIndex(i, _dirichletNodeIds);
				VecSetValue(_Tn,i,_VV(globalIndex), INSERT_VALUES);
			}
		else
			for(int i = 0; i<_Nmailles; i++)
				VecSetValue( _Tn, i, _VV(i), INSERT_VALUES);
	VecAssemblyBegin(_Tn);
	VecAssemblyEnd(_Tn);
		
	/* Matrix creation */
   	MatCreateAIJ(PETSC_COMM_WORLD, _localNbUnknowns, _localNbUnknowns, _globalNbUnknowns, _globalNbUnknowns, _d_nnz, PETSC_NULL, _o_nnz, PETSC_NULL, &_A);

	/* Local sequential vector creation */
	if(_mpi_size>1 && _mpi_rank == 0)
		VecCreateSeq(PETSC_COMM_SELF, _globalNbUnknowns, &_Tn_seq);//For saving results on proc 0
	VecScatterCreateToZero(_Tn,&_scat,&_Tn_seq);

	//Linear solver
	KSPCreate(PETSC_COMM_WORLD, &_ksp);
	KSPSetType(_ksp, _ksptype);
	KSPGetPC(_ksp, &_pc);
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
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);

	_initializedMemory=true;
	save();//save initial data
}

double DiffusionEquation::computeTimeStep(bool & stop){
	if(!_diffusionMatrixSet)//The diffusion matrix is computed once and for all time steps
		_dt_diffusion=computeDiffusionMatrix(stop);

    //reset right hand side
   	VecCopy(_b0,_b);

	_dt_src=computeRHS(stop);

	VecAssemblyBegin(_b);          
	VecAssemblyEnd(  _b);

    if(_verbose or _system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Right hand side of the linear system\n");
        VecView(_b,PETSC_VIEWER_STDOUT_WORLD);
	}

	stop=false;
	return min(_dt_diffusion,_dt_src);
}

double DiffusionEquation::computeDiffusionMatrix(bool & stop)
{
    double result;
    
	MatZeroEntries(_A);
	VecZeroEntries(_b0);

    if(_FECalculation)
        result=computeDiffusionMatrixFE(stop);
    else
        result=computeDiffusionMatrixFV(stop);

	PetscPrintf(PETSC_COMM_WORLD,"Maximum diffusivity is %.2e, CFL = %.2f, Delta x = %.2e\n",_maxvp,_cfl,_minl);

    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);          
	VecAssemblyEnd(  _b0);


    //Contribution from the solid/fluid heat exchange with assumption of constant heat transfer coefficient
    //update value here if variable  heat transfer coefficient
    if(_timeScheme == Implicit and _heatTransfertCoeff/(_rho*_cp)>_precision)
        MatShift(_A,_heatTransfertCoeff/(_rho*_cp));//Contribution from the liquit/solid heat transfer
        
    if(_verbose or _system)
        MatView(_A,PETSC_VIEWER_STDOUT_WORLD);

    return  result;
}

double DiffusionEquation::computeDiffusionMatrixFE(bool & stop){

	if(_mpi_rank == 0)
	{
		Cell Cj;
		string nameOfGroup;
		double coeff;//Diffusion coefficients between nodes i and j
	    
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
	                        MatSetValue(_A,i_int,j_int,(_DiffusionTensor*GradShapeFuncs[idim])*GradShapeFuncs[jdim]/Cj.getMeasure(), ADD_VALUES);
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
	                        coeff =-1.*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
	                        VecSetValue(_b0,i_int,coeff, ADD_VALUES);                        
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
	                    j_int=unknownNodeIndex(Fi.getNodeId(j), _dirichletNodeIds);//indice du noeud j en tant que noeud inconnu
	                    if( _neumannValuesSet )
	                        coeff =Fi.getMeasure()/_Ndim*_neumannBoundaryValues[Fi.getNodeId(j)];
	                    else
	                        coeff =Fi.getMeasure()/_Ndim*_limitField[_mesh.getNode(Fi.getNodeId(j)).getGroupName()].normalFlux;
	                    VecSetValue(_b0, j_int, coeff, ADD_VALUES);
	                }
	            }
	        }
	    }
	}

	_diffusionMatrixSet=true;

	_maxvp=_diffusivity;//To do : optimise value with the mesh while respecting stability
	if(fabs(_maxvp)<_precision)
		throw CdmathException("DiffusionEquation::computeDiffusionMatrixFE(): Error computing time step ! Maximum diffusivity is zero => division by zero");
	else
		return _cfl*_minl*_minl/_maxvp;
}

double DiffusionEquation::computeDiffusionMatrixFV(bool & stop){
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
			if(fabs(dn)>_maxvp)
				_maxvp=fabs(dn);
	
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
				nameOfGroup = Fj.getGroupName();
	
				if (_limitField[nameOfGroup].bcType==NeumannDiffusion){
	                VecSetValue(_b0,idm,   -dn*inv_dxi*_limitField[nameOfGroup].normalFlux, ADD_VALUES);
				}
				else if(_limitField[nameOfGroup].bcType==DirichletDiffusion){
					barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
					MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                           , ADD_VALUES);
					VecSetValue(_b0,idm,    dn*inv_dxi/barycenterDistance*_limitField[nameOfGroup].T, ADD_VALUES);
				}
				else {
	                stop=true ;
					cout<<"!!!!!!!!!!!!!!!!! Error DiffusionEquation::computeDiffusionMatrixFV !!!!!!!!!!"<<endl;
	                cout<<"Boundary condition not accepted for boundary named "<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
					cout<<"Accepted boundary conditions are NeumannDiffusion "<<NeumannDiffusion<< " and DirichletDiffusion "<<DirichletDiffusion<<endl;
	                *_runLogFile<<"Boundary condition not accepted for boundary named "<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
					throw CdmathException("Boundary condition not accepted");
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
				throw CdmathException("DiffusionEquation::computeDiffusionMatrixFV(): incompatible number of cells around a face");
		}
	}

	_diffusionMatrixSet=true;

	MPI_Bcast(&_maxvp, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);//The determination of _maxvp is optimal, unlike in the FE case

	if(fabs(_maxvp)<_precision)
		throw CdmathException("DiffusionEquation::computeDiffusionMatrixFV(): Error computing time step ! Maximum diffusivity is zero => division by zero");
	else
		return _cfl*_minl*_minl/_maxvp;
}

double DiffusionEquation::computeRHS(bool & stop){//Contribution of the PDE RHS to the linear systemm RHS (boundary conditions do contribute to the system RHS via the function computeDiffusionMatrix

	if(_mpi_rank == 0)
	{
	    double Ti;  
	    if(!_FECalculation)
	        for (int i=0; i<_Nmailles;i++)
	            //Contribution due to fluid/solide heat exchange + Contribution of the volumic heat power
	            if(_timeScheme == Explicit)
	            {
	                VecGetValues(_Tn, 1, &i, &Ti);
	                VecSetValue(_b,i,(_heatTransfertCoeff*(_fluidTemperatureField(i)-Ti)+_heatPowerField(i))/(_rho*_cp),ADD_VALUES);
	            }
	            else//Implicit scheme    
	                VecSetValue(_b,i,(_heatTransfertCoeff* _fluidTemperatureField(i)    +_heatPowerField(i))/(_rho*_cp)    ,ADD_VALUES);
	    else
	        {
	            Cell Ci;
	            std::vector< int > nodesId;
	            for (int i=0; i<_Nmailles;i++)
	            {
	                Ci=_mesh.getCell(i);
	                nodesId=Ci.getNodesId();
	                for (int j=0; j<nodesId.size();j++)
	                    if(!_mesh.isBorderNode(nodesId[j])) //or for better performance nodeIds[idim]>dirichletNodes.upper_bound()
	                    {
	                        double coeff = (_heatTransfertCoeff*_fluidTemperatureField(nodesId[j]) + _heatPowerField(nodesId[j]))/(_rho*_cp);
	                        VecSetValue(_b,unknownNodeIndex(nodesId[j], _dirichletNodeIds), coeff*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
	                    }
	            }
	        }
	}

    if(_heatTransfertCoeff>_precision)
        return _rho*_cp/_heatTransfertCoeff;
    else
        return INFINITY;
}

bool DiffusionEquation::initTimeStep(double dt){

	/* tricky because of code coupling */
    if(_dt>0 and dt>0)//Previous time step was set and used
    {
        //Remove the contribution from dt to prepare for new time step. The diffusion matrix is not recomputed
        if(_timeScheme == Implicit)
            MatShift(_A,-1/_dt+1/dt);
        //No need to remove the contribution to the right hand side since it is recomputed from scratch at each time step
    }
    else if(dt>0)//_dt==0, first time step
    {
        if(_timeScheme == Implicit)
            MatShift(_A,1/dt);        
    }
    else//dt<=0
    {
        PetscPrintf(PETSC_COMM_WORLD,"DiffusionEquation::initTimeStep %.2e = \n",dt);
        throw CdmathException("Error DiffusionEquation::initTimeStep : cannot set time step to zero");        
    }
    //At this stage _b contains _b0 + power + heat exchange
    VecAXPY(_b, 1/dt, _Tn);        

	_dt = dt;

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Matrix of the linear system\n");
		MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
	}
	
	return _dt>0;
}

void DiffusionEquation::abortTimeStep(){
    //Remove contribution of dt to the RHS
	VecAXPY(_b,  -1/_dt, _Tn);
    //Remove contribution of dt to the matrix
	if(_timeScheme == Implicit)
		MatShift(_A,-1/_dt);
	_dt = 0;
}

bool DiffusionEquation::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		_maxvp=0;
		computeTimeStep(stop);
	}
	if(stop){
		converged=false;
		return false;
	}

	if(_timeScheme == Explicit)
	{
		MatMult(_A, _Tn, _Tk);
		VecAXPY(_Tk, -1, _b);
		VecScale(_Tk, -_dt);

		converged = true;
	}
	else
	{

#if PETSC_VERSION_GREATER_3_5
		KSPSetOperators(_ksp, _A, _A);
#else
		KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

		if(_conditionNumber)
			KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
		KSPSolve(_ksp, _b, _Tk);

		KSPGetIterationNumber(_ksp, &_PetscIts);
		if( _MaxIterLinearSolver < _PetscIts)
			_MaxIterLinearSolver = _PetscIts;
		if(_PetscIts>=_maxPetscIts)
		{
			PetscPrintf(PETSC_COMM_WORLD,"Systeme lineaire : pas de convergence de Petsc. Itérations maximales %d atteintes \n",_maxPetscIts);
			*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
			converged=false;
			return false;
		}
		else
		{
			VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
			VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1
			VecNorm(_deltaT,NORM_INFINITY,&_erreur_rel);
			converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
		}
	}

	VecCopy(_Tk, _Tkm1);

	return true;
}
void DiffusionEquation::validateTimeStep()
{
	VecCopy(_Tk, _deltaT);//ici Tk=Tnp1 donc on a deltaT=Tnp1
	VecAXPY(_deltaT,  -1, _Tn);//On obtient deltaT=Tnp1-Tn
	VecNorm(_deltaT,NORM_INFINITY,&_erreur_rel);

	_isStationary =(_erreur_rel <_precision);

	VecCopy(_Tk, _Tn);
	VecCopy(_Tk, _Tkm1);

	_time+=_dt;
	_nbTimeStep++;
	
	if ((_nbTimeStep-1)%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
		save();
}

void DiffusionEquation::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string resultFile(_path+"/DiffusionEquation");//Results

	resultFile+="_";
	resultFile+=_fileName;

	if(_mpi_size>1){
		VecScatterBegin(_scat,_Tn,_Tn_seq,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(  _scat,_Tn,_Tn_seq,INSERT_VALUES,SCATTER_FORWARD);
	}
	
    if(_verbose or _system)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Unknown of the linear system :\n");
        VecView(_Tn,PETSC_VIEWER_STDOUT_WORLD);
	}

	if(_mpi_rank==0){
	    //On remplit le champ
	    double Ti;
	    if(!_FECalculation)
	        for(int i =0; i<_Nmailles;i++)
	        {
				if(_mpi_size>1)
					VecGetValues(_Tn_seq, 1, &i, &Ti);
				else
					VecGetValues(_Tn    , 1, &i, &Ti);
					
	            _VV(i)=Ti;
	        }
	    else
	    {
	        int globalIndex;
	        for(int i=0; i<_NunknownNodes; i++)
	        {
				if(_mpi_size>1)
					VecGetValues(_Tn_seq, 1, &i, &Ti);
				else
					VecGetValues(_Tk    , 1, &i, &Ti);
	            globalIndex = globalNodeIndex(i, _dirichletNodeIds);
	            _VV(globalIndex)=Ti;
	        }
	
	        Node Ni;
	        string nameOfGroup;
	        for(int i=0; i<_NdirichletNodes; i++)//Assumes node numbering starts with border nodes
	        {
	            Ni=_mesh.getNode(_dirichletNodeIds[i]);
	            nameOfGroup = Ni.getGroupName();
	            _VV(_dirichletNodeIds[i])=_limitField[nameOfGroup].T;
	        }
	    }
		_VV.setTime(_time,_nbTimeStep);
	
		// create mesh and component info
		if (_nbTimeStep ==0 || _restartWithNewFileName){
			if (_restartWithNewFileName)
				_restartWithNewFileName=false;
			string suppress ="rm -rf "+resultFile+"_*";
			system(suppress.c_str());//Nettoyage des précédents calculs identiques
	        
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
		else{	// do not create mesh
			switch(_saveFormat)
			{
			case VTK :
				_VV.writeVTK(resultFile,false);
				break;
			case MED :
				_VV.writeMED(resultFile,false);
				break;
			case CSV :
				_VV.writeCSV(resultFile);
				break;
			}
		}
	    
	    if(_isStationary)
		{
	        resultFile+="_Stat";
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
}

void DiffusionEquation::terminate(){
	VecDestroy(&_Tn);
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b0);
	VecDestroy(&_b);
	MatDestroy(&_A);
	if(_mpi_size>1 && _mpi_rank == 0)
		VecDestroy(&_Tn_seq);
}

void 
DiffusionEquation::setDirichletValues(map< int, double> dirichletBoundaryValues)
{
    _dirichletValuesSet=true;
    _dirichletBoundaryValues=dirichletBoundaryValues;
}

void 
DiffusionEquation::setNeumannValues(map< int, double> neumannBoundaryValues)
{
    _neumannValuesSet=true;
    _neumannBoundaryValues=neumannBoundaryValues;
}


Field& 
DiffusionEquation::getOutputTemperatureField()
{
    if(!_initializedMemory)
        throw("Computation not initialized. No temperature field available");
    else
        return _VV;
}

Field& 
DiffusionEquation::getRodTemperatureField()
{
   return getOutputTemperatureField();
}

vector<string> 
DiffusionEquation::getInputFieldsNames()
{
	vector<string> result(2);
	
	result[0]="FluidTemperature";
	result[1]="HeatPower";
	
	return result;
}
vector<string> 
DiffusionEquation::getOutputFieldsNames()
{
	vector<string> result(1);
	
	result[0]="RodTemperature";
	
	return result;
}

Field& 
DiffusionEquation::getOutputField(const string& nameField )
{
	if(nameField=="RodTemperature" || nameField=="RODTEMPERATURE" || nameField=="TEMPERATURECOMBUSTIBLE" || nameField=="TemperatureCombustible" )
		return getRodTemperatureField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first" << endl;
        throw CdmathException("DiffusionEquation::getOutputField error : Unknown Field name");
    }
}

void
DiffusionEquation::setInputField(const string& nameField, Field& inputField )
{
	if(!_initialDataSet)
		throw CdmathException("!!!!!!!! DiffusionEquation::setInputField set initial field first");

	if(nameField=="FluidTemperature" || nameField=="FLUIDTEMPERATURE" || nameField=="TemperatureFluide" || nameField=="TEMPERATUREFLUIDE")
		return setFluidTemperatureField( inputField) ;
	else if(nameField=="HeatPower" || nameField=="HEATPOWER" || nameField=="PuissanceThermique" || nameField=="PUISSANCETHERMIQUE" )
		return setHeatPowerField( inputField );
	else
    {
        cout<<"Error : Field name "<< nameField << " is not an input field name, call getInputFieldsNames first" << endl;
        throw CdmathException("DiffusionEquation::setInputField error : Unknown Field name");
    }
}

void 
DiffusionEquation::setFluidTemperatureField(Field coupledTemperatureField){
	if(!_initialDataSet)
		throw CdmathException("!!!!!!!! DiffusionEquation::setFluidTemperatureField set initial field first");

	coupledTemperatureField.getMesh().checkFastEquivalWith(_mesh);
	_fluidTemperatureField=coupledTemperatureField;
	_fluidTemperatureFieldSet=true;
};
