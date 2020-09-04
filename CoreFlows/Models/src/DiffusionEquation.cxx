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
    else //unknownNodeMax>=unknownNodeIndex) hence our node global number is between dirichletNodes[j-1] and dirichletNodes[j]
        return unknownNodeIndex - unknownNodeMax + dirichletNodes[j]-1;
}

Vector DiffusionEquation::gradientNodal(Matrix M, vector< double > values){
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

DiffusionEquation::DiffusionEquation(int dim, bool FECalculation,double rho,double cp, double lambda){
    /* Control input value are acceptable */
    if(rho<_precision or cp<_precision)
    {
        std::cout<<"rho="<<rho<<", cp= "<<cp<< ", precision= "<<_precision;
        throw CdmathException("Error : parameters rho and cp should be strictly positive");
    }
    if(lambda < 0.)
    {
        std::cout<<"conductivity="<<lambda<<endl;
        throw CdmathException("Error : conductivity parameter lambda cannot  be negative");
    }
    if(dim<=0)
    {
        std::cout<<"space dimension="<<dim<<endl;
        throw CdmathException("Error : parameter dim cannot  be negative");
    }

    cout<<"Diffusion problem with density "<<rho<<", specific heat "<< cp<<", conductivity "<< lambda;
    if(FECalculation)
        cout<<" and finite elements method"<<endl;
    else
        cout<<" and finite volumes method"<<endl;
    
    _FECalculation=FECalculation;
    
    /* Finite element data */
    _neibMaxNbNodes=0;    
    _boundaryNodeIds=std::vector< int >(0);
    _dirichletNodeIds=std::vector< int >(0);
    _NboundaryNodes=0;
    _NdirichletNodes=0;
    _NunknownNodes=0;

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

	_fileName = "CoreFlowsDiffusionProblem";

	_runLogFile=new ofstream;

    /* Default diffusion tensor is identity matrix */
   	_DiffusionTensor=Matrix(_Ndim);
	for(int idim=0;idim<_Ndim;idim++)
		_DiffusionTensor(idim,idim)=1;
}

void DiffusionEquation::initialize()
{
	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation

	if(_Ndim != _mesh.getSpaceDimension() or _Ndim!=_mesh.getMeshDimension())//for the moment we must have space dim=mesh dim
	{
        cout<< "Problem : dimension defined is "<<_Ndim<< " but mesh dimension= "<<_mesh.getMeshDimension()<<", and space dimension is "<<_mesh.getSpaceDimension()<<endl;
		*_runLogFile<< "Problem : dim = "<<_Ndim<< " but mesh dim= "<<_mesh.getMeshDimension()<<", mesh space dim= "<<_mesh.getSpaceDimension()<<endl;
		*_runLogFile<<"DiffusionEquation::initialize: mesh has incorrect dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("DiffusionEquation::initialize: mesh has incorrect  dimension");
	}

	if(!_initialDataSet)
		throw CdmathException("DiffusionEquation::initialize() set initial data first");
	else
        {
            cout<<"Initialising the diffusion of a solid temperature using ";
            *_runLogFile<<"Initialising the diffusion of a solid temperature using ";
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
				*_runLogFile<<"DiffusionEquation::initialize mesh has incorrect dimension"<<endl;
				_runLogFile->close();
				throw CdmathException("DiffusionEquation::initialize: mesh has incorrect cell types");
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
				*_runLogFile<<"DiffusionEquation::initialize mesh has incorrect dimension"<<endl;
				_runLogFile->close();
				throw CdmathException("DiffusionEquation::initialize: mesh has incorrect cell types");
			}
        }
			
        _boundaryNodeIds = _mesh.getBoundaryNodeIds();
        _NboundaryNodes=_boundaryNodeIds.size();

        if(_NboundaryNodes==_Nnodes)
            cout<<"!!!!! Warning : all nodes are boundary nodes !!!!!"<<endl;

        for(int i=0; i<_NboundaryNodes; i++)
            if(_limitField[(_mesh.getNode(_boundaryNodeIds[i])).getGroupName()].bcType==DirichletDiffusion)
                _dirichletNodeIds.push_back(_boundaryNodeIds[i]);
        _NdirichletNodes=_dirichletNodeIds.size();
        _NunknownNodes=_Nnodes - _NdirichletNodes;
        cout<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
        *_runLogFile<<"Number of unknown nodes " << _NunknownNodes <<", Number of boundary nodes " << _NboundaryNodes<< ", Number of Dirichlet boundary nodes " << _NdirichletNodes <<endl<<endl;
    }

	//creation de la matrice
    if(!_FECalculation)
        MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles, _Nmailles, (1+_neibMaxNb), PETSC_NULL, &_A);
    else
        MatCreateSeqAIJ(PETSC_COMM_SELF, _NunknownNodes, _NunknownNodes, (1+_neibMaxNbNodes), PETSC_NULL, &_A);

	VecCreate(PETSC_COMM_SELF, &_Tk);

    if(!_FECalculation)
        VecSetSizes(_Tk,PETSC_DECIDE,_Nmailles);
    else
        VecSetSizes(_Tk,PETSC_DECIDE,_NunknownNodes);

	VecSetFromOptions(_Tk);
	VecDuplicate(_Tk, &_Tn);
	VecDuplicate(_Tk, &_Tkm1);
	VecDuplicate(_Tk, &_deltaT);
	VecDuplicate(_Tk, &_b);//RHS of the linear system: _b=Tn/dt + _b0 + puisance volumique + couplage thermique avec le fluide
	VecDuplicate(_Tk, &_b0);//part of the RHS that comes from the boundary conditions. Computed only once at the first time step

    if(!_FECalculation)
        for(int i =0; i<_Nmailles;i++)
            VecSetValue(_Tn,i,_VV(i), INSERT_VALUES);
    else
        for(int i =0; i<_Nnodes;i++)
            VecSetValue(_Tn,i,_VV(i), INSERT_VALUES);

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

	_initializedMemory=true;
	save();//save initial data
}

double DiffusionEquation::computeTimeStep(bool & stop){
	if(!_diffusionMatrixSet)//The diffusion matrix is computed once and for all time steps
		_dt_diffusion=computeDiffusionMatrix(stop);

    //reset right hand side
   	VecCopy(_b0,_b);

	_dt_src=computeRHS(stop);

	stop=false;
	return min(_dt_diffusion,_dt_src);
}

double DiffusionEquation::computeDiffusionMatrix(bool & stop)
{
    double result;
    
    if(_FECalculation)
        result=computeDiffusionMatrixFE(stop);
    else
        result=computeDiffusionMatrixFV(stop);

    //Contribution from the solid/fluid heat exchange with assumption of constant heat transfer coefficient
    //update value here if variable  heat transfer coefficient
    if(_timeScheme == Implicit and _heatTransfertCoeff/(_rho*_cp)>_precision)
        MatShift(_A,_heatTransfertCoeff/(_rho*_cp));//Contribution from the liquit/solid heat transfer
        
    if(_verbose or _system)
        MatView(_A,PETSC_VIEWER_STDOUT_SELF);

    return  result;
}

double DiffusionEquation::computeDiffusionMatrixFE(bool & stop){
	Cell Cj;
	string nameOfGroup;
	double dij;//Diffusion coefficients between nodes i and j
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
                        dij=_conductivity*(_DiffusionTensor*GradShapeFuncs[idim])*GradShapeFuncs[jdim]/Cj.getMeasure();
                        MatSetValue(_A,i_int,j_int,dij, ADD_VALUES);
                        if(fabs(dij)>_maxvp)
                            _maxvp=fabs(dij);
                    }
                    else if (!dirichletCell_treated)
                    {//Second node of the edge is a Dirichlet node
                        dirichletCell_treated=true;
                        for (int kdim=0; kdim<_Ndim+1;kdim++)
                        {
                            if(find(_dirichletNodeIds.begin(),_dirichletNodeIds.end(),nodeIds[kdim])!=_dirichletNodeIds.end())
                                valuesBorder[kdim]=_limitField[nameOfGroup].T;
                            else
                                valuesBorder[kdim]=0;                            
                        }
                        GradShapeFuncBorder=gradientNodal(M,valuesBorder)/fact(_Ndim);
                        dij =-_conductivity*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
                        VecSetValue(_b,i_int,dij, ADD_VALUES);                        
                    }
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

	cout<<"_maxvp= "<<_maxvp<< " cfl= "<<_cfl<<" minl= "<<_minl<<endl;
	if(fabs(_maxvp)<_precision)
		throw CdmathException("DiffusionEquation::computeDiffusionMatrix(): maximum eigenvalue for time step is zero");
	else
		return _cfl/_maxvp;
}

double DiffusionEquation::computeDiffusionMatrixFV(bool & stop){
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
	VecZeroEntries(_b0);
    
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
                VecSetValue(_b,idm,   -dn*inv_dxi*_limitField[nameOfGroup].normalFlux, ADD_VALUES);
			}
			else if(_limitField[nameOfGroup].bcType==DirichletDiffusion){
				barycenterDistance=Cell1.getBarryCenter().distance(Fj.getBarryCenter());
				MatSetValue(_A,idm,idm,dn*inv_dxi/barycenterDistance                           , ADD_VALUES);
				VecSetValue(_b,idm,    dn*inv_dxi/barycenterDistance*_limitField[nameOfGroup].T, ADD_VALUES);
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

	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);
	VecAssemblyEnd(_b0);

	_diffusionMatrixSet=true;
    stop=false ;

	cout<<"_maxvp= "<<_maxvp<< " cfl= "<<_cfl<<" minl= "<<_minl<<endl;
	if(fabs(_maxvp)<_precision)
		throw CdmathException("DiffusionEquation::computeDiffusionMatrixFV(): maximum eigenvalue for time step is zero");
	else
		return _cfl*_minl*_minl/_maxvp;
}

double DiffusionEquation::computeRHS(bool & stop){
	VecAssemblyBegin(_b);          
    double Ti;  
    if(!_FECalculation)
        for (int i=0; i<_Nmailles;i++)
        {
            VecSetValue(_b,i,_heatPowerField(i)/(_rho*_cp),ADD_VALUES);//Contribution of the volumic heat power
            //Contribution due to fluid/solide heat exchange
            if(_timeScheme == Explicit)
            {
                VecGetValues(_Tn, 1, &i, &Ti);
                VecSetValue(_b,i,_heatTransfertCoeff/(_rho*_cp)*(_fluidTemperatureField(i)-Ti),ADD_VALUES);
            }
            else//Implicit scheme    
                VecSetValue(_b,i,_heatTransfertCoeff/(_rho*_cp)* _fluidTemperatureField(i)    ,ADD_VALUES);
        }
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
                        double coeff = _heatTransfertCoeff*_fluidTemperatureField(nodesId[j]) + _heatPowerField(nodesId[j]);
                        VecSetValue(_b,unknownNodeIndex(nodesId[j], _dirichletNodeIds), coeff*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                    }
            }
        }
	VecAssemblyEnd(_b);

    if(_verbose or _system)
        VecView(_b,PETSC_VIEWER_STDOUT_SELF);

    stop=false ;
    if(_heatTransfertCoeff>_precision)
        return _rho*_cp/_heatTransfertCoeff;
    else
        return INFINITY;
}

bool DiffusionEquation::initTimeStep(double dt){

    if(_dt>0 and dt>0)
    {
        //Remove the contribution from dt to prepare for new initTimeStep. The diffusion matrix is not recomputed
        if(_timeScheme == Implicit)
            MatShift(_A,-1/_dt+1/dt);
        //No need to remove the contribution to the right hand side since it is recomputed from scratch at each time step
    }
    else if(dt>0)//_dt==0, first time step
    {
        MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
        if(_timeScheme == Implicit)
            MatShift(_A,1/_dt);        
    }
    else//dt<=0
    {
        cout<<"DiffusionEquation::initTimeStep dt= "<<dt<<endl;
        throw CdmathException("Error DiffusionEquation::initTimeStep : cannot set time step to zero");        
    }
    //At this stage _b contains _b0 + power + heat exchange
    VecAXPY(_b, 1/_dt, _Tn);        

	_dt = dt;

	if(_verbose && _nbTimeStep%_freqSave ==0)
		MatView(_A,PETSC_VIEWER_STDOUT_SELF);

	return _dt>0;
}

void DiffusionEquation::abortTimeStep(){
    //Remove contribution od dt to the RHS
	VecAXPY(_b,  -1/_dt, _Tn);
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
    //Remove contribution od dt to the matrix
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
		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

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
			cout<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
			*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
			converged=false;
			return false;
		}
		else{
			VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
			VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1
			_erreur_rel= 0;
			double Ti, dTi;

            if(!_FECalculation)
                for(int i=0; i<_Nmailles; i++)
                {
                    VecGetValues(_deltaT, 1, &i, &dTi);
                    VecGetValues(_Tk, 1, &i, &Ti);
                    if(_erreur_rel < fabs(dTi/Ti))
                        _erreur_rel = fabs(dTi/Ti);
                }
            else
                for(int i=0; i<_Nnodes; i++)
                {
                    VecGetValues(_deltaT, 1, &i, &dTi);
                    VecGetValues(_Tk, 1, &i, &Ti);
                    if(_erreur_rel < fabs(dTi/Ti))
                        _erreur_rel = fabs(dTi/Ti);
                }
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

	_erreur_rel= 0;
	double Ti, dTi;

    if(!_FECalculation)
        for(int i=0; i<_Nmailles; i++)
        {
            VecGetValues(_deltaT, 1, &i, &dTi);
            VecGetValues(_Tk, 1, &i, &Ti);
            if(_erreur_rel < fabs(dTi/Ti))
                _erreur_rel = fabs(dTi/Ti);
        }
    else
        for(int i=0; i<_Nnodes; i++)
        {
            VecGetValues(_deltaT, 1, &i, &dTi);
            VecGetValues(_Tk, 1, &i, &Ti);
            if(_erreur_rel < fabs(dTi/Ti))
                _erreur_rel = fabs(dTi/Ti);
        }

	_isStationary =(_erreur_rel <_precision);

	VecCopy(_Tk, _Tn);
	VecCopy(_Tk, _Tkm1);

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout <<"Valeur propre locale max: " << _maxvp << endl;

	_time+=_dt;
	_nbTimeStep++;
	if (_nbTimeStep%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
        save();
}

void DiffusionEquation::save(){
    cout<< "Saving numerical results"<<endl<<endl;
    *_runLogFile<< "Saving numerical results"<< endl<<endl;

	string resultFile(_path+"/DiffusionEquation");//Results

	resultFile+="_";
	resultFile+=_fileName;

    if(_verbose or _system)
        VecView(_Tk,PETSC_VIEWER_STDOUT_SELF);

    //On remplit le champ
    double Ti;
    if(!_FECalculation)
        for(int i =0; i<_Nmailles;i++)
        {
            VecGetValues(_Tn, 1, &i, &Ti);
            _VV(i)=Ti;
        }
    else
    {
        int globalIndex;
        for(int i=0; i<_NunknownNodes; i++)
        {
            VecGetValues(_Tk, 1, &i, &Ti);
            globalIndex = globalNodeIndex(i, _dirichletNodeIds);
            _VV(globalIndex)=Ti;//Assumes node numbering starts with border nodes
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

void DiffusionEquation::terminate(){
	VecDestroy(&_Tn);
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b0);
	VecDestroy(&_b);
	MatDestroy(&_A);
}

vector<string> DiffusionEquation::getOutputFieldsNames()
{
	vector<string> result(2);
	
	result[0]="FluidTemperature";
	result[1]="RodTemperature";
	
	return result;
}

Field& DiffusionEquation::getOutputField(const string& nameField )
{
	if(nameField=="FluidTemperature" || nameField=="FLUIDTEMPERATURE" || nameField=="TemperatureFluide" || nameField=="TEMPERATUREFLUIDE" )
		return getFluidTemperatureField();
	else if(nameField=="RodTemperature" || nameField=="RODTEMPERATURE" || nameField=="TEMPERATURECOMBUSTIBLE" || nameField=="TemperatureCombustible" )
		return getRodTemperatureField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first" << endl;
        throw CdmathException("DiffusionEquation::getOutputField error : Unknown Field name");
    }
}

