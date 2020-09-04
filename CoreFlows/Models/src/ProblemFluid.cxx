#include "ProblemFluid.hxx"
#include "Fluide.h"
#include "math.h"

using namespace std;

ProblemFluid::ProblemFluid(void){
	_latentHeat=1e30;
	_Tsat=1e30;
	_Psat=-1e30;
	_dHsatl_over_dp=0;
	_porosityFieldSet=false;
	_pressureLossFieldSet=false;
	_sectionFieldSet=false;
	_GravityField3d=vector<double>(3,0);
	_gravityReferencePoint=vector<double>(3,0);
	_Uroe=NULL;_Udiff=NULL;_temp=NULL;_l=NULL;_vec_normal=NULL;;
	_idm=NULL;_idn=NULL;
	_saveVelocity=false;
	_saveConservativeField=false;
	_usePrimitiveVarsInNewton=false;
	_saveInterfacialField=false;
	_err_press_max=0; _part_imag_max=0; _nbMaillesNeg=0; _nbVpCplx=0;_minm1=1e30;_minm2=1e30;//Pour affichage paramètres diphasiques dans IterateTimeStep
	_isScaling=false;
	_entropicCorrection=false;
	_pressureCorrectionOrder=2;
	_nonLinearFormulation=Roe;
	_maxvploc=0.;
	_spaceScheme=upwind;
}

SpaceScheme ProblemFluid::getSpaceScheme()
{
	return _spaceScheme;
}
void ProblemFluid::setNumericalScheme(SpaceScheme spaceScheme, TimeScheme timeScheme)
{
	if( _nbTimeStep>0 && timeScheme!=_timeScheme)//This is a change of time scheme during a simulation
		_restartWithNewTimeScheme=true;
	_timeScheme = timeScheme;
	_spaceScheme = spaceScheme;
}

void ProblemFluid::initialize()
{
	if(!_initialDataSet){
		*_runLogFile<<"ProblemFluid::initialize() set initial data first"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemFluid::initialize() set initial data first");
	}
	cout << "Number of Phases = " << _nbPhases << " mesh dimension = "<<_Ndim<<" number of variables = "<<_nVar<<endl;
	*_runLogFile << "Number of Phases = " << _nbPhases << " spaceDim= "<<_Ndim<<" number of variables= "<<_nVar<<endl;

	/********* local arrays ****************/
	_AroePlus = new PetscScalar[_nVar*_nVar];
	_Diffusion = new PetscScalar[_nVar*_nVar];
	_AroeMinus = new PetscScalar[_nVar*_nVar];
	_Aroe = new PetscScalar[_nVar*_nVar];
	_absAroe = new PetscScalar[_nVar*_nVar];
	_signAroe = new PetscScalar[_nVar*_nVar];
	_invAroe = new PetscScalar[_nVar*_nVar];
	_AroeImplicit = new PetscScalar[_nVar*_nVar];
	_AroeMinusImplicit = new PetscScalar[_nVar*_nVar];
	_AroePlusImplicit = new PetscScalar[_nVar*_nVar];
	_absAroeImplicit = new PetscScalar[_nVar*_nVar];
	_primToConsJacoMat = new PetscScalar[_nVar*_nVar];
	_phi = new PetscScalar[_nVar];
	_Jcb = new PetscScalar[_nVar*_nVar];
	_JcbDiff = new PetscScalar[_nVar*_nVar];
	_a = new PetscScalar[_nVar*_nVar];
	_externalStates = new PetscScalar[_nVar];
	_Ui = new PetscScalar[_nVar];
	_Uj = new PetscScalar[_nVar];
	_Vi = new PetscScalar[_nVar];
	_Vj = new PetscScalar[_nVar];
	_Uij = new PetscScalar[_nVar];//used for VFRoe scheme
	_Vij = new PetscScalar[_nVar];//used for VFRoe scheme
	_Si = new PetscScalar[_nVar];
	_Sj = new PetscScalar[_nVar];
	_pressureLossVector= new PetscScalar[_nVar];
	_porosityGradientSourceVector= new PetscScalar[_nVar];

	_l = new double[_nVar];
	_r = new double[_nVar];
	_Udiff = new double[_nVar];
	_temp = new double[_nVar*_nVar];

	_idm = new int[_nVar];
	_idn = new int[_nVar];
	_vec_normal = new double[_Ndim];

	for(int k=0;k<_nVar*_nVar;k++)
	{
		_Diffusion[k]=0;
		_JcbDiff[k]=0;
	}
	for(int k=0; k<_nVar; k++){
		_idm[k] = k;
		_idn[k] = k;
	}

	/**************** Field creation *********************/
	if(!_heatPowerFieldSet){
		_heatPowerField=Field("Heat Power",CELLS,_mesh,1);
		for(int i =0; i<_Nmailles; i++)
			_heatPowerField(i) = _heatSource;
	}
	if(_Psat>-1e30)
		_dp_over_dt=Field("dP/dt",CELLS,_mesh,1);
	if(!_porosityFieldSet){
		_porosityField=Field("Porosity field",CELLS,_mesh,1);
		for(int i =0; i<_Nmailles; i++)
			_porosityField(i) = 1;
		_porosityFieldSet=true;
	}

	//conservative field used only for saving results
	_UU=Field ("Conservative vec", CELLS, _mesh, _nVar);
	if(_saveInterfacialField && _nonLinearFormulation==VFRoe)
	{
		_UUstar=Field ("Interfacial U", CELLS, _mesh, _nVar);
		_VVstar=Field ("Interfacial V", CELLS, _mesh, _nVar);
	}

	//Construction des champs primitifs et conservatifs initiaux comme avant dans ParaFlow
	double * initialFieldPrim = new double[_nVar*_Nmailles];
	for(int i =0; i<_Nmailles;i++)
		for(int j =0; j<_nVar;j++)
			initialFieldPrim[i*_nVar+j]=_VV(i,j);
	double *initialFieldCons = new double[_nVar*_Nmailles];
	for(int i=0; i<_Nmailles; i++)
		primToCons(initialFieldPrim, i, initialFieldCons, i);
	for(int i =0; i<_Nmailles;i++)
		for(int j =0; j<_nVar;j++)
			_UU(i,j)=initialFieldCons[i*_nVar+j];

	/**********Petsc structures:  ****************/

	//creation de la matrice
	if(_timeScheme == Implicit)
		MatCreateSeqBAIJ(PETSC_COMM_SELF, _nVar, _nVar*_Nmailles, _nVar*_Nmailles, (1+_neibMaxNb), PETSC_NULL, &_A);

	//creation des vecteurs
	VecCreateSeq(PETSC_COMM_SELF, _nVar, &_Uext);
	VecCreateSeq(PETSC_COMM_SELF, _nVar, &_Vext);
	VecCreateSeq(PETSC_COMM_SELF, _nVar, &_Uextdiff);
	//	  VecCreateSeq(PETSC_COMM_SELF, _nVar*_Nmailles, &_conservativeVars);
	VecCreate(PETSC_COMM_SELF, &_conservativeVars);
	VecSetSizes(_conservativeVars,PETSC_DECIDE,_nVar*_Nmailles);
	VecSetBlockSize(_conservativeVars,_nVar);
	VecSetFromOptions(_conservativeVars);
	VecDuplicate(_conservativeVars, &_old);
	VecDuplicate(_conservativeVars, &_newtonVariation);
	VecDuplicate(_conservativeVars, &_b);
	VecDuplicate(_conservativeVars, &_primitiveVars);

	if(_isScaling)
	{
		VecDuplicate(_conservativeVars, &_vecScaling);
		VecDuplicate(_conservativeVars, &_invVecScaling);
		VecDuplicate(_conservativeVars, &_bScaling);

		_blockDiag = new PetscScalar[_nVar];
		_invBlockDiag = new PetscScalar[_nVar];
	}


	int *indices = new int[_Nmailles];
	for(int i=0; i<_Nmailles; i++)
		indices[i] = i;
	VecSetValuesBlocked(_conservativeVars, _Nmailles, indices, initialFieldCons, INSERT_VALUES);
	VecAssemblyBegin(_conservativeVars);
	VecAssemblyEnd(_conservativeVars);
	VecCopy(_conservativeVars, _old);
	VecAssemblyBegin(_old);
	VecAssemblyEnd(_old);
	VecSetValuesBlocked(_primitiveVars, _Nmailles, indices, initialFieldPrim, INSERT_VALUES);
	VecAssemblyBegin(_primitiveVars);
	VecAssemblyEnd(_primitiveVars);
	if(_system)
	{
		cout << "Variables primitives initiales : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
		cout<<"Variables conservatives initiales : "<<endl;
		VecView(_conservativeVars,  PETSC_VIEWER_STDOUT_SELF);
	}

	delete[] initialFieldPrim;
	delete[] initialFieldCons;
	delete[] indices;

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

bool ProblemFluid::initTimeStep(double dt){
	_dt = dt;
	return _dt>0;//No need to call MatShift as the linear system matrix is filled at each Newton iteration (unlike linear problem)
}

bool ProblemFluid::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		_maxvp=0.;
		computeTimeStep(stop);//This compute timestep is just to update the linear system. The time step was imposed befor starting the Newton iterations
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
	else{//Implicit scheme

		KSPGetIterationNumber(_ksp, &_PetscIts);
		if( _MaxIterLinearSolver < _PetscIts)//save the maximum number of iterations needed during the newton scheme
			_MaxIterLinearSolver = _PetscIts;

		KSPConvergedReason reason;
		KSPGetConvergedReason(_ksp,&reason);

		if(reason<0)//solving the linear system failed
		{
			cout<<"Systeme lineaire : pas de convergence de Petsc. Raison PETSC numéro "<<reason<<endl;
			cout<<"Nombre d'itérations effectuées "<< _PetscIts<<" nombre maximal Itérations autorisées "<<_maxPetscIts<<endl;
			*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Raison PETSC numéro "<<reason<<endl;
			*_runLogFile<<"Nombre d'itérations effectuées "<< _PetscIts<<" nombre maximal Itérations autorisées "<<_maxPetscIts<<endl;
			converged=false;
			return false;
		}
		else{//solving the linear system succeeded
			//Calcul de la variation relative Uk+1-Uk
			_erreur_rel = 0.;
			double x, dx;
			int I;
			for(int j=1; j<=_Nmailles; j++)
			{
				for(int k=0; k<_nVar; k++)
				{
					I = (j-1)*_nVar + k;
					VecGetValues(_newtonVariation, 1, &I, &dx);
					VecGetValues(_conservativeVars, 1, &I, &x);
					if (fabs(x)*fabs(x)< _precision)
					{
						if(_erreur_rel < fabs(dx))
							_erreur_rel = fabs(dx);
					}
					else if(_erreur_rel < fabs(dx/x))
						_erreur_rel = fabs(dx/x);
				}
			}
		}
		converged = _erreur_rel <= _precision_Newton;
	}

	double relaxation=1;//Uk+1=Uk+relaxation*daltaU

	VecAXPY(_conservativeVars,  relaxation, _newtonVariation);

	//mise a jour du champ primitif
	updatePrimitives();

	if(_nbPhases==2 && fabs(_err_press_max) > _precision)//la pression n'a pu être calculée en diphasique à partir des variables conservatives
	{
		cout<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
		*_runLogFile<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
		converged=false;
		return false;
	}
	if(_system)
	{
		cout<<"Vecteur Ukp1-Uk "<<endl;
		VecView(_newtonVariation,  PETSC_VIEWER_STDOUT_SELF);
		cout << "Nouvel etat courant Uk de l'iteration Newton: " << endl;
		VecView(_conservativeVars,  PETSC_VIEWER_STDOUT_SELF);
	}

	if(_nbPhases==2 && _nbTimeStep%_freqSave ==0){
		if(_minm1<-_precision || _minm2<-_precision)
		{
			cout<<"!!!!!!!!! WARNING masse partielle negative sur " << _nbMaillesNeg << " faces, min m1= "<< _minm1 << " , minm2= "<< _minm2<< " precision "<<_precision<<endl;
			*_runLogFile<<"!!!!!!!!! WARNING masse partielle negative sur " << _nbMaillesNeg << " faces, min m1= "<< _minm1 << " , minm2= "<< _minm2<< " precision "<<_precision<<endl;
		}

		if (_nbVpCplx>0){
			cout << "!!!!!!!!!!!!!!!!!!!!!!!! Complex eigenvalues on " << _nbVpCplx << " cells, max imag= " << _part_imag_max << endl;
			*_runLogFile << "!!!!!!!!!!!!!!!!!!!!!!!! Complex eigenvalues on " << _nbVpCplx << " cells, max imag= " << _part_imag_max << endl;
		}
	}
	_minm1=1e30;
	_minm2=1e30;
	_nbMaillesNeg=0;
	_nbVpCplx =0;
	_part_imag_max=0;

	return true;
}

double ProblemFluid::computeTimeStep(bool & stop){
	VecZeroEntries(_b);

	if(_restartWithNewTimeScheme)//This is a change of time scheme during a simulation
	{
		if(_timeScheme == Implicit)
			MatCreateSeqBAIJ(PETSC_COMM_SELF, _nVar, _nVar*_Nmailles, _nVar*_Nmailles, (1+_neibMaxNb), PETSC_NULL, &_A);			
		else
			MatDestroy(&_A);
		_restartWithNewTimeScheme=false;
	}
	if(_timeScheme == Implicit)
		MatZeroEntries(_A);

	VecAssemblyBegin(_b);
	VecAssemblyBegin(_conservativeVars);
	std::vector< int > idCells;
	PetscInt idm, idn, size = 1;

	long nbFaces = _mesh.getNumberOfFaces();
	Face Fj;
	Cell Ctemp1,Ctemp2;
	string nameOfGroup;

	for (int j=0; j<nbFaces;j++){
		Fj = _mesh.getFace(j);
		_isBoundary=Fj.isBorder();
		idCells = Fj.getCellsId();

		// If Fj is on the boundary
		if (_isBoundary)
		{
			for(int k=0;k<Fj.getNumberOfCells();k++)//there will be at most two neighours in the case of an inner wall
			{
				// compute the normal vector corresponding to face j : from Ctemp1 outward
				Ctemp1 = _mesh.getCell(idCells[k]);//origin of the normal vector
				if (_Ndim >1){
					for(int l=0; l<Ctemp1.getNumberOfFaces(); l++)
					{//we look for l the index of the face Fj for the cell Ctemp1
						if (j == Ctemp1.getFacesId()[l])
						{
							for (int idim = 0; idim < _Ndim; ++idim)
								_vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
							break;
						}
					}
				}else{ // _Ndim = 1, build normal vector (bug cdmath)
					if(!_sectionFieldSet)
					{
						if (Fj.x()<Ctemp1.x())
							_vec_normal[0] = -1;
						else
							_vec_normal[0] = 1;
					}
					else
					{
						if(idCells[0]==0)
							_vec_normal[0] = -1;
						else//idCells[0]==31
							_vec_normal[0] = 1;
					}
				}
				if(_verbose && _nbTimeStep%_freqSave ==0)
				{
					cout << "face numero " << j << " cellule frontiere " << idCells[k] << " ; vecteur normal=(";
					for(int p=0; p<_Ndim; p++)
						cout << _vec_normal[p] << ",";
					cout << "). "<<endl;
				}
				nameOfGroup = Fj.getGroupName();
				_porosityi=_porosityField(idCells[k]);
				_porosityj=_porosityi;
				setBoundaryState(nameOfGroup,idCells[k],_vec_normal);
				convectionState(idCells[k],0,true);
				convectionMatrices();
				diffusionStateAndMatrices(idCells[k], 0, true);
				// compute 1/dxi
				if (_Ndim > 1)
					_inv_dxi = Fj.getMeasure()/Ctemp1.getMeasure();
				else
					_inv_dxi = 1/Ctemp1.getMeasure();

				addConvectionToSecondMember(idCells[k],-1,true,nameOfGroup);
				addDiffusionToSecondMember(idCells[k],-1,true);
				addSourceTermToSecondMember(idCells[k],(_mesh.getCell(idCells[k])).getNumberOfFaces(),-1, -1,true,j,_inv_dxi*Ctemp1.getMeasure());

				if(_timeScheme == Implicit){
					for(int l=0; l<_nVar*_nVar;l++){
						_AroeMinusImplicit[l] *= _inv_dxi;
						_Diffusion[l] *=_inv_dxi*_inv_dxi;
					}

					jacobian(idCells[k],nameOfGroup,_vec_normal);
					jacobianDiff(idCells[k],nameOfGroup);
					if(_verbose && _nbTimeStep%_freqSave ==0){
						cout << "Matrice Jacobienne CL convection:" << endl;
						for(int p=0; p<_nVar; p++){
							for(int q=0; q<_nVar; q++)
								cout << _Jcb[p*_nVar+q] << "\t";
							cout << endl;
						}
						cout << endl;
						cout << "Matrice Jacobienne CL diffusion:" << endl;
						for(int p=0; p<_nVar; p++){
							for(int q=0; q<_nVar; q++)
								cout << _JcbDiff[p*_nVar+q] << "\t";
							cout << endl;
						}
						cout << endl;
					}
					idm = idCells[k];
					Polynoms Poly;
					//calcul et insertion de A^-*Jcb
					Poly.matrixProduct(_AroeMinusImplicit, _nVar, _nVar, _Jcb, _nVar, _nVar, _a);
					MatSetValuesBlocked(_A, size, &idm, size, &idm, _a, ADD_VALUES);

					if(_verbose)
						displayMatrix(_a, _nVar, "produit A^-*Jcb pour CL");

					//insertion de -A^-
					for(int k=0; k<_nVar*_nVar;k++){
						_AroeMinusImplicit[k] *= -1;
					}
					MatSetValuesBlocked(_A, size, &idm, size, &idm, _AroeMinusImplicit, ADD_VALUES);
					if(_verbose)
						displayMatrix(_AroeMinusImplicit, _nVar,"-_AroeMinusImplicit: ");

					//calcul et insertion de D*JcbDiff
					Poly.matrixProduct(_Diffusion, _nVar, _nVar, _JcbDiff, _nVar, _nVar, _a);
					MatSetValuesBlocked(_A, size, &idm, size, &idm, _a, ADD_VALUES);
					for(int k=0; k<_nVar*_nVar;k++)
						_Diffusion[k] *= -1;
					MatSetValuesBlocked(_A, size, &idm, size, &idm, _Diffusion, ADD_VALUES);
				}
			}
		} else 	if (Fj.getNumberOfCells()==2 ){	// Fj is inside the domain and has two neighours (no junction)
			// compute the normal vector corresponding to face j : from Ctemp1 to Ctemp2
			Ctemp1 = _mesh.getCell(idCells[0]);//origin of the normal vector
			Ctemp2 = _mesh.getCell(idCells[1]);
			if (_Ndim >1){
				for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
					if (j == Ctemp1.getFacesId()[l]){
						for (int idim = 0; idim < _Ndim; ++idim)
							_vec_normal[idim] = Ctemp1.getNormalVector(l,idim);
						break;
					}
				}
			}else{ // _Ndim = 1, build normal vector (bug cdmath)
				if(!_sectionFieldSet)
				{
					if (Fj.x()<Ctemp1.x())
						_vec_normal[0] = -1;
					else
						_vec_normal[0] = 1;
				}
				else
				{
					if(idCells[0]>idCells[1])
						_vec_normal[0] = -1;
					else
						_vec_normal[0] = 1;
				}
			}
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "face numero " << j << " cellule gauche " << idCells[0] << " cellule droite " << idCells[1];
				cout<<" Normal vector= ";
				for (int idim = 0; idim < _Ndim; ++idim)
					cout<<_vec_normal[idim]<<", ";
				cout<<endl;
			}
			_porosityi=_porosityField(idCells[0]);
			_porosityj=_porosityField(idCells[1]);
			convectionState(idCells[0],idCells[1],false);
			convectionMatrices();
			diffusionStateAndMatrices(idCells[0], idCells[1], false);

			// compute 1/dxi and 1/dxj
			if (_Ndim > 1)
			{
				_inv_dxi = Fj.getMeasure()/Ctemp1.getMeasure();
				_inv_dxj = Fj.getMeasure()/Ctemp2.getMeasure();
			}
			else
			{
				_inv_dxi = 1/Ctemp1.getMeasure();
				_inv_dxj = 1/Ctemp2.getMeasure();
			}

			addConvectionToSecondMember( idCells[0],idCells[1], false);
			addDiffusionToSecondMember( idCells[0],idCells[1], false);
			addSourceTermToSecondMember(idCells[0], Ctemp1.getNumberOfFaces(),idCells[1], _mesh.getCell(idCells[1]).getNumberOfFaces(),false,j,_inv_dxi*Ctemp1.getMeasure());

			if(_timeScheme == Implicit){
				for(int k=0; k<_nVar*_nVar;k++)
				{
					_AroeMinusImplicit[k] *= _inv_dxi;
					_Diffusion[k] *=_inv_dxi*2/(1/_inv_dxi+1/_inv_dxj);
				}
				idm = idCells[0];
				idn = idCells[1];
				//cout<<"idm= "<<idm<<"idn= "<<idn<<"nbvoismax= "<<_neibMaxNb<<endl;
				MatSetValuesBlocked(_A, size, &idm, size, &idn, _AroeMinusImplicit, ADD_VALUES);
				MatSetValuesBlocked(_A, size, &idm, size, &idn, _Diffusion, ADD_VALUES);

				if(_verbose){
					displayMatrix(_AroeMinusImplicit, _nVar, "+_AroeMinusImplicit: ");
					displayMatrix(_Diffusion, _nVar, "+_Diffusion: ");
				}
				for(int k=0;k<_nVar*_nVar;k++){
					_AroeMinusImplicit[k] *= -1;
					_Diffusion[k] *= -1;
				}
				MatSetValuesBlocked(_A, size, &idm, size, &idm, _AroeMinusImplicit, ADD_VALUES);
				MatSetValuesBlocked(_A, size, &idm, size, &idm, _Diffusion, ADD_VALUES);
				if(_verbose){
					displayMatrix(_AroeMinusImplicit, _nVar, "-_AroeMinusImplicit: ");
					displayMatrix(_Diffusion, _nVar, "-_Diffusion: ");
				}
				for(int k=0; k<_nVar*_nVar;k++)
				{
					_AroePlusImplicit[k]  *= _inv_dxj;
					_Diffusion[k] *=_inv_dxj/_inv_dxi;
				}
				MatSetValuesBlocked(_A, size, &idn, size, &idn, _AroePlusImplicit, ADD_VALUES);
				MatSetValuesBlocked(_A, size, &idn, size, &idn, _Diffusion, ADD_VALUES);
				if(_verbose)
					displayMatrix(_AroePlusImplicit, _nVar, "+_AroePlusImplicit: ");

				for(int k=0;k<_nVar*_nVar;k++){
					_AroePlusImplicit[k] *= -1;
					_Diffusion[k] *= -1;
				}
				MatSetValuesBlocked(_A, size, &idn, size, &idm, _AroePlusImplicit, ADD_VALUES);
				MatSetValuesBlocked(_A, size, &idn, size, &idm, _Diffusion, ADD_VALUES);

				if(_verbose)
					displayMatrix(_AroePlusImplicit, _nVar, "-_AroePlusImplicit: ");
			}
		}
		else if( Fj.getNumberOfCells()>2 && _Ndim==1 ){//inner face with more than two neighbours
			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"lattice mesh junction at face "<<j<<" nbvoismax= "<<_neibMaxNb<<endl;
			*_runLogFile<<"Warning: treatment of a junction node"<<endl;

			if(!_sectionFieldSet)
			{
				_runLogFile->close();
				throw CdmathException("ProblemFluid::ComputeTimeStep(): pipe network requires section field");
			}
			int largestSectionCellIndex=0;
			for(int i=1;i<Fj.getNumberOfCells();i++){
				if(_sectionField(idCells[i])>_sectionField(idCells[largestSectionCellIndex]))
					largestSectionCellIndex=i;
			}
			idm = idCells[largestSectionCellIndex];
			Ctemp1 = _mesh.getCell(idm);//origin of the normal vector
			_porosityi=_porosityField(idm);

			if (j==15)// bug cdmath (Fj.x() > _mesh.getCell(idm).x())
				_vec_normal[0] = 1;
			else//j==16
				_vec_normal[0] = -1;
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout<<"Cell with largest section has index "<< largestSectionCellIndex <<" and number "<<idm<<endl;
				cout << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << _vec_normal[p] << ",";
				cout << "). "<<endl;
			}
			for(int i=0;i<Fj.getNumberOfCells();i++){
				if(i != largestSectionCellIndex){
					idn = idCells[i];
					Ctemp2 = _mesh.getCell(idn);
					_porosityj=_porosityField(idn);
					convectionState(idm,idn,false);
					convectionMatrices();
					diffusionStateAndMatrices(idm, idn,false);

					if(_verbose && _nbTimeStep%_freqSave ==0)
						cout<<"Neighbour index "<<i<<" cell number "<< idn<<endl;

					_inv_dxi = _sectionField(idn)/_sectionField(idm)/Ctemp1.getMeasure();
					_inv_dxj = 1/Ctemp2.getMeasure();

					addConvectionToSecondMember(idm,idn, false);
					_inv_dxi = sqrt(_sectionField(idn)/_sectionField(idm))/Ctemp1.getMeasure();
					addDiffusionToSecondMember(idm,idn, false);
					_inv_dxi = _sectionField(idn)/_sectionField(idm)/Ctemp1.getMeasure();
					addSourceTermToSecondMember(idm, Ctemp1.getNumberOfFaces()*(Fj.getNumberOfCells()-1),idn, Ctemp2.getNumberOfFaces(),false,j,_inv_dxi*Ctemp1.getMeasure());

					if(_timeScheme == Implicit){
						for(int k=0; k<_nVar*_nVar;k++)
						{
							_AroeMinusImplicit[k] *= _inv_dxi;
							_Diffusion[k] *=_inv_dxi*2/(1/_inv_dxi+1/_inv_dxj);//use sqrt as above
						}
						MatSetValuesBlocked(_A, size, &idm, size, &idn, _AroeMinusImplicit, ADD_VALUES);
						MatSetValuesBlocked(_A, size, &idm, size, &idn, _Diffusion, ADD_VALUES);

						if(_verbose){
							displayMatrix(_AroeMinusImplicit, _nVar, "+_AroeMinusImplicit: ");
							displayMatrix(_Diffusion, _nVar, "+_Diffusion: ");
						}
						for(int k=0;k<_nVar*_nVar;k++){
							_AroeMinusImplicit[k] *= -1;
							_Diffusion[k] *= -1;
						}
						MatSetValuesBlocked(_A, size, &idm, size, &idm, _AroeMinusImplicit, ADD_VALUES);
						MatSetValuesBlocked(_A, size, &idm, size, &idm, _Diffusion, ADD_VALUES);
						if(_verbose){
							displayMatrix(_AroeMinusImplicit, _nVar, "-_AroeMinusImplicit: ");
							displayMatrix(_Diffusion, _nVar, "-_Diffusion: ");
						}
						for(int k=0; k<_nVar*_nVar;k++)
						{
							_AroePlusImplicit[k] *= _inv_dxj;
							_Diffusion[k] *=_inv_dxj/_inv_dxi;//use sqrt as above
						}
						MatSetValuesBlocked(_A, size, &idn, size, &idn, _AroePlusImplicit, ADD_VALUES);
						MatSetValuesBlocked(_A, size, &idn, size, &idn, _Diffusion, ADD_VALUES);
						if(_verbose)
							displayMatrix(_AroePlusImplicit, _nVar, "+_AroePlusImplicit: ");

						for(int k=0;k<_nVar*_nVar;k++){
							_AroePlusImplicit[k] *= -1;
							_Diffusion[k] *= -1;
						}
						MatSetValuesBlocked(_A, size, &idn, size, &idm, _AroePlusImplicit, ADD_VALUES);
						MatSetValuesBlocked(_A, size, &idn, size, &idm, _Diffusion, ADD_VALUES);

						if(_verbose)
							displayMatrix(_AroePlusImplicit, _nVar, "-_AroePlusImplicit: ");
					}
				}
			}
		}
		else
		{
			cout<< "Face j="<<j<< " is not a boundary face and has "<<Fj.getNumberOfCells()<< " neighbour cells"<<endl;
			_runLogFile->close();
			throw CdmathException("ProblemFluid::ComputeTimeStep(): incompatible number of cells around a face");
		}

	}
	VecAssemblyEnd(_conservativeVars);
	VecAssemblyEnd(_b);

	if(_timeScheme == Implicit){
		for(int imaille = 0; imaille<_Nmailles; imaille++)
			MatSetValuesBlocked(_A, size, &imaille, size, &imaille, _GravityImplicitationMatrix, ADD_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
			displayMatrix(_GravityImplicitationMatrix,_nVar,"Gravity matrix:");

		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
		if(_verbose && _nbTimeStep%_freqSave ==0){
			cout << "Fin ProblemFluid.cxx : matrice implicite"<<endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
	}

	stop=false;


	/*
	if(_nbTimeStep+1<_cfl)
		return (_nbTimeStep+1)*_minl/_maxvp;
	else
	 */
	return _cfl*_minl/_maxvp;
}

void ProblemFluid::computeNewtonVariation()
{
	if(_verbose)
	{
		cout<<"Vecteur courant Uk "<<endl;
		VecView(_conservativeVars,PETSC_VIEWER_STDOUT_SELF);
		cout << endl;
	}
	if(_timeScheme == Explicit)
	{
		VecCopy(_b,_newtonVariation);
		VecScale(_newtonVariation, _dt);
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"Vecteur _newtonVariation =_b*dt"<<endl;
			VecView(_newtonVariation,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
	}
	else
	{
		if(_verbose)
		{
			cout << "Matrice du système linéaire avant contribution delta t" << endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Second membre du système linéaire avant contribution delta t" << endl;
			VecView(_b, PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

		VecAXPY(_b, 1/_dt, _old);
		VecAXPY(_b, -1/_dt, _conservativeVars);
		MatShift(_A, 1/_dt);

#if PETSC_VERSION_GREATER_3_5
		KSPSetOperators(_ksp, _A, _A);
#else
		KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

		if(_verbose)
		{
			cout << "Matrice du système linéaire après contribution delta t" << endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Second membre du système linéaire après contribution delta t" << endl;
			VecView(_b, PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}

		if(_conditionNumber)
			KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
		if(!_isScaling)
		{
			KSPSolve(_ksp, _b, _newtonVariation);
		}
		else
		{
			computeScaling(_maxvp);
			int indice;
			VecAssemblyBegin(_vecScaling);
			VecAssemblyBegin(_invVecScaling);
			for(int imaille = 0; imaille<_Nmailles; imaille++)
			{
				indice = imaille;
				VecSetValuesBlocked(_vecScaling,1 , &indice, _blockDiag, INSERT_VALUES);
				VecSetValuesBlocked(_invVecScaling,1,&indice,_invBlockDiag, INSERT_VALUES);
			}
			VecAssemblyEnd(_vecScaling);
			VecAssemblyEnd(_invVecScaling);

			if(_system)
			{
				cout << "Matrice avant le preconditionneur des vecteurs propres " << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
				cout << "Second membre avant le preconditionneur des vecteurs propres " << endl;
				VecView(_b, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
			MatDiagonalScale(_A,_vecScaling,_invVecScaling);
			if(_system)
			{
				cout << "Matrice apres le preconditionneur des vecteurs propres " << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
			VecCopy(_b,_bScaling);
			VecPointwiseMult(_b,_vecScaling,_bScaling);
			if(_system)
			{
				cout << "Produit du second membre par le preconditionneur bloc diagonal  a gauche" << endl;
				VecView(_b, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}

			KSPSolve(_ksp,_b, _bScaling);
			VecPointwiseMult(_newtonVariation,_invVecScaling,_bScaling);
		}
		if(_verbose)
		{
			cout << "solution du systeme lineaire local:" << endl;
			VecView(_newtonVariation, PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
	}
}

void ProblemFluid::validateTimeStep()
{
	if(_system)
	{
		cout<<" Vecteur Un"<<endl;
		VecView(_old,  PETSC_VIEWER_STDOUT_WORLD);
		cout<<" Vecteur Un+1"<<endl;
		VecView(_conservativeVars,  PETSC_VIEWER_STDOUT_WORLD);
	}
	VecAXPY(_old,  -1, _conservativeVars);//old contient old-courant

	//Calcul de la variation Un+1-Un
	_erreur_rel= 0;
	double x, dx;
	int I;

	for(int j=1; j<=_Nmailles; j++)
	{
		for(int k=0; k<_nVar; k++)
		{
			I = (j-1)*_nVar + k;
			VecGetValues(_old, 1, &I, &dx);
			VecGetValues(_conservativeVars, 1, &I, &x);
			if (fabs(x)< _precision)
			{
				if(_erreur_rel < fabs(dx))
					_erreur_rel = fabs(dx);
			}
			else if(_erreur_rel < fabs(dx/x))
				_erreur_rel = fabs(dx/x);
		}
	}

	_isStationary =_erreur_rel <_precision;

	VecCopy(_conservativeVars, _old);

	if(_verbose && _nbTimeStep%_freqSave ==0){
		if(!_usePrimitiveVarsInNewton)
			testConservation();
		cout <<"Valeur propre locale max: " << _maxvp << endl;
	}

	if(_nbPhases==2 && _nbTimeStep%_freqSave ==0){
		//Find minimum and maximum void fractions
		double alphamin=1e30;
		double alphamax=-1e30;
		double T, Tmax=-1e30;
		int J=_nVar-1;
		I = 0;
		for(int j=0; j<_Nmailles; j++)
		{
			VecGetValues(_primitiveVars, 1, &I, &x);//extract void fraction
			if(x>alphamax)
				alphamax=x;
			if(x<alphamin)
				alphamin=x;
			I += _nVar;
			VecGetValues(_primitiveVars, 1, &J, &T);//extract void fraction
			if(T>Tmax)
				Tmax=T;
			J += _nVar;
		}
		cout<<"Alpha min = " << alphamin << ", Alpha max = " << alphamax<<", Tmax= "<<Tmax<<" K"<<endl;
		cout<<endl;
		*_runLogFile<<"Alpha min = " << alphamin << " Alpha max = " << alphamax<<", Tmax= "<<Tmax<<" K"<<endl;
	}

	_time+=_dt;
	_nbTimeStep++;
	if (_nbTimeStep%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
		save();
}

void ProblemFluid::abortTimeStep(){
	_dt = 0;
}

void ProblemFluid::addConvectionToSecondMember
(		const int &i,
		const int &j, bool isBord, string groupname
)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"ProblemFluid::addConvectionToSecondMember start"<<endl;

	//extraction des valeurs
	for(int k=0; k<_nVar; k++)
		_idm[k] = _nVar*i + k;
	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

	if(!isBord){
		for(int k=0; k<_nVar; k++)
			_idn[k] = _nVar*j + k;
		VecGetValues(_conservativeVars, _nVar, _idn, _Uj);
		VecGetValues(_primitiveVars, _nVar, _idn, _Vj);
	}
	else{
		for(int k=0; k<_nVar; k++)
			_idn[k] = k;
		VecGetValues(_Uext, _nVar, _idn, _Uj);
		consToPrim(_Uj, _Vj,_porosityj);
	}
	_idm[0] = i;
	_idn[0] = j;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "addConvectionToSecondMember : état i= " << i << ", _Ui=" << endl;
		for(int q=0; q<_nVar; q++)
			cout << _Ui[q] << endl;
		cout << endl;
		cout << "addConvectionToSecondMember : état j= " << j << ", _Uj=" << endl;
		for(int q=0; q<_nVar; q++)
			cout << _Uj[q] << endl;
		cout << endl;
	}
	if(_nonLinearFormulation!=reducedRoe){//VFRoe, Roe or VFFC
		Vector Ui(_nVar), Uj(_nVar), Vi(_nVar), Vj(_nVar), Fi(_nVar), Fj(_nVar),Fij(_nVar);
		for(int i1=0;i1<_nVar;i1++)
		{
			Ui(i1)=_Ui[i1];
			Uj(i1)=_Uj[i1];
			Vi(i1)=_Vi[i1];
			Vj(i1)=_Vj[i1];
		}
		Vector normale(_Ndim);
		for(int i1=0;i1<_Ndim;i1++)
			normale(i1)=_vec_normal[i1];

		Matrix signAroe(_nVar);
		for(int i1=0;i1<_nVar;i1++)
			for(int i2=0;i2<_nVar;i2++)
				signAroe(i1,i2)=_signAroe[i1*_nVar+i2];

		Matrix absAroe(_nVar);
		for(int i1=0;i1<_nVar;i1++)
			for(int i2=0;i2<_nVar;i2++)
				absAroe(i1,i2)=_absAroe[i1*_nVar+i2];

		if(_nonLinearFormulation==VFRoe)//VFRoe
		{
			Vector Uij(_nVar),Vij(_nVar);
			double porosityij=sqrt(_porosityi*_porosityj);

			Uij=(Ui+Uj)/2+signAroe*(Ui-Uj)/2;

			for(int i1=0;i1<_nVar;i1++)
				_Uij[i1]=Uij(i1);

			consToPrim(_Uij, _Vij,porosityij);

			applyVFRoeLowMachCorrections(isBord, groupname);

			for(int i1=0;i1<_nVar;i1++)
			{
				Uij(i1)=_Uij[i1];
				Vij(i1)=_Vij[i1];
			}

			Fij=convectionFlux(Uij,Vij,normale,porosityij);

			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout<<"Etat interfacial conservatif "<<i<<", "<<j<< endl;
				cout<<Uij<<endl;
				cout<<"Etat interfacial primitif "<<i<<", "<<j<< endl;
				cout<<Vij<<endl;
			}
		}
		else //Roe or VFFC
		{
			if(_spaceScheme==staggered && _nonLinearFormulation==VFFC)//special case
				Fij=staggeredVFFCFlux();
			else
			{
				Fi=convectionFlux(Ui,Vi,normale,_porosityi);
				Fj=convectionFlux(Uj,Vj,normale,_porosityj);
				if(_nonLinearFormulation==VFFC)//VFFC
					Fij=(Fi+Fj)/2+signAroe*(Fi-Fj)/2;
				else if(_nonLinearFormulation==Roe)//Roe
					Fij=(Fi+Fj)/2+absAroe*(Ui-Uj)/2;

				if(_verbose && _nbTimeStep%_freqSave ==0)
				{
					cout<<"Flux cellule "<<i<<" = "<< endl;
					cout<<Fi<<endl;
					cout<<"Flux cellule "<<j<<" = "<< endl;
					cout<<Fj<<endl;
				}
			}
		}
		for(int i1=0;i1<_nVar;i1++)
			_phi[i1]=-Fij(i1)*_inv_dxi;
		VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Ajout convection au 2nd membre pour les etats " << i << "," << j << endl;
			cout<<"Flux interfacial "<<i<<", "<<j<< endl;
			cout<<Fij<<endl;
			cout << "Contribution convection à " << i << ", -Fij(i1)*_inv_dxi= "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}
		if(!isBord){
			for(int i1=0;i1<_nVar;i1++)
				_phi[i1]*=-_inv_dxj/_inv_dxi;
			VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution convection à " << j << ", Fij(i1)*_inv_dxj= "<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _phi[q] << endl;
				cout << endl;
			}
		}
	}else //_nonLinearFormulation==reducedRoe)
	{
		for(int k=0; k<_nVar; k++)
			_temp[k]=(_Ui[k] - _Uj[k])*_inv_dxi;//(Ui-Uj)*_inv_dxi
		Polynoms Poly;
		Poly.matrixProdVec(_AroeMinus, _nVar, _nVar, _temp, _phi);//phi=A^-(U_i-U_j)/dx
		VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Ajout convection au 2nd membre pour les etats " << i << "," << j << endl;
			cout << "(Ui - Uj)*_inv_dxi= "<<endl;;
			for(int q=0; q<_nVar; q++)
				cout << _temp[q] << endl;
			cout << endl;
			cout << "Contribution convection à " << i << ", A^-*(Ui - Uj)*_inv_dxi= "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}

		if(!isBord)
		{
			for(int k=0; k<_nVar; k++)
				_temp[k]*=_inv_dxj/_inv_dxi;//(Ui-Uj)*_inv_dxj
			Poly.matrixProdVec(_AroePlus, _nVar, _nVar, _temp, _phi);//phi=A^+(U_i-U_j)/dx
			VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);

			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution convection à  " << j << ", A^+*(Ui - Uj)*_inv_dxi= "<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _phi[q] << endl;
				cout << endl;
			}
		}
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"ProblemFluid::addConvectionToSecondMember end : matrices de décentrement cellules i= " << i << ", et j= " << j<< "):"<<endl;
		displayMatrix(_absAroe, _nVar,"Valeur absolue matrice de Roe");
		displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus");
		displayMatrix(_AroePlus, _nVar,"Matrice _AroePlus");
		displayMatrix(_signAroe, _nVar,"Signe de la matrice de Roe");
	}
}

void ProblemFluid::addSourceTermToSecondMember
(	const int i, int nbVoisinsi,
		const int j, int nbVoisinsj,
		bool isBord, int ij, double mesureFace)//To do : generalise to unstructured meshes
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"ProblemFluid::addSourceTerm cell i= "<<i<< " cell j= "<< j<< " isbord "<<isBord<<endl;

	_idm[0] = i*_nVar;
	for(int k=1; k<_nVar;k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	sourceVector(_Si,_Ui,_Vi,i);

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Terme source S(Ui), i= " << i<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Si[q] << endl;
		cout << endl;
	}
	if(!isBord){
		for(int k=0; k<_nVar; k++)
			_idn[k] = _nVar*j + k;
		VecGetValues(_conservativeVars, _nVar, _idn, _Uj);
		VecGetValues(_primitiveVars, _nVar, _idn, _Vj);
		sourceVector(_Sj,_Uj,_Vj,j);
	}else
	{
		for(int k=0; k<_nVar; k++)
			_idn[k] = k;
		VecGetValues(_Uext, _nVar, _idn, _Uj);
		consToPrim(_Uj, _Vj,_porosityj);
		sourceVector(_Sj,_Uj,_Vj,i);
	}
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		if(!isBord)
			cout << "Terme source S(Uj), j= " << j<<endl;
		else
			cout << "Terme source S(Uj), cellule fantôme "<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Sj[q] << endl;
		cout << endl;
	}

	//Compute pressure loss vector
	double K;
	if(_pressureLossFieldSet){
		K=_pressureLossField(ij);
		pressureLossVector(_pressureLossVector, K,_Ui, _Vi, _Uj, _Vj);	
	}
	else{
		K=0;
		for(int k=0; k<_nVar;k++)
			_pressureLossVector[k]=0;
	}
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"interface i= "<<i<<", j= "<<j<<", ij="<<ij<<", K="<<K<<endl;
		for(int k=0; k<_nVar;k++)
			cout<< _pressureLossVector[k]<<", ";
		cout<<endl;
	}
	//Contribution of the porosityField gradient:
	if(!isBord)
		porosityGradientSourceVector();
	else{
		for(int k=0; k<_nVar;k++)
			_porosityGradientSourceVector[k]=0;
	}

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		if(!isBord)
			cout<<"interface i= "<<i<<", j= "<<j<<", ij="<<ij<<", dxi= "<<1/_inv_dxi<<", dxj= "<<1/_inv_dxj<<endl;
		else
			cout<<"interface frontière i= "<<i<<", ij="<<ij<<", dxi= "<<1/_inv_dxi<<endl;
		cout<<"Gradient de porosite à l'interface"<<endl;
		for(int k=0; k<_nVar;k++)
			cout<< _porosityGradientSourceVector[k]<<", ";
		cout<<endl;
	}
	Polynoms Poly;
	if(!isBord){
		if(_wellBalancedCorrection){
			for(int k=0; k<_nVar;k++)
				_phi[k]=(_Si[k]+_Sj[k])/2+_pressureLossVector[k]+_porosityGradientSourceVector[k];
			Poly.matrixProdVec(_signAroe, _nVar, _nVar, _phi, _l);
			for(int k=0; k<_nVar;k++){
				_Si[k]=(_phi[k]-_l[k])*mesureFace/_perimeters(i);///nbVoisinsi;
				_Sj[k]=(_phi[k]+_l[k])*mesureFace/_perimeters(j);///nbVoisinsj;
			}
			if (_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution au terme source Si de la cellule i= " << i<<" venant  (après décentrement) de la face (i,j), j="<<j<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Si[q] << endl;
				cout << "Contribution au terme source Sj de la cellule j= " << j<<" venant  (après décentrement) de la face (i,j), i="<<i<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Sj[q] << endl;
				cout << endl;
				cout<<"ratio surface sur volume i = "<<mesureFace/_perimeters(i)<<" perimeter = "<< _perimeters(i)<<endl;
				cout<<"ratio surface sur volume j = "<<mesureFace/_perimeters(j)<<" perimeter = "<< _perimeters(j)<<endl;
				cout << endl;
			}
		}
		else{
			for(int k=0; k<_nVar;k++){
				_Si[k]=_Si[k]/nbVoisinsi+_pressureLossVector[k]/2+_porosityGradientSourceVector[k]/2;//mesureFace/_perimeters(i)
				_Sj[k]=_Sj[k]/nbVoisinsj+_pressureLossVector[k]/2+_porosityGradientSourceVector[k]/2;//mesureFace/_perimeters(j)
			}
			if (_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution au terme source Si de la cellule i = " << i<<" venant  de la face (i,j), j="<<j<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Si[q] << endl;
				cout << "Contribution au terme source Sj de la cellule j = " << j<<" venant  de la face (i,j), i="<<i <<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Sj[q] << endl;
				cout << endl;
			}
		}
		_idn[0] = j;
		VecSetValuesBlocked(_b, 1, _idn, _Sj, ADD_VALUES);
	}else{
		if(_wellBalancedCorrection){
			for(int k=0; k<_nVar;k++)
				_phi[k]=(_Si[k]+_Sj[k])/2+_pressureLossVector[k]+_porosityGradientSourceVector[k];
			Poly.matrixProdVec(_signAroe, _nVar, _nVar, _phi, _l);
			for(int k=0; k<_nVar;k++)
				_Si[k]=(_phi[k]-_l[k])*mesureFace/_perimeters(i);///nbVoisinsi;
			if (_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution au terme source Si de la cellule i= " << i<<" venant  (après décentrement) de la face (i,bord)"<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Si[q] << endl;
				cout<<"ratio surface sur volume i ="<<mesureFace/_perimeters(i)<<" perimeter = "<< _perimeters(i)<<endl;
				cout << endl;
			}
		}
		else
		{
			for(int k=0; k<_nVar;k++)
				_Si[k]=_Si[k]/nbVoisinsi+_pressureLossVector[k]/2+_porosityGradientSourceVector[k]/2;//mesureFace/_perimeters(i);//
			if (_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution au terme source Si de la cellule i = " << i<<" venant de la face (i,bord) "<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _Si[q] << endl;
				cout << endl;
			}
		}
	}
	_idm[0] = i;
	VecSetValuesBlocked(_b, 1, _idm, _Si, ADD_VALUES);

	if(_verbose && _nbTimeStep%_freqSave ==0 && _wellBalancedCorrection)
		displayMatrix( _signAroe,_nVar,"Signe matrice de Roe");
}

void ProblemFluid::updatePrimitives()
{
	VecAssemblyBegin(_primitiveVars);
	for(int i=1; i<=_Nmailles; i++)
	{
		_idm[0] = _nVar*( (i-1));
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "ProblemFluid::updatePrimitives() cell " << i-1 << endl;
			cout << "Ui = ";
			for(int q=0; q<_nVar; q++)
				cout << _Ui[q]  << "\t";
			cout << endl;
			cout << "Vi = ";
			for(int q=0; q<_nVar; q++)
				cout << _Vi[q]  << "\t";
			cout << endl;
		}

		consToPrim(_Ui,_Vi,_porosityField(i-1));
		if(_nbPhases==2 && _Psat>-1e30){//Cas simulation flashing
			double pressure;
			VecGetValues(_primitiveVars, 1, _idm+1, &pressure);
			_dp_over_dt(i-1)=(_Vi[1]-pressure)/_dt;//pn+1-pn
		}
		_idm[0] = i-1;
		VecSetValuesBlocked(_primitiveVars, 1, _idm, _Vi, INSERT_VALUES);
	}
	VecAssemblyEnd(_primitiveVars);

	if(_system)
	{
		cout << "Nouvelles variables primitives : " << endl;
		VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}
}

void ProblemFluid::updateConservatives()
{
	VecAssemblyBegin(_conservativeVars);
	for(int i=1; i<=_Nmailles; i++)
	{
		_idm[0] = _nVar*( (i-1));
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

		primToCons(_Vi,0,_Ui,0);
		_idm[0] = i-1;
		VecSetValuesBlocked(_conservativeVars, 1, _idm, _Ui, INSERT_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "ProblemFluid::updateConservatives() cell " << i-1 << endl;
			cout << "Vi = ";
			for(int q=0; q<_nVar; q++)
				cout << _Vi[q]  << "\t";
			cout << endl;
			cout << "Ui = ";
			for(int q=0; q<_nVar; q++)
				cout << _Ui[q]  << "\t";
			cout << endl;
		}
	}
	VecAssemblyEnd(_conservativeVars);

	if(_system)
	{
		cout << "Nouvelles variables conservatives : " << endl;
		VecView(_conservativeVars,  PETSC_VIEWER_STDOUT_WORLD);
		cout << endl;
	}
}

vector< complex<double> > ProblemFluid::getRacines(vector< double > pol_car){
	int degre_polynom = pol_car.size()-1;
	double tmp;
	vector< complex< double > > valeurs_propres (degre_polynom);
	vector< double > zeror(degre_polynom);
	vector< double > zeroi(degre_polynom);
	for(int j=0; j<(degre_polynom+1)/2; j++){//coefficients in order of decreasing powers for rpoly
		tmp=pol_car[j];
		pol_car[j] =pol_car[degre_polynom-j];
		pol_car[degre_polynom-j]=tmp;
	}

	//Calcul des racines du polynome
	roots_polynoms roots;
	int size_vp = roots.rpoly(&pol_car[0],degre_polynom,&zeror[0],&zeroi[0]);

	//On ordonne les valeurs propres
	if(zeror[1]<zeror[0])
	{
		tmp=zeror[0];
		zeror[0]=zeror[1];
		zeror[1]=tmp;
		tmp=zeroi[0];
		zeroi[0]=zeroi[1];
		zeroi[1]=tmp;
	}

	if(size_vp == degre_polynom) {
		for(int ct =0; ct<degre_polynom; ct++)
			valeurs_propres[ct] = complex<double> (zeror[ct],zeroi[ct]);  //vp non triviales
	}
	else {
		cout << " Pb, found only " << size_vp << " eigenvalues over " << degre_polynom << endl;
		*_runLogFile << " Pb, found only " << size_vp << " eigenvalues over " << degre_polynom << endl;
		cout<<"Coefficients polynome caracteristique: "<<endl;
		for(int ct =0; ct<degre_polynom+1; ct++)
			cout<<pol_car[ct]<< " , " <<endl;

		*_runLogFile<<"getRacines computation failed"<<endl;
		_runLogFile->close();
		throw CdmathException("getRacines computation failed");
	}

	return valeurs_propres;
}

void ProblemFluid::AbsMatriceRoe(vector< complex<double> > valeurs_propres_dist)
{
	Polynoms Poly;
	int nbVp_dist=valeurs_propres_dist.size();
	vector< complex< double > > y (nbVp_dist,0);
	for( int i=0 ; i<nbVp_dist ; i++)
		y[i] = Poly.abs_generalise(valeurs_propres_dist[i]);
	Poly.abs_par_interp_directe(nbVp_dist,valeurs_propres_dist, _Aroe, _nVar,_precision, _absAroe,y);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<< endl<<"ProblemFluid::AbsMatriceRoe: Valeurs propres :" << nbVp_dist<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<valeurs_propres_dist[ct].real()<< ", " <<valeurs_propres_dist[ct].imag() <<")  ";
		cout<< endl;
		cout<<"ProblemFluid::AbsMatriceRoe: Valeurs à atteindre pour l'interpolation"<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<y[ct].real()<< ", " <<y[ct].imag() <<")  ";
		cout<< endl;
		displayMatrix(_absAroe,_nVar,"Valeur absolue de la matrice de Roe");
	}
}

void ProblemFluid::SigneMatriceRoe(vector< complex<double> > valeurs_propres_dist)
{
	int nbVp_dist=valeurs_propres_dist.size();
	vector< complex< double > > y (nbVp_dist,0);
	for( int i=0 ; i<nbVp_dist ; i++)
	{
		if(valeurs_propres_dist[i].real()>0)
			y[i] = 1;
		else if(valeurs_propres_dist[i].real()<0)
			y[i] = -1;
		else
			y[i] = 0;
	}
	Polynoms Poly;
	Poly.abs_par_interp_directe(nbVp_dist,valeurs_propres_dist, _Aroe, _nVar,_precision, _signAroe,y);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<< endl<<"ProblemFluid::SigneMatriceRoe: Valeurs propres :" << nbVp_dist<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<valeurs_propres_dist[ct].real()<< ", " <<valeurs_propres_dist[ct].imag() <<")  ";
		cout<< endl;
		cout<<"ProblemFluid::SigneMatriceRoe: Valeurs à atteindre pour l'interpolation"<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<y[ct].real()<< ", " <<y[ct].imag() <<")  ";
		cout<< endl;
		displayMatrix(_signAroe,_nVar,"Signe matrice de Roe");
	}
}
void ProblemFluid::InvMatriceRoe(vector< complex<double> > valeurs_propres_dist)
{
	int nbVp_dist=valeurs_propres_dist.size();
	vector< complex< double > > y (nbVp_dist,0);
	Polynoms Poly;
	for( int i=0 ; i<nbVp_dist ; i++)
	{
		if(Poly.module(valeurs_propres_dist[i])>_precision)
			y[i] = 1./valeurs_propres_dist[i];
		else
			y[i] = 1./_precision;
	}
	Poly.abs_par_interp_directe(nbVp_dist,valeurs_propres_dist, _Aroe, _nVar,_precision, _invAroe,y);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<< endl<<"ProblemFluid::InvMatriceRoe : Valeurs propres :" << nbVp_dist<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<valeurs_propres_dist[ct].real()<< ", " <<valeurs_propres_dist[ct].imag() <<")  ";
		cout<< endl;
		cout<<"ProblemFluid::InvMatriceRoe : Valeurs à atteindre pour l'interpolation"<<endl;
		for(int ct =0; ct<nbVp_dist; ct++)
			cout<< "("<<y[ct].real()<< ", " <<y[ct].imag() <<")  ";
		cout<< endl;
		displayMatrix(_invAroe,_nVar,"Inverse matrice de Roe");
	}
}

void ProblemFluid::terminate(){

	delete[] _AroePlus;
	delete[] _Diffusion;
	delete[] _GravityImplicitationMatrix;
	delete[] _AroeMinus;
	delete[] _Aroe;
	delete[] _absAroe;
	delete[] _signAroe;
	delete[] _invAroe;
	delete[] _AroeImplicit;
	delete[] _AroeMinusImplicit;
	delete[] _AroePlusImplicit;
	delete[] _absAroeImplicit;
	delete[] _phi;
	delete[] _Jcb;
	delete[] _JcbDiff;
	delete[] _a;
	delete[] _primToConsJacoMat;

	delete[] _l;
	delete[] _r;
	delete[] _Uroe;
	delete[] _Udiff;
	delete[] _temp;
	delete[] _externalStates;
	delete[] _idm;
	delete[] _idn;
	delete[] _vec_normal;
	delete[] _Ui;
	delete[] _Uj;
	delete[] _Vi;
	delete[] _Vj;
	if(_nonLinearFormulation==VFRoe){
		delete[] _Uij;
		delete[] _Vij;
	}
	delete[] _Si;
	delete[] _Sj;
	delete[] _pressureLossVector;
	delete[] _porosityGradientSourceVector;
	if(_isScaling)
	{
		delete[] _blockDiag;
		delete[] _invBlockDiag;

		VecDestroy(&_vecScaling);
		VecDestroy(&_invVecScaling);
		VecDestroy(&_bScaling);
	}

	VecDestroy(&_conservativeVars);
	VecDestroy(&_newtonVariation);
	VecDestroy(&_b);
	VecDestroy(&_primitiveVars);
	VecDestroy(&_Uext);
	VecDestroy(&_Vext);
	VecDestroy(&_Uextdiff);

	// 	PCDestroy(_pc);
	KSPDestroy(&_ksp);
	for(int i=0;i<_nbPhases;i++)
		delete _fluides[i];
}
