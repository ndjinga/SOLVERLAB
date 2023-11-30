/*
 * WaveStaggered.cxx
 */

#include "WaveStaggered.hxx"
#include "StiffenedGas.hxx"

using namespace std;

void
computeVelocityMCells(const Field& velocity,
                      Field& velocityMCells)
{
    Mesh myMesh=velocity.getMesh();
    int nbCells=myMesh.getNumberOfCells();

    for(int i=0;i<nbCells;i++)
    {
        std::vector< int > facesId=myMesh.getCell(i).getFacesId();
        velocityMCells(i)=(velocity(facesId[0])+velocity(facesId[1]))/2.;
    }
}

WaveStaggered::WaveStaggered(phaseType fluid, pressureEstimate pEstimate, int dim, bool useDellacherieEOS){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1;
	_dragCoeffs=vector<double>(1,0);
	_fluides.resize(1);
	_useDellacherieEOS=useDellacherieEOS;

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		if(fluid==Gas){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(1.4,743,300,2.22e5);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, e=2.22e5, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_fluides[0] = new StiffenedGas(996,1e5,300,1.12e5,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C, e=1.12e5, c_v=4130
		}
	}
	else{//EOS at 155 bars and 618K 
		if(fluid==Gas){
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_fluides[0] = new StiffenedGas(102,1.55e7,618,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
		}
		else{//To do : change to normal regime: 155 bars and 573K
			cout<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 618 K (345°C)"<<endl;
			if(_useDellacherieEOS)
				_fluides[0]= new StiffenedGasDellacherie(2.35,1e9,-1.167e6,1816); //stiffened gas law for water from S. Dellacherie
			else
				_fluides[0]= new StiffenedGas(594.,1.55e7,618.,1.6e6, 621.,3100.);  //stiffened gas law for water at pressure 155 bar, and temperature 345°C
		}
	}
}
void WaveStaggered::initialize(){
	cout<<"\n Initialising the Wave System model\n"<<endl;
	*_runLogFile<<"\n Initialising the Wave Sytem model\n"<<endl;

	_globalNbUnknowns = (_nVar-1)*_Nmailles + _Nfaces;//Staggered discretisation : velocity is on faces

	ProblemFluid::initialize();
}

double WaveStaggered::computeTimeStep(bool & stop){//dt is not known and will not contribute to the Newton scheme

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "ProblemFluid::computeTimeStep : Début calcul matrice implicite et second membre"<<endl;
		cout << endl;
	}
	if(_restartWithNewTimeScheme)//This is a change of time scheme during a simulation
	{
		if(_timeScheme == Implicit)
			MatCreateSeqBAIJ(PETSC_COMM_SELF, _nVar, _nVar*_Nmailles, _nVar*_Nmailles, (1+_neibMaxNbCells), PETSC_NULL, &_A);			
		else
			MatDestroy(&_A);
		_restartWithNewTimeScheme=false;
	}
	if(_timeScheme == Implicit)
		MatZeroEntries(_A);

	VecAssemblyBegin(_b);
	VecZeroEntries(_b);

	std::vector< int > idCells(2);
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
					//calcul et insertion de A^-*Jcb
					Polynoms::matrixProduct(_AroeMinusImplicit, _nVar, _nVar, _Jcb, _nVar, _nVar, _a);
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
					Polynoms::matrixProduct(_Diffusion, _nVar, _nVar, _JcbDiff, _nVar, _nVar, _a);
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

			addConvectionToSecondMember(idCells[0],idCells[1], false);
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
				//cout<<"idm= "<<idm<<"idn= "<<idn<<"nbvoismax= "<<_neibMaxNbCells<<endl;
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
				cout<<"lattice mesh junction at face "<<j<<" nbvoismax= "<<_neibMaxNbCells<<endl;
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
	VecAssemblyEnd(_b);

	if(_timeScheme == Implicit){
		for(int imaille = 0; imaille<_Nmailles; imaille++)
			MatSetValuesBlocked(_A, size, &imaille, size, &imaille, _GravityImplicitationMatrix, ADD_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
			displayMatrix(_GravityImplicitationMatrix,_nVar,"Gravity matrix:");

		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
		if(_verbose && _nbTimeStep%_freqSave ==0){
			cout << "ProblemFluid::computeTimeStep : Fin calcul matrice implicite et second membre"<<endl;
			cout << "ProblemFluid::computeTimeStep : Matrice implicite :"<<endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << "ProblemFluid::computeTimeStep : Second membre :"<<endl;
			VecView(_b,  PETSC_VIEWER_STDOUT_WORLD);
			cout << endl;
		}
	}

	stop=false;

	/*
	if(_nbTimeStep+1<_cfl)
		return (_nbTimeStep+1)*_minl/_maxvp;
	else
	 */
	if(_maxvp>0)
		return _cfl*_minl/_maxvp;
	else//case of incompressible fluid at rest. Use a velocity of 1
		return _cfl*_minl;
}
bool WaveStaggered::iterateTimeStep(bool &converged)
{
	if(_timeScheme == Explicit || !_usePrimitiveVarsInNewton)
		ProblemFluid::iterateTimeStep(converged);
	else
	{
		bool stop=false;

		if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
			_maxvp=0;
			computeTimeStep(stop);
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
			if(_PetscIts>=_maxPetscIts)//solving the linear system failed
			{
				cout<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
				*_runLogFile<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
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
						VecGetValues(_primitiveVars, 1, &I, &x);
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

		double relaxation=1;//Vk+1=Vk+relaxation*deltaV

		VecAXPY(_primitiveVars,  relaxation, _newtonVariation);

		//mise a jour du champ primitif
		updateConservatives();

		if(_nbPhases==2 && fabs(_err_press_max) > _precision)//la pression n'a pu être calculée en diphasique à partir des variables conservatives
		{
			cout<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
			*_runLogFile<<"Warning consToPrim: nbiter max atteint, erreur relative pression= "<<_err_press_max<<" precision= " <<_precision<<endl;
			converged=false;
			return false;
		}
		if(_system)
		{
			cout<<"Vecteur Vkp1-Vk "<<endl;
			VecView(_newtonVariation,  PETSC_VIEWER_STDOUT_SELF);
			cout << "Nouvel etat courant Vk de l'iteration Newton: " << endl;
			VecView(_primitiveVars,  PETSC_VIEWER_STDOUT_SELF);
		}

		if(_nbPhases==2 && (_nbTimeStep-1)%_freqSave ==0){
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
}
void WaveStaggered::computeNewtonVariation()
{
	if(!_usePrimitiveVarsInNewton)
		ProblemFluid::computeNewtonVariation();
	else
	{
		if(_verbose)
		{
			cout<<"Vecteur courant Vk "<<endl;
			VecView(_primitiveVars,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Matrice du système linéaire avant contribution delta t" << endl;
			MatView(_A,PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
			cout << "Second membre du système linéaire avant contribution delta t" << endl;
			VecView(_b, PETSC_VIEWER_STDOUT_SELF);
			cout << endl;
		}
		if(_timeScheme == Explicit)
		{
			VecCopy(_b,_newtonVariation);
			VecScale(_newtonVariation, _dt);
			if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
			{
				cout<<"Vecteur _newtonVariation =_b*dt"<<endl;
				VecView(_newtonVariation,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
		}
		else
		{
			MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);

			VecAXPY(_b, 1/_dt, _old);
			VecAXPY(_b, -1/_dt, _conservativeVars);

			for(int imaille = 0; imaille<_Nmailles; imaille++){
				_idm[0] = _nVar*imaille;
				for(int k=1; k<_nVar; k++)
					_idm[k] = _idm[k-1] + 1;
				VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
				primToConsJacobianMatrix(_Vi);
				for(int k=0; k<_nVar*_nVar; k++)
					_primToConsJacoMat[k]*=1/_dt;
				MatSetValuesBlocked(_A, 1, &imaille, 1, &imaille, _primToConsJacoMat, ADD_VALUES);
			}
			MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

#if PETSC_VERSION_GREATER_3_5
			KSPSetOperators(_ksp, _A, _A);
#else
			KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

			if(_verbose)
			{
				cout << "Matrice du système linéaire" << endl;
				MatView(_A,PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
				cout << "Second membre du système linéaire" << endl;
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
			if(_system)
			{
				cout << "solution du systeme lineaire local:" << endl;
				VecView(_newtonVariation, PETSC_VIEWER_STDOUT_SELF);
				cout << endl;
			}
		}
	}
}
void WaveStaggered::convectionState( const long &i, const long &j, const bool &IsBord){
}

void WaveStaggered::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	int k;
	double v2=0, q_n=0;//q_n=quantité de mouvement normale à la face frontière;
	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+1)]*normale[k];

	double porosityj=_porosityField(j);

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< " face unit normal vector "<<endl;
		for(k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall){
		//Pour la convection, inversion du sens de la vitesse normale
		for(k=0; k<_Ndim; k++)
			_externalStates[(k+1)]-= 2*q_n*normale[k];

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		//Pour la diffusion, paroi à vitesse et temperature imposees
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double pression=_externalStates[0];
		double T=_limitField[nameOfGroup].T;
		double rho=_fluides[0]->getDensity(pression,T);

		_externalStates[0]=porosityj*rho;
		_externalStates[1]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
		v2 +=_limitField[nameOfGroup].v_x[0]*_limitField[nameOfGroup].v_x[0];
		if(_Ndim>1)
		{
			v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
			if(_Ndim==3)
			{
				_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
				v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
			}
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){

		if(q_n<=0){
			VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_fluides[0]->getDensity(pression,T);

			_externalStates[0]=porosityj*rho;
			_externalStates[1]=_externalStates[0]*(_limitField[nameOfGroup].v_x[0]);
			v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
			if(_Ndim>1)
			{
				v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
				_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
				if(_Ndim==3)
				{
					_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
					v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		}
		else if((_nbTimeStep-1)%_freqSave ==0)
			cout<< "Warning : fluid going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		//Building the external state
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		if(q_n<=0){
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress,_limitField[nameOfGroup].T);
		}
		else{
			if((_nbTimeStep-1)%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature"<<endl;
			_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		}

		for(k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);


		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		if(q_n<=0 &&  (_nbTimeStep-1)%_freqSave ==0)
			cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature"<<endl;

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
		{
			cout<<"Cond lim outlet densite= "<<_externalStates[0]<<" gravite= "<<_GravityField3d[0]<<" Cj.x()= "<<Cj.x()<<endl;
			cout<<"Cond lim outlet pression ref= "<<_limitField[nameOfGroup].p<<" pression hydro= "<<hydroPress<<" total= "<<_limitField[nameOfGroup].p+hydroPress<<endl;
		}
		//Building the external state
		_idm[0] = _nVar*j;// Kieu
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		_externalStates[0]=porosityj*_fluides[0]->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		for(k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}
}

void WaveStaggered::convectionMatrices()
{
	//entree: URoe = rho, u, H
	//sortie: matrices Roe+  et Roe-

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
		cout<<"WaveStaggered::convectionMatrices()"<<endl;

	double u_n=0, u_2=0;//vitesse normale et carré du module

	for(int i=0;i<_Ndim;i++)
	{
		u_2 += _Uroe[1+i]*_Uroe[1+i];
		u_n += _Uroe[1+i]*_vec_normal[i];
	}

	vector<complex<double> > vp_dist(3,0);

	if(_spaceScheme==staggered && _nonLinearFormulation==VFFC)//special case
	{
		staggeredVFFCMatricesConservativeVariables(u_n);//Computation of classical upwinding matrices
		if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)//For use in implicit matrix
			staggeredVFFCMatricesPrimitiveVariables(u_n);
	}
	else
	{
		Vector vitesse(_Ndim);
		for(int idim=0;idim<_Ndim;idim++)
			vitesse[idim]=_Uroe[1+idim];

		double  c, H, K, k;
		/***********Calcul des valeurs propres ********/
		H = _Uroe[_nVar-1];
		c = _fluides[0]->vitesseSonEnthalpie(H-u_2/2);//vitesse du son a l'interface
		k = _fluides[0]->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
		K = u_2*k/2; //g-1/2 *|u|²

		vp_dist[0]=u_n-c;vp_dist[1]=u_n;vp_dist[2]=u_n+c;

		_maxvploc=fabs(u_n)+c;
		if(_maxvploc>_maxvp)
			_maxvp=_maxvploc;

		if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
			cout<<"WaveStaggered::convectionMatrices Eigenvalues "<<u_n-c<<" , "<<u_n<<" , "<<u_n+c<<endl;

		RoeMatrixConservativeVariables( u_n, H,vitesse,k,K);

		/******** Construction des matrices de decentrement ********/
		if( _spaceScheme ==centered){
			if(_entropicCorrection)
			{
				*_runLogFile<<"WaveStaggered::convectionMatrices: entropy scheme not available for centered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("WaveStaggered::convectionMatrices: entropy scheme not available for centered scheme");
			}

			for(int i=0; i<_nVar*_nVar;i++)
				_absAroe[i] = 0;
		}
		else if(_spaceScheme == upwind || _spaceScheme ==pressureCorrection || _spaceScheme ==lowMach){
			if(_entropicCorrection)
				entropicShift(_vec_normal);
			else
				_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

			vector< complex< double > > y (3,0);
			for( int i=0 ; i<3 ; i++)
				y[i] = Polynoms::abs_generalise(vp_dist[i])+_entropicShift[i];
			Polynoms::abs_par_interp_directe(3,vp_dist, _Aroe, _nVar,_precision, _absAroe,y);

			if( _spaceScheme ==pressureCorrection)
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(1+i)*_nVar+1+j]-=(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			else if( _spaceScheme ==lowMach){
				double M=sqrt(u_2)/c;
				for( int i=0 ; i<_Ndim ; i++)
					for( int j=0 ; j<_Ndim ; j++)
						_absAroe[(1+i)*_nVar+1+j]-=(1-M)*(vp_dist[2].real()-vp_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
			}
		}
		else if( _spaceScheme ==staggered ){
			if(_entropicCorrection)//To do: study entropic correction for staggered
			{
				*_runLogFile<<"WaveStaggered::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("WaveStaggered::convectionMatrices: entropy scheme not available for staggered scheme");
			}

			staggeredRoeUpwindingMatrixConservativeVariables( u_n, H, vitesse, k, K);
		}
		else
		{
			*_runLogFile<<"WaveStaggered::convectionMatrices: scheme not treated"<<endl;
			_runLogFile->close();
			throw CdmathException("WaveStaggered::convectionMatrices: scheme not treated");
		}

		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinus[i] = (_Aroe[i]-_absAroe[i])/2;
			_AroePlus[i]  = (_Aroe[i]+_absAroe[i])/2;
		}
		if(_timeScheme==Implicit)
		{
			if(_usePrimitiveVarsInNewton)//Implicitation using primitive variables
			{
				_Vij[0]=_fluides[0]->getPressureFromEnthalpy(_Uroe[_nVar-1]-u_2/2, _Uroe[0]);//pressure
				_Vij[_nVar-1]=_fluides[0]->getTemperatureFromPressure( _Vij[0], _Uroe[0]);//Temperature
				for(int idim=0;idim<_Ndim; idim++)
					_Vij[1+idim]=_Uroe[1+idim];
				primToConsJacobianMatrix(_Vij);
				Polynoms::matrixProduct(_AroeMinus, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroeMinusImplicit);
				Polynoms::matrixProduct(_AroePlus,  _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _AroePlusImplicit);
			}
			else
				for(int i=0; i<_nVar*_nVar;i++)
				{
					_AroeMinusImplicit[i] = _AroeMinus[i];
					_AroePlusImplicit[i]  = _AroePlus[i];
				}
		}
		if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
		{
			displayMatrix(_Aroe, _nVar,"Matrice de Roe");
			displayMatrix(_absAroe, _nVar,"Valeur absolue matrice de Roe");
			displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus");
			displayMatrix(_AroePlus, _nVar,"Matrice _AroePlus");
		}
	}

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0 && _timeScheme==Implicit)
	{
		displayMatrix(_AroeMinusImplicit, _nVar,"Matrice _AroeMinusImplicit");
		displayMatrix(_AroePlusImplicit,  _nVar,"Matrice _AroePlusImplicit");
	}

	/*********Calcul de la matrice signe pour VFFC, VFRoe et décentrement des termes source*****/
	if(_entropicCorrection)
	{
		InvMatriceRoe( vp_dist);
		Polynoms::matrixProduct(_absAroe, _nVar, _nVar, _invAroe, _nVar, _nVar, _signAroe);
	}
	else if (_spaceScheme==upwind || (_spaceScheme ==pressureCorrection ) || (_spaceScheme ==lowMach ))//upwind sans entropic
		SigneMatriceRoe( vp_dist);
	else if (_spaceScheme==centered)//centre  sans entropic
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
	else if( _spaceScheme ==staggered )//à tester
	{
		double signu;
		if(u_n>0)
			signu=1;
		else if (u_n<0)
			signu=-1;
		else
			signu=0;
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
		_signAroe[0] = signu;
		for(int i=1; i<_nVar-1;i++)
			_signAroe[i*_nVar+i] = -signu;
		_signAroe[_nVar*(_nVar-1)+_nVar-1] = signu;
	}
	else
	{
		*_runLogFile<<"WaveStaggered::convectionMatrices: well balanced option not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("WaveStaggered::convectionMatrices: well balanced option not treated");
	}
}
void WaveStaggered::computeScaling(double maxvp)
{
	_blockDiag[0]=1;
	_invBlockDiag[0]=1;
	for(int q=1; q<_nVar-1; q++)
	{
		_blockDiag[q]=1./maxvp;//
		_invBlockDiag[q]= maxvp;//1.;//
	}
	_blockDiag[_nVar - 1]=(_fluides[0]->constante("gamma")-1)/(maxvp*maxvp);//1
	_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
}


void WaveStaggered::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{}

void WaveStaggered::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{}

void WaveStaggered::porosityGradientSourceVector()
{}

void WaveStaggered::jacobian(const int &j, string nameOfGroup,double * normale)
{}


Vector WaveStaggered::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
	{
		cout<<"WaveStaggered::convectionFlux start"<<endl;
		cout<<"Ucons"<<endl;
		cout<<U<<endl;
		cout<<"Vprim"<<endl;
		cout<<V<<endl;
	}

	double phirho=U(0);//phi rho
	Vector phiq(_Ndim);//phi rho u
	for(int i=0;i<_Ndim;i++)
		phiq(i)=U(1+i);

	double pression=V(0);
	Vector vitesse(_Ndim);
	for(int i=0;i<_Ndim;i++)
		vitesse(i)=V(1+i);
	double Temperature= V(1+_Ndim);

	double vitessen=vitesse*normale;
	double rho=phirho/porosity;
	double e_int=_fluides[0]->getInternalEnergy(Temperature,rho);

	Vector F(_nVar);
	F(0)=phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F(1+i)=phirho*vitessen*vitesse(i)+pression*normale(i)*porosity;
	F(1+_Ndim)=phirho*(e_int+0.5*vitesse*vitesse+pression/rho)*vitessen;

	if(_verbose && (_nbTimeStep-1)%_freqSave ==0)
	{
		cout<<"WaveStaggered::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void WaveStaggered::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
{}

void WaveStaggered::getDensityDerivatives( double pressure, double temperature, double v2)
{}
void WaveStaggered::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/WaveStaggeredPrim_");///Results
	string cons(_path+"/WaveStaggeredCons_");
	prim+=_fileName;
	cons+=_fileName;

	PetscInt Ii;
	for (long i = 0; i < _Nmailles; i++){
		// j = 0 : pressure; j = _nVar - 1: temperature; j = 1,..,_nVar-2: velocity
		for (int j = 0; j < _nVar; j++){
			Ii = i*_nVar +j;
			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
		}
	}
	if(_saveConservativeField){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : density; j = _nVar - 1 : energy j = 1,..,_nVar-2: momentum
			for (int j = 0; j < _nVar; j++){
				Ii = i*_nVar +j;
				VecGetValues(_conservativeVars,1,&Ii,&_UU(i,j));
			}
		}
		_UU.setTime(_time,_nbTimeStep);
	}
	_VV.setTime(_time,_nbTimeStep);

	// create mesh and component info
	if (_nbTimeStep ==0){
		string prim_suppress ="rm -rf "+prim+"_*";
		string cons_suppress ="rm -rf "+cons+"_*";

		system(prim_suppress.c_str());//Nettoyage des précédents calculs identiques
		system(cons_suppress.c_str());//Nettoyage des précédents calculs identiques

		if(_saveConservativeField){
			_UU.setInfoOnComponent(0,"Density_(kg/m^3)");
			_UU.setInfoOnComponent(1,"Momentum_x");// (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(2,"Momentum_y");// (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(3,"Momentum_z");// (kg/m^2/s)

			_UU.setInfoOnComponent(_nVar-1,"Energy_(J/m^3)");

			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons);
				break;
			case MED :
				_UU.writeMED(cons);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}
		_VV.setInfoOnComponent(0,"Pressure_(Pa)");
		_VV.setInfoOnComponent(1,"Velocity_x_(m/s)");
		if (_Ndim>1)
			_VV.setInfoOnComponent(2,"Velocity_y_(m/s)");
		if (_Ndim>2)
			_VV.setInfoOnComponent(3,"Velocity_z_(m/s)");
		_VV.setInfoOnComponent(_nVar-1,"Temperature_(K)");

		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim);
			break;
		case MED :
			_VV.writeMED(prim);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}
	}
	// do not create mesh
	else{
		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim,false);
			break;
		case MED :
			_VV.writeMED(prim,false);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}
		if(_saveConservativeField){
			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons,false);
				break;
			case MED :
				_UU.writeMED(cons,false);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}
	}
	if(_saveVelocity){
		for (long i = 0; i < _Nmailles; i++){
			// j = 0 : pressure; j = _nVar - 1: temperature; j = 1,..,_nVar-2: velocity
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +1+j;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse(i,j));
			}
			for (int j = _Ndim; j < 3; j++)//On met à zero les composantes de vitesse si la dimension est <3
				_Vitesse(i,j)=0;
		}
		_Vitesse.setTime(_time,_nbTimeStep);
		if (_nbTimeStep ==0){
			_Vitesse.setInfoOnComponent(0,"Velocity_x_(m/s)");
			_Vitesse.setInfoOnComponent(1,"Velocity_y_(m/s)");
			_Vitesse.setInfoOnComponent(2,"Velocity_z_(m/s)");

			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity",false);
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity",false);
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
	}
	if(_isStationary)
	{
		prim+="_Stat";
		cons+="_Stat";

		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(prim);
			break;
		case MED :
			_VV.writeMED(prim);
			break;
		case CSV :
			_VV.writeCSV(prim);
			break;
		}

		if(_saveConservativeField){
			switch(_saveFormat)
			{
			case VTK :
				_UU.writeVTK(cons);
				break;
			case MED :
				_UU.writeMED(cons);
				break;
			case CSV :
				_UU.writeCSV(cons);
				break;
			}
		}

		if(_saveVelocity){
			switch(_saveFormat)
			{
			case VTK :
				_Vitesse.writeVTK(prim+"_Velocity");
				break;
			case MED :
				_Vitesse.writeMED(prim+"_Velocity");
				break;
			case CSV :
				_Vitesse.writeCSV(prim+"_Velocity");
				break;
			}
		}
	}
}
