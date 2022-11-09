/*
 * IsothermalSinglePhase.cxx
 *
 *  Created on: July 7 2022
 */

#include<bits/stdc++.h>
#include "IsothermalSinglePhase.hxx"

using namespace std;

IsothermalSinglePhase::IsothermalSinglePhase(phaseType fluid, pressureEstimate pEstimate, int dim, bool isCompressibleFluid){
	_Ndim=dim;
	_nVar=_Ndim+1;
	_nbPhases = 1;
	_dragCoeffs=vector<double>(1,0);
	_isCompressibleFluid = isCompressibleFluid;
	_fluides.resize(1);

	if (pEstimate==around1bar300K)//EOS at 1 bar and 300K
	{
		_Temperature=300;//Constant temperature of the model
		if(fluid==Gas)
		{
			if(isCompressibleFluid)
			{
				cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
				*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
				_internalEnergy=2.22e5;//nitrogen internal energy at 1bar, 300K
				_fluides[0] = new StiffenedGas(1.4,743,_Temperature,_internalEnergy);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, c_v=743
			}
			else
			{
				cout<<"Fluid is incompressible air around 1 bar and 300 K (27°C)"<<endl;
				*_runLogFile<<"Fluid is incompressible air around 1 bar and 300 K (27°C)"<<endl;
				_fluides[0] = new IncompressibleFluid(1.12);  
			}
		}
		else
		{
			if(isCompressibleFluid)
			{
				cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
				*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
				_internalEnergy=1.12e5;//water internal energy at 1 bar, 300K
				_fluides[0] = new StiffenedGas(996.56,1e5,_Temperature,_internalEnergy,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C
			}
			else
			{
				cout<<"Fluid is incompressible water around 1 bar and 300 K (27°C)"<<endl;
				*_runLogFile<<"Fluid is incompressible water around 1 bar and 300 K (27°C)"<<endl;
				_fluides[0] = new IncompressibleFluid(996.56);  
			}
		}
	}
	else
	{
		if(fluid==Gas)//EOS at 155 bars and 618K
		{
			if(isCompressibleFluid)
			{
				cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
				*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
				_Temperature=618;//Constant temperature of the model
				_internalEnergy=2.44e6;//Gas internal energy at saturation at 155 bar
				_fluides[0] = new StiffenedGas(102,1.55e7,_Temperature,_internalEnergy, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C:
			}
			else
			{
				cout<<"Fluid is incompressible air around 155 bars and 618 K (345°C)"<<endl;
				*_runLogFile<<"Fluid is incompressible air around 155 bars and 618 K (345°C)"<<endl;
				_fluides[0] = new IncompressibleFluid(102);  
			}
		}
		else//EOS at 155 bars and 573K
		{
			if(isCompressibleFluid)
			{
				cout<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
				*_runLogFile<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
				_Temperature=573;//Constant temperature of the model
				_internalEnergy=1.3e6;//water internal energy at saturation at 155 bar
				_fluides[0] = new StiffenedGas(726.82,1.55e7,_Temperature,_internalEnergy, 971.,5454.);  //stiffened gas law for water at pressure 155 bar and temperature 300°C:
			}
			else
			{
				cout<<"Fluid is incompressible water around 155 bars and 573 K (300°C)"<<endl;
				*_runLogFile<<"Fluid is incompressible water around 155 bars and 573 K (300°C)"<<endl;
				_fluides[0] = new IncompressibleFluid(594);  
			}
		}
	}

	_compressibleFluid=dynamic_cast< CompressibleFluid * >(_fluides[0] );//To be used when the numerical method requires a compressible fluid e.g. upwind
    _usePrimitiveVarsInNewton=true;//This class is designed only to solve linear system in primitive variables
    _Vdiff=NULL;
    _saveAllFields = false;
	_nonLinearFormulation=reducedRoe;//Only case implemented is reduced roe

	_fileName = "SolverlabIsothermalSinglePhase";
    PetscPrintf(PETSC_COMM_WORLD,"\n Isothermal single phase problem \n");
    }

void IsothermalSinglePhase::initialize(){
	cout<<"\n Initialising the isothermal single phase model (memory allocations for matrices, vectors and fields) \n"<<endl;
	*_runLogFile<<"\n Initialising the isothermal single phase model (memory allocations for matrices, vectors and fields) \n"<<endl;

	_Uroe = new double[_nVar+1];//Deleted in ProblemFluid::terminate()
	_Vextdiff= new double[_nVar];
	_Vext= new double[_nVar];

	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+1]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];//Deleted in ProblemFluid::terminate()

	_Vdiff = new double[_nVar];
	
	if(_saveVelocity || _saveAllFields)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant

	if(_saveAllFields)
	{
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		_Pressure.setInfoOnComponent(0,"Pressure_(Pa)");
		_Density=Field("Density",CELLS,_mesh,1);
		_Density.setInfoOnComponent(0,"Density_(kg/m^3)");
		_MachNumber=Field("MachNumber",CELLS,_mesh,1);
		_MachNumber.setInfoOnComponent(0,"Mach number");
		_VitesseX=Field("Velocity x",CELLS,_mesh,1);
		_VitesseX.setInfoOnComponent(0,"Velocity_x_(m/s)");
		if(_Ndim>1)
		{
			_VitesseY=Field("Velocity y",CELLS,_mesh,1);
			_VitesseY.setInfoOnComponent(0,"Velocity_y_(m/s)");
			if(_Ndim>2)
				_VitesseZ=Field("Velocity z",CELLS,_mesh,1);
				_VitesseZ.setInfoOnComponent(0,"Velocity_z_(m/s)");
		}
	}

	/* In case of an incompressible fluid, pressure may be defined up to a constant hence a singular inear system (no uniqueness unless outlet boundary condition) */
	if(_isCompressibleFluid)
		_isSingularSystem=false;
	else
	{
		_isSingularSystem=true;//unless outlet boundary condition
		/* Check the presence of an unless outlet boundary condition */
		for(  map<string, LimitField>::iterator it = _limitField.begin() ; it!=_limitField.end() ; it++ )
		{
			if( (it->second).bcType == Outlet or (it->second).bcType == InletPressure )
			{
				_isSingularSystem=false;
				break;
			}
		}
	if(_verbose && _nbTimeStep%_freqSave ==0)	
		if(_isSingularSystem)
			cout<<"!!!!!!!########## singular system ###########!!!!!!!!!!"<<endl;
		else
			cout<<"!!!!!!!########## non singular system ###########!!!!!!!!!!"<<endl;
	}
	
	ProblemFluid::initialize();

	/* Deal with the particular case of singular systems */
	if( _isSingularSystem )
	{
		cout<<"No pressure imposed at the boundary : the PDEs solution pressure is defined up to a constant."<<endl;
		cout<<"Singular linear system : the discrete pressure will be computed to have a zero mean."<<endl;
		/* Build vector in the kernel of the system matrix */
		VecDuplicate(_conservativeVars, &_constantPressureVector);//Vector _constantPressureVector has same parallel structure as _conservativeVars
		VecZeroEntries(_constantPressureVector);
		IS pressureCoeffsIS;
		ISCreateStride(PETSC_COMM_SELF, _Nmailles,0,_nVar, &pressureCoeffsIS);
		VecISSet(_constantPressureVector, pressureCoeffsIS,1./sqrt(_Nmailles));
		VecAssemblyBegin(_constantPressureVector);
		VecAssemblyEnd(_constantPressureVector);
		ISDestroy(&pressureCoeffsIS);
		
		/* Give kernel vector to the system matrix */
		MatNullSpace nullsp;
		MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_FALSE, 1, &_constantPressureVector, &nullsp);//Declaration of a kernel containing a vector of constant pressure but NOT containing constant vectors
		//MatSetTransposeNullSpace(_A, nullsp);//To be used if the linear solver needs the kernel of the transpose matrix
		MatSetNullSpace(_A, nullsp);
		MatNullSpaceDestroy(&nullsp);
	}
	//PCFactorSetShiftType(_pc,MAT_SHIFT_INBLOCKS);//To be used in case of a zero pivot
	//PCFactorSetShiftAmount(_pc,1e-10);//To be used in case of a zero pivot
	//PCFactorSetColumnPivot(_pc,0.80);
	//PCFactorSetZeroPivot(_pc,1e-5);
	//PCFactorReorderForNonzeroDiagonal(_pc,1e-5);
}

void IsothermalSinglePhase::terminate(){
	delete[] _Vdiff,_Vextdiff,_Vext;
	
	if( _isSingularSystem )
		VecDestroy(&_constantPressureVector);

	ProblemFluid::terminate();
}

void IsothermalSinglePhase::convectionState( const long &i, const long &j, const bool &IsBord){
	//entree: indices des cellules Vi et Vj.
	//Attention : Vj peut être un état fantôme
	//sortie: URoe = rho, u, 1/c^2

	//Extraction of primitive states
	_idm[0] = _nVar*i; 
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Convection state: variables primitives maille " << i<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q] << endl;
		cout << endl;
	}

	if(!IsBord ){
		for(int k=0; k<_nVar; k++)
			_idn[k] = _nVar*j + k;
		VecGetValues(_primitiveVars, _nVar, _idn, _Vj);
	}
	else
		for(int k=0; k<_nVar; k++)
			_Vj[k]=_Vext[k];

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Convection state: variables primitives maille " <<j <<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q] << endl;
		cout << endl;
	}
	assert(_Vi[0]>0);
	assert(_Vj[0]>0);

	//Computation of conservative states Ui and Uj
	primToCons(_Vi,0,_Ui,0);
	primToCons(_Vj,0,_Uj,0);

	//Computation of Roe density and velocity
	PetscScalar ri, rj, xi, xj, pi, pj;
	PetscInt Ii;
	ri = sqrt(_Ui[0]);//racine carre de phi_i rho_i
	rj = sqrt(_Uj[0]);
	_Uroe[0] = ri*rj;	//moyenne geometrique des densites
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Densité moyenne Roe  gauche " << i << ": " << _Ui[0] << ", droite " << j << ": " << _Uj[0] << "->" << _Uroe[0] << endl;
	for(int k=0;k<_Ndim;k++){
		xi = _Ui[k+1];
		xj = _Uj[k+1];
		_Uroe[1+k] = (xi/ri + xj/rj)/(ri + rj);
		//"moyenne" des vitesses
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout << "Vitesse de Roe composante "<< k<<"  gauche " << i << ": " << xi/_Ui[0] << ", droite " << j << ": " << xj/_Uj[0] << "->" << _Uroe[k+1] << endl;
	}

	//Computation of 1/c²// Todo :  add porosity in the sound speed formula
	if(abs(_Vi[0]-_Vj[0])<_precision*max(abs(_Vi[0]),abs(_Vj[0])))//need to use EOS
	{
		if(not _isCompressibleFluid)//case of an incompressible fluid
		    _Uroe[_nVar] = 0;
		else
		{
		    _Uroe[_nVar]  = _compressibleFluid->getInverseSoundSpeed(max(_Vi[0],_Vj[0]), _Temperature);//store 1/c
		    _Uroe[_nVar] *= _Uroe[_nVar];//store 1/c^2
		}
    }    
	else
	    _Uroe[_nVar] = (_Ui[0]-_Uj[0])/(_Vi[0]-_Vj[0]);//_Uroe has size _nVar+1 !!!
	
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar+1;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
}

void IsothermalSinglePhase::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//entree: indices des cellules Vi et Vj.
	//Attention : Vj peut être un état fantôme
	//sortie: matrices et etat de diffusion (rho, q) ou (p, v) ?
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

	if(IsBord)
		for(int k=0; k<_nVar; k++)
			_Vj[k]=_Vextdiff[k];
	else
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices cellule gauche = " << i << endl;
		cout << "Vi = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q]  << "\t";
		cout << endl;
		cout << "IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices cellule droite =" << j << endl;
		cout << "Vj = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q]  << "\t";
		cout << endl;
	}

	for(int k=0; k<_nVar; k++)
		_Vdiff[k] = (_Vi[k]/_porosityi+_Vj[k]/_porosityj)/2;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices primitive diffusion state" << endl;
		cout << "_Vdiff = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vdiff[q]  << "\t";
		cout << endl;
		cout << "porosite gauche= "<<_porosityi<< ", porosite droite= "<<_porosityj<<endl;
	}

	if(_timeScheme==Implicit)
	{
		double mu = _fluides[0]->getViscosity(_Temperature);
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		for(int i=1;i<_nVar;i++)
			_Diffusion[i*_nVar+i] = mu;
	}
}

void IsothermalSinglePhase::convectionMatrices()
{
	//entree: URoe = rho, u, 1/c^2
	//sortie: matrices Roe+  et Roe-

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"IsothermalSinglePhase::convectionMatrices()"<<endl;

	double u_n=0, u_2=0;//vitesse normale et carré du module
	Vector vitesse(_Ndim);

	for(int idim=0; idim <_Ndim; idim++)
	{
		u_2 += _Uroe[1+idim]*_Uroe[1+idim];
		u_n += _Uroe[1+idim]*_vec_normal[idim];
		vitesse[idim]=_Uroe[1+idim];
	}

	vector<std::complex<double>> vp_dist;
	int nb_vp_dist;
	double c;//store the sound speed
	if(_Uroe[_nVar]==0.)//infinite sound speed
		if(_timeScheme==Explicit)
			throw CdmathException("Explicit scheme cannot be used for incompressible fluids since dt=0");
		else
		{
			nb_vp_dist=1;
			vp_dist.resize(nb_vp_dist);
			c=0.;//The velocity will be used to determine the time step
			vp_dist[0]=u_n;//store the fluid velocity for calculation of the time step
		}
	else
	{
		c=1./sqrt(_Uroe[_nVar]);
		if(_Ndim==1)//In 1D two acoustic eigenvalues
		{
			nb_vp_dist=2;
			vp_dist.resize(nb_vp_dist);
			vp_dist[0]=u_n-c;
			vp_dist[1]=u_n+c;
		}
		else//In 2D and 3D extra shear eigenvalue u_n
		{
			nb_vp_dist=3;
			vp_dist.resize(nb_vp_dist);
			vp_dist[0]=u_n-c;
			vp_dist[1]=u_n;
			vp_dist[2]=u_n+c;
		}
	}
	
	_maxvploc=fabs(u_n)+c;
	if(_maxvploc>_maxvp)
		_maxvp=_maxvploc;

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"IsothermalSinglePhase::convectionMatrices, "<< nb_vp_dist<< " eigenvalues, u_n= "<<u_n<<", c= "<<c<<endl;

	if(_entropicCorrection)
	{
		*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: entropy scheme not available for IsothermalSinglePhase"<<endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::convectionMatrices: entropy scheme not available for IsothermalSinglePhase");
	}

	/******** Construction des matrices de decentrement ********/
	if( _spaceScheme ==centered){
		if(_usePrimitiveVarsInNewton)//We use primitive variables in Newton iterations
		{
			convectionMatrixPrimitiveVariables(u_n);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
			for(int i=0; i<_nVar*_nVar;i++)
				_absAroeImplicit[i] = 0;
		}
		else//We use conservative variables in Newton iterations
		{
			convectionMatrixConservativeVariables(u_n);//Ici on calcule Acons et on le stocke dans _Aroe
			for(int i=0; i<_nVar*_nVar;i++)
				_absAroe[i] = 0;
		}
	}
	else if(_spaceScheme == upwind )
	{
		if(_Uroe[_nVar]==0.)//infinite sound speed
			throw CdmathException("Upwind scheme cannot be used with incompressible fluids (infinite sound speed->infinite upwinding)");

		/* Calcul de Aprim */
		convectionMatrixPrimitiveVariables(u_n);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
		/* Calcul de Acons (first step in the the computaton of upwinding matrix) */
		convectionMatrixConservativeVariables(u_n);//Ici on calcule Acons et on le stocke dans _Aroe
		/* Calcul de |Acons| */
		vector< complex< double > > y (	nb_vp_dist,0);
		for( int i=0 ; i<nb_vp_dist ; i++)
			y[i] = Polynoms::abs_generalise(vp_dist[i]);
		Polynoms::abs_par_interp_directe( nb_vp_dist, vp_dist, _Aroe, _nVar,_precision, _absAroe,y);//Ici on calcule |Acons| et on le stocke dans _absAroe
		
		/* Calcul de JacU(rho, velocity, 1/c^2) */
		primToConsJacobianMatrix(_Uroe[0], _Uroe+1,_Uroe[_nVar]);//Ici on calcule la jacobienne nabla U en l'état de Roe et on la stocke dans _primConsJacoMat
		//Calcul du produit |Acons|JacU
		Polynoms::matrixProduct(_absAroe, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _absAroeImplicit);
	}
	else if( _spaceScheme ==staggered ){
		if(_entropicCorrection)//To do: study entropic correction for staggered
		{
			*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
			_runLogFile->close();
			throw CdmathException("IsothermalSinglePhase::convectionMatrices: entropy scheme not available for staggered scheme");
		}
		/* Calcul de Aprim */
		convectionMatrixPrimitiveVariables(u_n);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
		/* Calcul du décentrement staggered */
		staggeredRoeUpwindingMatrixPrimitiveVariables( u_n);//Ici on calcule le décentrement staggered et on le stocke dans _absAroeImplicit
	}
	else
	{
		*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: scheme not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::convectionMatrices: scheme not treated");
	}

	if(_usePrimitiveVarsInNewton)//We use primitive variables in Newton iterations
		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinusImplicit[i] = (_AroeImplicit[i]-_absAroeImplicit[i])/2;
			_AroePlusImplicit[i]  = (_AroeImplicit[i]+_absAroeImplicit[i])/2;
		}
	else//We use conservative variables in contribution to right hand side
		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinus[i] = (_Aroe[i]-_absAroe[i])/2;
			_AroePlus[i]  = (_Aroe[i]+_absAroe[i])/2;
			_AroeMinusImplicit[i] = _AroeMinus[i];
			_AroePlusImplicit[i]  = _AroePlus[i];
		}
	
	if(_verbose && _nbTimeStep%_freqSave ==0)
		if(_usePrimitiveVarsInNewton)//We use primitive variables in Newton iterations
		{
			displayMatrix(_AroeImplicit,      _nVar,"Matrice de Roe en variables primitives");
			displayMatrix(_absAroeImplicit,   _nVar,"Décentrement en variables primitives");
			displayMatrix(_AroeMinusImplicit, _nVar,"Matrice _AroeMinus en variables primitives");
			displayMatrix(_AroePlusImplicit,  _nVar,"Matrice _AroePlus en variables primitives");
		}
		else//We use conservative variables in Newton iterations
		{
			displayMatrix(_Aroe,      _nVar,"Matrice de Roe en variables conservatives");
			displayMatrix(_absAroe,   _nVar,"Décentrement en variables conservatives");
			displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus en variables conservatives");
			displayMatrix(_AroePlus,  _nVar,"Matrice _AroePlus en variables conservatives");
		}
}

void IsothermalSinglePhase::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< " face unit normal vector "<<endl;
		for(int k=0; k<_Ndim; k++)
			cout<<normale[k]<<", ";
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne primitif
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< ", inner primitive state Vi = "<<endl;
			for(int k=0; k<_nVar; k++)
				cout<<_externalStates[k]<<", ";
			cout<<endl;
		}
		double u_n=0;//u_n=vitesse normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_n+=_externalStates[(k+1)]*normale[k];
			
		//Pour la convection, inversion du sens de la vitesse normale
		for(int k=0; k<_Ndim; k++)
			_externalStates[(k+1)]-= 2*u_n*normale[k];
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On prend l'état fantôme égal à l'état interne (conditions limites de Neumann)
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne 
		double u_int_n=0;//u_int_n=composante normale de la vitesse à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_int_n+=_externalStates[(k+1)]*normale[k];//On calcule la vitesse normale sortante

		double u_ext_n=_limitField[nameOfGroup].v_x[0]*normale[0];
		if(_Ndim>1)
			{
				u_ext_n+=_limitField[nameOfGroup].v_y[0]*normale[1];
				if(_Ndim>2)
					u_ext_n+=_limitField[nameOfGroup].v_z[0]*normale[2];
			}

		if(u_int_n+u_ext_n<=0)
		{//Interfacial velocity goes inward
		    _externalStates[1] = _limitField[nameOfGroup].v_x[0];
			if(_Ndim>1)
			{
      		   _externalStates[2] = _limitField[nameOfGroup].v_y[0];
                if(_Ndim>2)
			        _externalStates[3] = _limitField[nameOfGroup].v_z[0];
            }			
		}
		else if(_nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid possibly going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état primitif interne

		double u_n=0;//u_n=vitesse normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_n+=_externalStates[(k+1)]*normale[k];
        
		if(u_n<=0)
		{
			//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
			Cell Cj=_mesh.getCell(j);
			double hydroPress=Cj.x()*_GravityField3d[0];
			if(_Ndim>1){
				hydroPress+=Cj.y()*_GravityField3d[1];
				if(_Ndim>2)
					hydroPress+=Cj.z()*_GravityField3d[2];
			}
			hydroPress*= _fluides[0]->getDensity(_limitField[nameOfGroup].p, _Temperature) ;//multiplication by rho the total density
			
			//First component : total pressure
			_externalStates[0] = hydroPress + _limitField[nameOfGroup].p;
			
	        //Contribution from the tangential velocity
	        if(_Ndim>1)
	        {
	            Vector omega(3);
	            omega[0]=_limitField[nameOfGroup].v_x[0];
	            omega[1]=_limitField[nameOfGroup].v_y[0];
	            
	            Vector Normale(3);
	            Normale[0]=normale[0];
	            Normale[1]=normale[1];
	
	            if(_Ndim==3)
	            {
	                omega[2]=_limitField[nameOfGroup].v_z[0];
	                Normale[2]=normale[2];
	            }
	            
	            Vector tangent_vel=omega%Normale;
	
				//Changing external state velocity
	            for(int k=0; k<_Ndim; k++)
	                _externalStates[(k+1)]=u_n*normale[k] + tangent_vel[k];
	        }
		}
		else{
			/*
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature (only pressure value is imposed as in outlet BC)."<<endl;
			_externalStates[0]=_porosityj*(_limitField[nameOfGroup].p+hydroPress);
			*/
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition."<<endl;
		}
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates); //On remplace l'état fantôme par l'état interne PRIMITIF
		double u_n=0;//u_n=vitesse normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
		    u_n+=_externalStates[(k+1)]*normale[k];

		if(u_n < -_precision &&  _nbTimeStep%_freqSave ==0)
		    cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<" with velocity "<< u_n<<endl;
        else
        {
			double hydroPress=0;
			if( _VV.getSpaceDimension()>1)
			{
				//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
				Cell Cj=_mesh.getCell(j);
				double hydroPress=(Cj.x()-_gravityReferencePoint[0])*_GravityField3d[0];
				if(_Ndim>1){
					hydroPress+=(Cj.y()-_gravityReferencePoint[1])*_GravityField3d[1];
					if(_Ndim>2)
						hydroPress+=(Cj.z()-_gravityReferencePoint[2])*_GravityField3d[2];
				}
				hydroPress*= _fluides[0]->getDensity(_limitField[nameOfGroup].p, _Temperature) ;//multiplication by rho the total density
		
				if(_verbose && _nbTimeStep%_freqSave ==0)
				{
					cout<<"Cond lim outlet pressure= "<<_externalStates[0]<<" gravite= "<<_GravityField3d[0]<<" Cj.x()= "<<Cj.x()<<endl;
					cout<<"Cond lim outlet reference pressure= "<<_limitField[nameOfGroup].p<<" pression hydro= "<<hydroPress<<" total= "<<_limitField[nameOfGroup].p+hydroPress<<endl;
				}
			}
			_externalStates[0] = hydroPress + _limitField[nameOfGroup].p;
		}
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<< " _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, InletPressure and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, InletPressure and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}

	for(int k=0; k<_nVar; k++){
		_Vext[k] = _externalStates[k];
		_Vextdiff[k] = _externalStates[k];//Will be changed later for wall boundary condition
	}

	if (_limitField[nameOfGroup].bcType==Wall)	//Pour la diffusion, paroi à vitesse et temperature imposees
	{
		_Vextdiff[0] =_externalStates[0];
		_Vextdiff[1] =_limitField[nameOfGroup].v_x[0];

		if(_Ndim>1)
		{
			_Vextdiff[2] = _limitField[nameOfGroup].v_y[0];
			if(_Ndim==3)
				_Vextdiff[3] = _limitField[nameOfGroup].v_z[0];
		}
	}
	
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< ", ghost primitive state Vj = "<<endl;
		for(int k=0; k<_nVar; k++)
			cout<<_externalStates[k]<<", ";
		cout<<endl;
	}
}

void IsothermalSinglePhase::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)

{
	double mu     = _fluides[0]->getViscosity(_Temperature);

	if(isBord)
	{
		for(int k=0; k<_nVar; k++)
			_Vj[k] = _Vextdiff[k];
		_inv_dxj=_inv_dxi;
	}
	//J'ai remplacé toutes les égalités du type Vj = Vext 
	//car si Vj pointe à la même adresse que Vext : lorsqu'on fait delete[] Vj on fait implicitement delete[] Vext => 
	//donc dès qu'il voit plus tard dans le code delete[] Vext il essaie de désallouer une deuxième fois

	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_nVar; k++)
		_phi[k] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*mu*(_porosityj*_Vj[k] - _porosityi*_Vi[k]);
	
	_idm[0] = i;
	VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion au 2nd membre pour la maille " << i << ": "<<endl;
		for(int q=0; q<_nVar; q++)
			cout << _phi[q] << endl;
		cout << endl;
	}

	if(!isBord)
	{
		//On change de signe pour l'autre contribution
		for(int k=0; k<_nVar; k++)
			_phi[k] *= -_inv_dxj/_inv_dxi;

		_idn[0] = j;
		VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Contribution diffusion au 2nd membre pour la maille  " << j << ": "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}
	}

	if(_verbose && _nbTimeStep%_freqSave ==0 && _timeScheme==Implicit)
	{
		cout << "Matrice de diffusion D, pour le couple (" << i << "," << j<< "):" << endl;
		for(int i=0; i<_nVar; i++)
		{
			for(int j=0; j<_nVar; j++)
				cout << _Diffusion[i*_nVar+j]<<", ";
			cout << endl;
		}
		cout << endl;
	}
}
void IsothermalSinglePhase::addSourceTermToSecondMember
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
			_Vj[k] = _Vext[k];
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

	if(!isBord){
		if(_wellBalancedCorrection){
			for(int k=0; k<_nVar;k++)
				_phi[k]=(_Si[k]+_Sj[k])/2+_pressureLossVector[k]+_porosityGradientSourceVector[k];
			Polynoms::matrixProdVec(_signAroe, _nVar, _nVar, _phi, _l);
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
			Polynoms::matrixProdVec(_signAroe, _nVar, _nVar, _phi, _l);
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
void IsothermalSinglePhase::sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i)
{
	double phirho=Ui[0];
	double norm_u=0;
	for(int k=0; k<_Ndim; k++)
		norm_u+=Vi[1+k]*Vi[1+k];
	norm_u=sqrt(norm_u);

	Si[0]=0;
	for(int k=1; k<_nVar; k++){
		Si[k]  =(_gravite[k]-_dragCoeffs[0]*norm_u*Vi[k])*phirho;
	}

	if(_timeScheme==Implicit)
	{
		for(int k=0; k<_nVar*_nVar;k++)
			_GravityImplicitationMatrix[k] = 0;
		double pression=Vi[0];
		getDensityDerivatives( pression);
		for(int k=0; k<_nVar;k++)
			_GravityImplicitationMatrix[k*_nVar+0] = -_gravite[k]*_drho_sur_dp;
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalSinglePhase::sourceVector"<<endl;
		cout<<"Ui="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Ui[k]<<", ";
		cout<<endl;
		cout<<"Vi="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vi[k]<<", ";
		cout<<endl;
		cout<<"Si="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Si[k]<<", ";
		cout<<endl;
		cout<<"_gravite="<<endl;
		for(int i=0; i<_nVar; i++)
			cout<<_gravite[i]<<", ";
		cout<<endl;
		if(_timeScheme==Implicit)
			displayMatrix(_GravityImplicitationMatrix, _nVar, "Gravity implicitation matrix");
	}
}

void IsothermalSinglePhase::getDensityDerivatives( double pressure)
{
	if(not _isCompressibleFluid)//Case of an incompressible fluid
		_drho_sur_dp = 0;
	else//Case of a compressible fluid
	{
		_drho_sur_dp=_compressibleFluid->getInverseSoundSpeed(pressure, _Temperature);
	
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<"_drho_sur_dp= "<<_drho_sur_dp<<endl;	
	}
}

void IsothermalSinglePhase::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{
	double norm_u=0, u_n=0, rho;
	for(int i=0;i<_Ndim;i++)
		u_n += _Uroe[1+i]*_vec_normal[i];

	pressureLoss[0]=0;
	if(u_n>0){
		for(int i=0;i<_Ndim;i++)
			norm_u += Vi[1+i]*Vi[1+i];
		norm_u=sqrt(norm_u);
		rho=Ui[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-1/2*K*rho*norm_u*Vi[1+i];
	}
	else{
		for(int i=0;i<_Ndim;i++)
			norm_u += Vj[1+i]*Vj[1+i];
		norm_u=sqrt(norm_u);
		rho=Uj[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-1/2*K*rho*norm_u*Vj[1+i];
	}
	pressureLoss[_nVar-1]=-1/2*K*rho*norm_u*norm_u*norm_u;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalSinglePhase::pressureLossVector K= "<<K<<endl;
		cout<<"Ui="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Ui[k]<<", ";
		cout<<endl;
		cout<<"Vi="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vi[k]<<", ";
		cout<<endl;
		cout<<"Uj="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Uj[k]<<", ";
		cout<<endl;
		cout<<"Vj="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Vj[k]<<", ";
		cout<<endl;
		cout<<"pressureLoss="<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<pressureLoss[k]<<", ";
		cout<<endl;
	}
}

void IsothermalSinglePhase::porosityGradientSourceVector()
{
	double u_ni=0, u_nj=0, rhoi,rhoj, pi=_Vi[0], pj=_Vj[0],pij;
	for(int i=0;i<_Ndim;i++) {
		u_ni += _Vi[1+i]*_vec_normal[i];
		u_nj += _Vj[1+i]*_vec_normal[i];
	}
	_porosityGradientSourceVector[0]=0;
	rhoj=_Uj[0]/_porosityj;
	rhoi=_Ui[0]/_porosityi;
	pij=(pi+pj)/2+rhoi*rhoj/2/(rhoi+rhoj)*(u_ni-u_nj)*(u_ni-u_nj);
	for(int i=0;i<_Ndim;i++)
		_porosityGradientSourceVector[1+i]=pij*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
}


void IsothermalSinglePhase::computeScaling(double maxvp)
{
	_blockDiag[0]=1;
	_invBlockDiag[0]=1;//_blockDiag[0];
	for(int q=1; q<_Ndim+1; q++)
	{
		_blockDiag[q]=1/(maxvp*maxvp);
		_invBlockDiag[q]=1;//_blockDiag[q];
		_blockDiag[q+1+_Ndim]=1/(maxvp*maxvp);
		_invBlockDiag[q+1+_Ndim]=1;//_blockDiag[q+1+_Ndim];
	}
}

void IsothermalSinglePhase::jacobian(const int &j, string nameOfGroup,double * normale)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite convection bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_Jcb[k] = 0;//No implicitation at this stage

	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
	double u_n=0;//vitesse normale à la paroi
	for(k=0; k<_Ndim; k++)
		u_n+=_externalStates[(k+1)]*normale[k];

	// switch on boundary types
	if (_limitField[nameOfGroup].bcType==Wall)
	{
		for(k=0; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
		for(k=1; k<_nVar;k++)
			for(int l=1; l<_nVar;l++)
				_Jcb[k*_nVar + l] -= 2*normale[k-1]*normale[l-1];
	}
	else if (_limitField[nameOfGroup].bcType==Inlet)
	{
		if(u_n<0)
		_Jcb[0] = 1;//Only the forst coefficient of the matrix is non zero
		else//Neumann BC
		{
			cout<<"Warning  IsothermalSinglePhase::jacobian : inlet boundary condition requested but fluid going outward. Applying Neumann boundary condition"<<endl;
			for(k=0;k<_nVar;k++)
				_Jcb[k*_nVar+k]=1;//Identity matrix
		}
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure && u_n<0){
		for(k=1; k<_nVar;k++)
			for(int l=1; l<_nVar;l++)
				_Jcb[k*_nVar + l] = normale[k-1]*normale[l-1];
	}
	// not wall, not inlet, not inletPressure
	else if(_limitField[nameOfGroup].bcType==Outlet || (_limitField[nameOfGroup].bcType==InletPressure && u_n>=0))
	{
		for(k=0; k<_nVar;k++)
			for(int l=0; l<_nVar;l++)
				_Jcb[k*_nVar + l] = 0;
		for(k=1; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
	}
	else  if (_limitField[nameOfGroup].bcType!=Neumann)// not wall, not inlet, not outlet
	{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::jacobian: This boundary condition is not treated");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" IsothermalSinglePhase::jacobian : Jacobienne condition limite convection bord "<< nameOfGroup<<endl;
		displayMatrix(_Jcb,_nVar," Jacobian matrix of convection BC : "+nameOfGroup);
	}
}

void IsothermalSinglePhase::jacobianDiff(const int &j, string nameOfGroup)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite diffusion bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;

	if (_limitField[nameOfGroup].bcType==Wall ||_limitField[nameOfGroup].bcType==Inlet ){
		_JcbDiff[0]=1;
	}
	else if ( _limitField[nameOfGroup].bcType==Neumann )
	{
		for(k=0;k<_nVar;k++)
			_JcbDiff[k*_nVar+k]=1;
	}
	else if (_limitField[nameOfGroup].bcType==Outlet || _limitField[nameOfGroup].bcType==InletPressure)
	{
		for(k=1;k<_nVar;k++)
			_JcbDiff[k*_nVar+k]=1;
	}
	else{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition. Known boundary types : Neumann, Wall, Inlet, Outlet, InletPresure." << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition. Known boundary types : Neumann, Wall, Inlet, Outlet, InletPresure." << endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::jacobianDiff: This boundary condition is not recognised. Known boundary types : Neumann, Wall, Inlet, Outlet, InletPresure.");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" IsothermalSinglePhase::jacobianDiff : Jacobienne condition limite diffusion bord "<< nameOfGroup<<endl;
		displayMatrix(_JcbDiff,_nVar," Jacobian matrix of diffusion BC : "+nameOfGroup);
	}
}

void IsothermalSinglePhase::primToCons(const double *P, const int &i, double *W, const int &j)
{   //We do not know the size of Wcons and Wprim 
	//Sometimes they have _nVar components, sometimes they have _Nmailles*_nVar

	assert( P[i*_nVar]>0);//Pressure should be positive

	double phi_rho =_porosityField(j)*_fluides[0]->getDensity(P[i*_nVar], _Temperature);
	W[j*_nVar] =  phi_rho;//phi*rho
	for(int k=0; k<_Ndim; k++)
		W[j*_nVar+(k+1)] = phi_rho*P[i*_nVar+(k+1)];//phi*rho*u

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalSinglePhase::primToCons Vecteur primitif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<P[k]<<endl;
		cout<<"IsothermalSinglePhase::primToCons Vecteur conservatif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<W[k]<<endl;
	}
	assert( W[0]>0);//Density should be positive
	assert( P[0]>0);//Pressure should be positive
}

void IsothermalSinglePhase::primToConsJacobianMatrix(double *V)
{//V vecteur primitif (p,u) de taille _nVar
	double pression=V[0];
	double rho=_fluides[0]->getDensity(pression,_Temperature);
	double invSoundSpeed = _fluides[0]->getInverseSoundSpeed(pression,_Temperature);
	
	primToConsJacobianMatrix( rho, V+1, invSoundSpeed*invSoundSpeed);
}

void IsothermalSinglePhase::primToConsJacobianMatrix(double rho, double* velocity, double invSoundSpeedsquared)
{	
	//Initialise all coeffs to zero
	for(int i=0;i<_nVar*_nVar;i++)
		_primToConsJacoMat[i] = 0;
	//Fill non zero coefficients
	_primToConsJacoMat[0] = invSoundSpeedsquared;//the (0,0) coefficient
	for(int idim=0;idim<_Ndim;idim++)
	{
		_primToConsJacoMat[(idim+1)*_nVar]=velocity[idim]*invSoundSpeedsquared;//the first column
		_primToConsJacoMat[(idim+1)*_nVar+idim+1]=rho;//the diagonal
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" IsothermalSinglePhase::primToConsJacobianMatrix : rho= "<<rho <<", invSoundSpeedsquared="<<invSoundSpeedsquared << endl;
		displayVector( velocity,_Ndim," velocity " );
		cout<<" Jacobian matrix primToCons: " << endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobian matrix primToCons: ");
	}
}

void IsothermalSinglePhase::consToPrim(const double *Wcons, double* Wprim,double porosity)//To do: treat porosity
{   //Function called only with explicit schemes
	//Wcons and Wprim are vectors with _nVar components

	if(not _isCompressibleFluid)
		throw CdmathException("IsothermalSinglePhase::consToPrim should not be used with incompressible fluids");
	else
	{		
		double rho=Wcons[0]/porosity;
		double e =_compressibleFluid->getInternalEnergy(_Temperature, rho);

		Wprim[0] =_compressibleFluid->getPressure(rho*e,rho);//pressure p
		for(int k=1;k<=_Ndim;k++)
			Wprim[k] = Wcons[k]/Wcons[0];//velocity u
	
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"IsothermalSinglePhase::ConsToPrim Vecteur conservatif"<<endl;
			for(int k=0;k<_nVar;k++)
				cout<<Wcons[k]<<endl;
			cout<<"IsothermalSinglePhase::ConsToPrim Vecteur primitif"<<endl;
			for(int k=0;k<_nVar;k++)
				cout<<Wprim[k]<<endl;
		}
		assert( Wcons[0]>0 );//Density should be positive
		assert( Wprim[0]>0 );//Pressure should be positive
		assert( e>0) ;//Internal energy should be positive
	}
}

void IsothermalSinglePhase::addConvectionToSecondMember
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
	VecGetValues(_primitiveVars,    _nVar, _idm, _Vi);

	if(!isBord){
		for(int k=0; k<_nVar; k++)
			_idn[k] = _nVar*j + k;
		VecGetValues(_conservativeVars, _nVar, _idn, _Uj);
		VecGetValues(_primitiveVars,    _nVar, _idn, _Vj);
	}
	else{
		for(int k=0; k<_nVar; k++)
			_Vj[k]= _Vext[k];
		primToCons(_Vj, 0, _Uj, 0);
	}
	_idm[0] = i;
	_idn[0] = j;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "addConvectionToSecondMember : état primitif i= "    << i << endl << " _Vi=" ;
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q] << ", ";
		cout << endl;
		cout << "addConvectionToSecondMember : état conservatif i= " << i << endl << " _Ui=" ;
		for(int q=0; q<_nVar; q++)
			cout << _Ui[q] << ", ";
		cout << endl;
		cout << "addConvectionToSecondMember : état primitif j= "    << j << endl << " _Vj=" ;
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q] <<  ", ";
		cout << endl;
		cout << "addConvectionToSecondMember : état conservatif j= " << j << endl << " _Uj=" ;
		for(int q=0; q<_nVar; q++)
			cout << _Uj[q] <<  ", ";
		cout << endl;
	}
	
	if(_usePrimitiveVarsInNewton)//We use primitive variables in Newton iterations
	{
		for(int k=0; k<_nVar; k++)
			_temp[k]=(_Vi[k] - _Vj[k])*_inv_dxi;//(Vi-Vj)*_inv_dxi
		Polynoms::matrixProdVec(_AroeMinusImplicit, _nVar, _nVar, _temp, _phi);//phi=A^-(V_i-V_j)/dx
		VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Ajout convection au 2nd membre pour les etats " << i << "," << j << endl;
			cout << "(Vi - Vj)*_inv_dxi= "<<endl;;
			for(int q=0; q<_nVar; q++)
				cout << _temp[q] << endl;
			cout << endl;
			cout << "Contribution convection à " << i << endl;
			cout << "A^-*(Vi - Vj)*_inv_dxi= "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}

		if(!isBord)
		{
			for(int k=0; k<_nVar; k++)
				_temp[k]*=_inv_dxj/_inv_dxi;//(Vi-Vj)*_inv_dxi
			Polynoms::matrixProdVec(_AroePlusImplicit, _nVar, _nVar, _temp, _phi);//phi=A^+(V_i-V_j)/dx
			VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);
	
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution convection à  " << j << endl;
				cout << "A^+*(Vi - Vj)*_inv_dxi= "<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _phi[q] << endl;
				cout << endl;
			}
		}
	}
	else//We use conservative variables in Newton iterations
	{
		for(int k=0; k<_nVar; k++)
			_temp[k]=(_Ui[k] - _Uj[k])*_inv_dxj;//(Ui-Uj)*_inv_dxj
		Polynoms::matrixProdVec(_AroeMinus, _nVar, _nVar, _temp, _phi);//phi=A^-(U_i-U_j)/dx
		VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout << "Ajout convection au 2nd membre pour les etats " << i << "," << j << endl;
			cout << "(Ui - Uj)*_inv_dxj= "<<endl;;
			for(int q=0; q<_nVar; q++)
				cout << _temp[q] << endl;
			cout << endl;
			cout << "Contribution convection à " << i << endl;
			cout << "A^-*(Ui - Uj)*_inv_dxj= "<<endl;
			for(int q=0; q<_nVar; q++)
				cout << _phi[q] << endl;
			cout << endl;
		}

		if(!isBord)
		{
			for(int k=0; k<_nVar; k++)
				_temp[k]*=_inv_dxj/_inv_dxi;//(Ui-Uj)*_inv_dxj
			Polynoms::matrixProdVec(_AroePlus, _nVar, _nVar, _temp, _phi);//phi=A^+(U_i-U_j)/dx
			VecSetValuesBlocked(_b, 1, _idn, _phi, ADD_VALUES);
	
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "Contribution convection à  " << j << endl;
				cout << "A^+*(Ui - Uj)*_inv_dxj= "<<endl;
				for(int q=0; q<_nVar; q++)
					cout << _phi[q] << endl;
				cout << endl;
			}
		}
	}
}

void IsothermalSinglePhase::convectionMatrixConservativeVariables(double u_n)
{
	assert( abs(_Uroe[_nVar])>0);
	double c_2 = 1./_Uroe[_nVar];
	
	/******** Construction de la matrice de Roe *********/
	//premiere ligne (masse)
	_Aroe[0]=0;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[1+idim]=_vec_normal[idim];

	//lignes intermédiaires(qdm)
	for(int idim=0; idim<_Ndim;idim++)
	{
		//premiere colonne
		_Aroe[(1+idim)*_nVar]=c_2*_vec_normal[idim] - u_n*_Uroe[1+idim];
		//colonnes intermediaires
		for(int jdim=0; jdim<_Ndim;jdim++)
			_Aroe[(1+idim)*_nVar + jdim + 1] = _Uroe[1+idim]*_vec_normal[jdim];
		//matrice identite
		_Aroe[(1+idim)*_nVar + idim + 1] += u_n;
	}
}

void IsothermalSinglePhase::convectionMatrixPrimitiveVariables(double u_n )
{
	//On remplit la matrice de Roe en variables primitives : F(V_L)-F(V_R)=Aroe (V_L-V_R)
	double rho = _Uroe[0 ];
	_drho_sur_dp = _Uroe[_nVar];
	
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]= rho*_vec_normal[i];
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*_Uroe[1+i]+_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*_Uroe[1+i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
	}
}

void IsothermalSinglePhase::staggeredRoeUpwindingMatrixPrimitiveVariables( double u_n)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//Calcul de décentrement de type décalé pour formulation Roe
	double rho = _Uroe[ 0 ];
	_drho_sur_dp = _Uroe[_nVar];
	
	_absAroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_absAroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	for(int i=0;i<_Ndim;i++)
	{
		_absAroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*_Uroe[1+i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_absAroeImplicit[(1+i)*_nVar+1+j]=rho*_Uroe[1+i]*_vec_normal[j];
		_absAroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
	}

	double signu=0;
	if(u_n>_precision)
		signu=1;
	else if (u_n<-_precision)
		signu=-1;

	for(int i=0; i<_nVar*_nVar;i++)
		_absAroeImplicit[i] *= signu;//atan(100*u_n);

}

Vector IsothermalSinglePhase::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalSinglePhase::convectionFlux start"<<endl;
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

	Vector F(_nVar);
	F(0)=phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F(1+i)=phirho*vitessen*vitesse(i)+pression*normale(i)*porosity;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalSinglePhase::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void IsothermalSinglePhase::applyVFRoeLowMachCorrections(bool isBord, string groupname)
{
	if(_nonLinearFormulation!=VFRoe)
	{
		*_runLogFile<< "IsothermalSinglePhase::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation" << endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation");
	}
	else//_nonLinearFormulation==VFRoe
	{
		if(_spaceScheme==lowMach){
			double u_2=0;
			for(int i=0;i<_Ndim;i++)
				u_2 += _Uroe[1+i]*_Uroe[1+i];
			double 	c = _maxvploc;//vitesse du son a l'interface
			double M=sqrt(u_2)/c;//Mach number
			_Vij[0]=M*_Vij[0]+(1-M)*(_Vi[0]+_Vj[0])/2;
		}
		else if(_spaceScheme==pressureCorrection)
		{//order 1 : no correction, oarder 2 : correction everywhere, order 3 : correction only inside, orders 4 and 5 : special correction at boundaries
			if(_pressureCorrectionOrder==2 || (!isBord && _pressureCorrectionOrder==3) || (!isBord && _pressureCorrectionOrder==4) || (!isBord && _pressureCorrectionOrder==5) )
			{
				double norm_uij=0, uij_n=0, ui_n=0, uj_n=0;
				for(int i=0;i<_Ndim;i++)
				{
					norm_uij += _Uroe[1+i]*_Uroe[1+i];
					uij_n += _Uroe[1+i]*_vec_normal[i];
					ui_n += _Vi[1+i]*_vec_normal[i];
					uj_n += _Vj[1+i]*_vec_normal[i];
				}
				norm_uij=sqrt(norm_uij);
				if(norm_uij>_precision)//avoid division by zero
					_Vij[0]=(_Vi[0]+_Vj[0])/2 + uij_n/norm_uij*(_Vj[0]-_Vi[0])/4 - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
				else
					_Vij[0]=(_Vi[0]+_Vj[0])/2                                    - _Uroe[0]*norm_uij*(uj_n-ui_n)/4;
			}
			else if(_pressureCorrectionOrder==4 && isBord)
				_Vij[0]=_Vi[0];
			else if(_pressureCorrectionOrder==5 && isBord)
			{
				double g_n=0;//scalar product of gravity and normal vector
				for(int i=0;i<_Ndim;i++)
					g_n += _GravityField3d[i]*_vec_normal[i];
				_Vij[0]=_Vi[0]- _Ui[0]*g_n/_inv_dxi/2;
			}
		}
		else if(_spaceScheme==staggered)
		{
			double uij_n=0;
			for(int i=0;i<_Ndim;i++)
				uij_n += _Uroe[1+i]*_vec_normal[i];

			if(uij_n>_precision){
				_Vij[0]=_Vj[0];
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=_Vi[1+i];
				_Vij[_nVar-1]=_Vi[_nVar-1];
			}
			else if(uij_n<-_precision){
				_Vij[0]=_Vi[0];
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=_Vj[1+i];
				_Vij[_nVar-1]=_Vj[_nVar-1];
			}
			else{
				_Vij[0]=(_Vi[0]+_Vi[0])/2;
				for(int i=0;i<_Ndim;i++)
					_Vij[1+i]=(_Vj[1+i]+_Vj[1+i])/2;
				_Vij[_nVar-1]=(_Vj[_nVar-1]+_Vj[_nVar-1])/2;
			}
		}
		primToCons(_Vij,0,_Uij,0);
	}
}

void IsothermalSinglePhase::testConservation()
{
	double SUM, DELTA, x;
	int I;
	for(int i=0; i<_nVar; i++)
	{
		{
			if(i == 0)
				cout << "Masse totale (kg): ";
			else
			{
				if(i == _nVar-1)
					cout << "Energie totale (J): ";
				else
					cout << "Quantite de mouvement totale (kg.m.s^-1): ";
			}
		}
		SUM = 0;
		I =  i;
		DELTA = 0;
		for(int j=0; j<_Nmailles; j++)
		{
			if(!_usePrimitiveVarsInNewton)
				VecGetValues(_conservativeVars, 1, &I, &x);//on recupere la valeur du champ
			else
				VecGetValues(_primitiveVars, 1, &I, &x);//on recupere la valeur du champ
			SUM += x*_mesh.getCell(j).getMeasure();
			VecGetValues(_newtonVariation, 1, &I, &x);//on recupere la variation du champ
			DELTA += x*_mesh.getCell(j).getMeasure();
			I += _nVar;
		}
		if(fabs(SUM)>_precision)
			cout << SUM << ", variation relative: " << fabs(DELTA /SUM)  << endl;
		else
			cout << " a une somme quasi nulle,  variation absolue: " << fabs(DELTA) << endl;
	}
}

void IsothermalSinglePhase::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/IsothermalSinglePhasePrim_");///Results
	string cons(_path+"/IsothermalSinglePhaseCons_");
	string allFields(_path+"/");
	prim+=_fileName;
	cons+=_fileName;
	allFields+=_fileName;

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

	if (_nbTimeStep ==0 || _restartWithNewFileName){	// write mesh and component info	
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
	else{// do not write mesh and component info
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
	if(_saveVelocity || _saveAllFields){
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
		if (_nbTimeStep ==0 || _restartWithNewFileName){// write mesh and component info		
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

	if(_saveAllFields)
	{
		double p,T,rho, h, vx,vy,vz,v2;
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&rho);
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar + _nVar - 1;
			VecGetValues(_primitiveVars,1,&Ii,&T);
			Ii = i*_nVar + 1;
			VecGetValues(_primitiveVars,1,&Ii,&vx);
			if(_Ndim>1)
			{
				Ii = i*_nVar + 2;
				VecGetValues(_primitiveVars,1,&Ii,&vy);
				if(_Ndim>2){
					Ii = i*_nVar + 3;
					VecGetValues(_primitiveVars,1,&Ii,&vz);
				}
			}

			_Density(i)=rho;
			_Pressure(i)=p;
			_VitesseX(i)=vx;
			v2=vx*vx;
			if(_Ndim>1)
			{
				_VitesseY(i)=vy;
				v2+=vy*vy;
				if(_Ndim>2)
				{
					_VitesseZ(i)=vz;
					v2+=vz*vz;
				}
			}
			_MachNumber(i)=sqrt(v2)*_fluides[0]->getInverseSoundSpeed(p,T);
		}
		_Density.setTime(_time,_nbTimeStep);
		_Pressure.setTime(_time,_nbTimeStep);
		_MachNumber.setTime(_time,_nbTimeStep);
		_VitesseX.setTime(_time,_nbTimeStep);
		if(_Ndim>1)
		{
			_VitesseY.setTime(_time,_nbTimeStep);
			if(_Ndim>2)
				_VitesseZ.setTime(_time,_nbTimeStep);
		}
		if (_nbTimeStep ==0 || _restartWithNewFileName){// write mesh and component info		
			switch(_saveFormat)
			{
			case VTK :
				_Density.writeVTK(allFields+"_Density");
				_Pressure.writeVTK(allFields+"_Pressure");
				_MachNumber.writeVTK(allFields+"_MachNumber");
				_VitesseX.writeVTK(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ");
				}
				break;
			case MED :
				_Density.writeMED(allFields+"_Density");
				_Pressure.writeMED(allFields+"_Pressure");
				_MachNumber.writeMED(allFields+"_MachNumber");
				_VitesseX.writeMED(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ");
				}
				break;
			case CSV :
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_MachNumber.writeCSV(allFields+"_MachNumber");
				_VitesseX.writeCSV(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeCSV(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeCSV(allFields+"_VelocityZ");
				}
				break;
			}
		}
		else{// do not write mesh and component info
			switch(_saveFormat)
			{
			case VTK :
				_Density.writeVTK(allFields+"_Density",false);
				_Pressure.writeVTK(allFields+"_Pressure",false);
				_MachNumber.writeVTK(allFields+"_MachNumber",false);
				_VitesseX.writeVTK(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeVTK(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeVTK(allFields+"_VelocityZ",false);
				}
				break;
			case MED :
				_Density.writeMED(allFields+"_Density",false);
				_Pressure.writeMED(allFields+"_Pressure",false);
				_MachNumber.writeMED(allFields+"_MachNumber",false);
				_VitesseX.writeMED(allFields+"_VelocityX",false);
				if(_Ndim>1)
				{
					_VitesseY.writeMED(allFields+"_VelocityY",false);
					if(_Ndim>2)
						_VitesseZ.writeMED(allFields+"_VelocityZ",false);
				}
				break;
			case CSV :
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_MachNumber.writeCSV(allFields+"_MachNumber");
				_VitesseX.writeCSV(allFields+"_VelocityX");
				if(_Ndim>1)
				{
					_VitesseY.writeCSV(allFields+"_VelocityY");
					if(_Ndim>2)
						_VitesseZ.writeCSV(allFields+"_VelocityZ");
				}
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

		if(_saveVelocity || _saveAllFields){
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

	if (_nbTimeStep ==0 || _restartWithNewFileName){	// delete mesh in memory	
		_VV.deleteMEDCouplingMesh();
		_UU.deleteMEDCouplingMesh();
	}
	if (_restartWithNewFileName)
		_restartWithNewFileName=false;
}

Vector IsothermalSinglePhase::staggeredVFFCFlux()
{
	*_runLogFile<< "IsothermalSinglePhase::staggeredVFFCFlux is not yet available "<<  endl;
	_runLogFile->close();
	throw CdmathException("IsothermalSinglePhase::staggeredVFFCFlux is not yet available");
}

void IsothermalSinglePhase::entropicShift(double* n)
{
	*_runLogFile<< "IsothermalSinglePhase::entropicShift is not yet available "<<  endl;
	_runLogFile->close();
	throw CdmathException("IsothermalSinglePhase::entropicShift is not yet available");
}

vector<string> IsothermalSinglePhase::getOutputFieldsNames()
{
	vector<string> result(7);
	
	result[0]="Pressure";
	result[1]="Velocity";
	result[2]="Density";
	result[3]="Momentum";
	result[4]="VelocityX";
	result[5]="VelocityY";
	result[6]="VelocityZ";
	
	return result;
}

Field& IsothermalSinglePhase::getOutputField(const string& nameField )
{
	string nameField_lower_case=nameField;
	transform(nameField_lower_case.begin(), nameField_lower_case.end(), nameField_lower_case.begin(), ::tolower);
	
	if(nameField_lower_case=="pressure" || nameField_lower_case=="pression" )
		return getPressureField();
	else if(nameField_lower_case=="velocity" || nameField_lower_case=="vitesse" )
		return getVelocityField();
	else if(nameField_lower_case== "velocityx" || nameField_lower_case=="vitessex" )
		return getVelocityXField();
	else if(nameField_lower_case=="velocityy" || nameField_lower_case=="vitessey" )
		return getVelocityYField();
	else if(nameField_lower_case=="velocityz" || nameField_lower_case=="vitessez" )
		return getVelocityZField();
	else if(nameField_lower_case=="density" || nameField_lower_case=="densite" )
		return getDensityField();
	else if(nameField_lower_case=="momentum" || nameField_lower_case=="qdm" )
		return getMomentumField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first to check" << endl;
        throw CdmathException("IsothermalSinglePhase::getOutputField error : Unknown Field name");
    }
}

Field& IsothermalSinglePhase::getPressureField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getPressureField. Call initialize first");

	if(!_saveAllFields)
	{
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&_Pressure(i));
		}
		_Pressure.setTime(_time,_nbTimeStep);
		_Pressure.setInfoOnComponent(0,"Pressure_(Pa)");
	}
	return _Pressure;
}

Field& IsothermalSinglePhase::getVelocityField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getVelocityField. Call initialize first");

	if(!_saveAllFields )
	{
		_Vitesse=Field("Vitesse",CELLS,_mesh,3);
		int Ii;
		for (long i = 0; i < _Nmailles; i++)
		{
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +1+j;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse(i,j));
			}
			for (int j = _Ndim; j < 3; j++)//On met à zero les composantes de vitesse si la dimension est <3
				_Vitesse(i,j)=0;
		}
		_Vitesse.setTime(_time,_nbTimeStep);
		_Vitesse.setInfoOnComponent(0,"Velocity_x_(m/s)");
		_Vitesse.setInfoOnComponent(1,"Velocity_y_(m/s)");
		_Vitesse.setInfoOnComponent(2,"Velocity_z_(m/s)");
	}
	
	return _Vitesse;
}

Field& IsothermalSinglePhase::getMachNumberField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getMachNumberField. Call initialize first");

	if(!_saveAllFields )
	{
		_MachNumber=Field("Mach number",CELLS,_mesh,1);
		int Ii;
		double p,T,rho,h, temp, u2=0;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&T);
			
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				int Ii = i*_nVar +1+j;
				VecGetValues(_primitiveVars,1,&Ii,&temp);
				u2+=temp*temp;
			}
	
			_MachNumber[i]  =sqrt(u2)*_fluides[0]->getInverseSoundSpeed(p,T);
			//cout<<"u="<<sqrt(u2)<<", 1/c= "<<_fluides[0]->getInverseSoundSpeed(p,T);<<", MachNumberField[i] = "<<MachNumberField[i] <<endl;
		}
		_MachNumber.setTime(_time,_nbTimeStep);
	}
	//cout<<", MachNumberField = "<<MachNumberField <<endl;

	return _MachNumber;
}

Field& IsothermalSinglePhase::getVelocityXField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getVelocityXField. Call initialize first");

	if(!_saveAllFields )
	{
		_VitesseX=Field("Velocity X",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++)
		{
			int Ii = i*_nVar +1;
			VecGetValues(_primitiveVars,1,&Ii,&_VitesseX(i));
		}
		_VitesseX.setTime(_time,_nbTimeStep);
		_VitesseX.setInfoOnComponent(0,"Velocity_x_(m/s)");
	}
	
	return _VitesseX;
}

Field& IsothermalSinglePhase::getVelocityYField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getVelocityYField. Call initialize first");

	if(_Ndim<2)
        throw CdmathException("IsothermalSinglePhase::getVelocityYField() error : dimension should be at least 2");	
	else
		if(!_saveAllFields )
		{
			_VitesseY=Field("Velocity Y",CELLS,_mesh,1);
			int Ii;
			for (long i = 0; i < _Nmailles; i++)
			{
				int Ii = i*_nVar +2;
				VecGetValues(_primitiveVars,1,&Ii,&_VitesseY(i));
			}
			_VitesseY.setTime(_time,_nbTimeStep);
			_VitesseY.setInfoOnComponent(0,"Velocity_y_(m/s)");
		}
		
		return _VitesseY;
}

Field& IsothermalSinglePhase::getVelocityZField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getVelocityZField. Call initialize first");

	if(_Ndim<3)
        throw CdmathException("IsothermalSinglePhase::getvelocityZField() error : dimension should be 3");	
	else
		if(!_saveAllFields )
		{
			_VitesseZ=Field("Velocity Z",CELLS,_mesh,1);
			int Ii;
			for (long i = 0; i < _Nmailles; i++)
			{
				int Ii = i*_nVar +3;
				VecGetValues(_primitiveVars,1,&Ii,&_VitesseZ(i));
			}
			_VitesseZ.setTime(_time,_nbTimeStep);
			_VitesseZ.setInfoOnComponent(0,"Velocity_z_(m/s)");
		}
		
		return _VitesseZ;
}

Field& IsothermalSinglePhase::getDensityField()
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getDensityField. Call initialize first");
		
	if(!_saveAllFields )
	{
		_Density=Field("Density",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&_Density(i));
		}
		_Density.setTime(_time,_nbTimeStep);
		_Density.setInfoOnComponent(0,"Density_(kg/m^3)");
	}
	return _Density;
}

Field& IsothermalSinglePhase::getMomentumField()//not yet managed by parameter _saveAllFields
{
	if(!_initializedMemory)
		throw CdmathException("IsothermalSinglePhase::getMomentumField. Call initialize first");

	_Momentum=Field("Momentum",CELLS,_mesh,_Ndim);
	int Ii;
	for (long i = 0; i < _Nmailles; i++)
		for (int j = 0; j < _Ndim; j++)//On récupère les composantes de qdm
		{
			int Ii = i*_nVar +1+j;
			VecGetValues(_conservativeVars,1,&Ii,&_Momentum(i,j));
		}
	_Momentum.setTime(_time,_nbTimeStep);
	_Momentum.setInfoOnComponent(0,"Momentum_x");// (kg/m^2/s)
	if (_Ndim>1)
		_Momentum.setInfoOnComponent(1,"Momentum_y");// (kg/m^2/s)
	if (_Ndim>2)
		_Momentum.setInfoOnComponent(2,"Momentum_z");// (kg/m^2/s)

	return _Momentum;
}
