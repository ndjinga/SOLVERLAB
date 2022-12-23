/*
 * SinglePhase.cxx
 *
 *  Created on: Sep 15, 2014
 *      Author: tn236279
 */

#include "NavierStokes.hxx"

using namespace std;

NavierStokes::NavierStokes(phaseType fluid, pressureEstimate pEstimate, int dim, bool isCompressibleFluid)){
	_Ndim=dim;
	_nVar=_Ndim+2;
	_nbPhases  = 1;
	_dragCoeffs=vector<double>(1,0);
	_useDellacherieEOS=useDellacherieEOS;
	_saveAllFields=false;

	if(pEstimate==around1bar300K){//EOS at 1 bar and 300K
		_Tref=300;
		_Pref=1e5;
		if(fluid==Gas){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_compressibleFluid = new StiffenedGas(1.4,743,_Tref,2.22e5);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, e=2.22e5, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_compressibleFluid = new StiffenedGas(996,_Pref,_Tref,1.12e5,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C, e=1.12e5, c_v=4130
		}
	}
	else{//EOS at 155 bars and 618K 
		_Tref=618;
		_Pref=155e5;
		if(fluid==Gas){
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_compressibleFluid = new StiffenedGas(102,_Pref,_Tref,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C
		}
		else{
			cout<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			if(_useDellacherieEOS)
				_compressibleFluid= new StiffenedGasDellacherie(2.35,1e9,-1.167e6,1816); //stiffened gas law for water from S. Dellacherie
			else
				_compressibleFluid= new StiffenedGas(726.82,_Pref,_Tref,1.3e6, 971.,5454.);  //stiffened gas law for water at pressure 155 bar, and temperature 345°C
		}
	}

	//Save into the fluid list
	_fluides.resize(1,_compressibleFluid);

	_fileName = "SolverlabSinglePhase";
    PetscPrintf(PETSC_COMM_WORLD,"\n Navier-Stokes equations for single phase flow\n");
}

void NavierStokes::initialize(){
	cout<<"\n Initialising the Navier-Stokes model\n"<<endl;
	*_runLogFile<<"\n Initialising the Navier-Stokes model\n"<<endl;

	_Uroe = new double[_nVar];//Deleted in ProblemFluid::terminate()
	_Vext, _Vdiff,_Vextdiff = new double[_nVar];

	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+1]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];//Deleted in ProblemFluid::terminate()

	if(_saveVelocity || _saveAllFields)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant

	if(_saveAllFields)
	{
		_Enthalpy=Field("Enthalpy",CELLS,_mesh,1);
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		_Pressure.setInfoOnComponent(0,"Pressure_(Pa)");
		_Density=Field("Density",CELLS,_mesh,1);
		_Density.setInfoOnComponent(0,"Density_(kg/m^3)");
		_Temperature=Field("Temperature",CELLS,_mesh,1);
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

	//* In case of an incompressible fluid, pressure may be defined up to a constant hence a singular inear system (no uniqueness unless outlet boundary condition) */
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

void Navierstokes::terminate(){
	delete[] _Vdiff,_Vextdiff,_Vext;
	
	if( _isSingularSystem )
		VecDestroy(&_constantPressureVector);

	ProblemFluid::terminate();
}

void NavierStokes::convectionState( const long &i, const long &j, const bool &IsBord){
	// Computation of "Roe values" rho, u, H
	//Extraction of primitive and conservative states
	_idm[0] = _nVar*i; 
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	if(IsBord)
		VecGetValues(_Uext, _nVar, _idm, _Uj);
		VecGetValues(_Vext, _nVar, _idm, _Vj);
	else
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection CONSERVATIVE Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Ui[k]<<endl;
		cout<<"Convection CONSERVATIVE Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Uj[k]<<endl;
		cout<<"Convection PRIMITIVE Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Vi[k]<<endl;
		cout<<"Convection PRIMITIVE Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Vj[k]<<endl;
	}
	assert(_Vj[0]>0);
	assert(_Vi[0]>0);

	//Computation of Roe state
	PetscScalar ri, rj, xi, xj, H_i, H_j, kinetic_i, kinetic_j;
	ri = sqrt(_Ui[0]);//racine carre de phi_i rho_i
	rj = sqrt(_Uj[0]);
	_Uroe[0] = ri*rj;	//moyenne geometrique des densites
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Densité moyenne Roe  gauche " << i << ": " << ri*ri << ", droite " << j << ": " << rj*rj << "->" << _Uroe[0] << endl;
	for(int k=0;k<_Ndim;k++){
		kinetic_i += (_Vi[ k+1 ] * _Vi[ k+1 ]);
		kinetic_j += (_Vj[ k+1 ] * _Vj[ k+1 ]);
		xi = _Ui[k+1];
		xj = _Uj[k+1];
		_Uroe[1+k] = (xi/ri + xj/rj)/(ri + rj);
		//"moyenne" des vitesses
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout << "Vitesse de Roe composante "<< k<<"  gauche " << i << ": " << xi/(ri*ri) << ", droite " << j << ": " << xj/(rj*rj) << "->" << _Uroe[k+1] << endl;
	}
	H_i = _Vi[ _nVar-1 ] + kinetic_i/2.0;
	H_j = _Vj[ _nVar-1 ] + kinetic_j/2.0;
	_Uroe[_nVar - 1] = (ri*H_i + H_j*rj)/(ri + rj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Enthalpie totale de Roe H  gauche " << i << ": " << H_i << ", droite " << j << ": " << H_j << "->" << _Uroe[_nVar-1] << endl;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
	
}

void NavierStokes::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat de diffusion (p, u,  h)
	__idm[0] = _nVar*i;
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
		cout << "NavierStokes::diffusionPrimitiveStateAndMatrices cellule gauche = " << i << endl;
		cout << "Vi = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q]  << "\t";
		cout << endl;
		cout << "NavierStokes::diffusionPrimitiveStateAndMatrices cellule droite =" << j << endl;
		cout << "Vj = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q]  << "\t";
		cout << endl;
	}

	for(int k=0; k<_nVar; k++)
		_Vdiff[k] = (_Vi[k]/_porosityi+_Vj[k]/_porosityj)/2;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "NavierStokes::diffusionPrimitiveStateAndMatrices primitive diffusion state" << endl;
		cout << "_Vdiff = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vdiff[q]  << "\t";
		cout << endl;
		cout << "porosite gauche= "<<_porosityi<< ", porosite droite= "<<_porosityj<<endl;
	}

	if(_timeScheme==Implicit)
	{
		double rho = _fluides[0]->getDensityFromEnthalpy(_Vdiff[0],_Vdiff[ _nVar-1 ] );
		double T = _fluides[0]->getTemperatureFromEnthalpy( _Vdiff[ _nVar -1], rho );
		
		double lambda = _fluides[0]->getConductivity( T, _Vdiff[0] );
		double mu     = _fluides[0]->getViscosity( T, rho  );
		_fluides[0]->getConductivity( T, _Vdiff[0] );
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		for(int i=1;i<_nVar-1;i++)
			_Diffusion[i*_nVar+i] = mu;
		//TODO à définir
		double dT_dp, dT_dh  = getTemperatureDerivatives(p,h);
		_Diffusion[(_nVar-1)*_nVar]= - lambda * dT_dp;
		_Diffusion[_nVar*_nVar-1]= - lambda * dT_dh;
	}
}

void NavierStokes::setBoundaryState(string nameOfGroup, const int &j,double *normale){
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
			_externalStates[_nVar - 1] = _limitField[nameOfGroup].h;
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

void NavierStokes::convectionMatrices()
{
	//entree: URoe = [rho, u, H],
	//sortie: matrices Roe+  et Roe-

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"NavierStokes::convectionMatrices()"<<endl;

	double u_n=0, u_2=0;//vitesse normale et carré du module

	for(int idim=0; idim <_Ndim; idim++)
	{
		u_2 += _Uroe[1+idim]*_Uroe[1+idim];
		u_n += _Uroe[1+idim]*_vec_normal[idim];
		vitesse[idim]=_Uroe[1+idim];
	}
	vector<std::complex<double> > vp_dist(3,0);

	double  invSound, H, c;
	/***********Calcul des valeurs propres ********/
	H = _Uroe[_nVar-1 ];
	invSound = 1/_compressibleFluid->vitesseSonEnthalpie(H-u_2/2).0;
	c = 1/sqrt(invSound).0;

	vp_dist[0]=u_n-c;vp_dist[1]=u_n;vp_dist[2]=u_n+c;

	if(invSound==0.){//infinite sound speed
		if(_timeScheme==Explicit)
			throw CdmathException("Explicit scheme cannot be used for incompressible fluids since dt=0");
		else
		{
			nb_vp_dist=1;
			vp_dist.resize(nb_vp_dist);
			c=0.;//The velocity will be used to determine the time step
			vp_dist[0]=u_n;//store the fluid velocity for calculation of the time step
		}
	}
	else
	{
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
		cout<<"NavierStokes::convectionMatrices, "<< nb_vp_dist<< " eigenvalues, u_n= "<<u_n<<", c= "<<c<<endl;

	if(_entropicCorrection)
	{
		*_runLogFile<<"NavierStokes::convectionMatrices: entropy scheme not available for NavierStokes"<<endl;
		_runLogFile->close();
		throw CdmathException("NavierStokes::convectionMatrices: entropy scheme not available for NavierStokes");
	}

	/******** Construction des matrices de decentrement ********/
	if( _spaceScheme ==centered){
		if(_usePrimitiveVarsInNewton)//We use primitive variables in Newton iterations
		{
			convectionMatrixPrimitiveVariables(_Uroe[0], u_n, _Uroe[2], vitesse);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
			for(int i=0; i<_nVar*_nVar;i++)
				_absAroeImplicit[i] = 0;
		}
		else//We use conservative variables in Newton iterations
		{
			//TODO calculer k, K pour une loi d'état générale
			k = _compressibleFluid->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
			K = u_2*k/2; //g-1/2 *|u|²
			convectionMatrixConservativeVariables(u_n, H, vitesse, k, K );//Ici on calcule Acons et on le stocke dans _Aroe
			for(int i=0; i<_nVar*_nVar;i++)
				_absAroe[i] = 0;
		}
	}
	else if(_spaceScheme == upwind )
	{
		if( invSound ==0.)//infinite sound speed
			throw CdmathException("Upwind scheme cannot be used with incompressible fluids (infinite sound speed->infinite upwinding)");

		/* Calcul de Aprim */
		convectionMatrixPrimitiveVariables(_Uroe[0], u_n, _Uroe[2], vitesse);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
		/* Calcul de Acons (first step in the the computaton of upwinding matrix) */
		convectionMatrixConservativeVariables(_Uroe[0], u_n, _Uroe[2], vitesse);//Ici on calcule Acons et on le stocke dans _Aroe
		/* Calcul de |Acons| */
		vector< complex< double > > y (	nb_vp_dist,0);
		for( int i=0 ; i<nb_vp_dist ; i++)
			y[i] = Polynoms::abs_generalise(vp_dist[i]);
		Polynoms::abs_par_interp_directe( nb_vp_dist, vp_dist, _Aroe, _nVar,_precision, _absAroe,y);//Ici on calcule |Acons| et on le stocke dans _absAroe
		
		//Calcul de JacU
		primToConsJacobianMatrix(_Uroe[0], _Uroe+1,_Uroe[_nVar]);//Ici on calcule la jacobienne nabla U en l'état de Roe et on la stocke dans _primConsJacoMat
		//Calcul du produit |Acons|JacU
		Polynoms::matrixProduct(_absAroe, _nVar, _nVar, _primToConsJacoMat, _nVar, _nVar, _absAroeImplicit);
	}
	else if( _spaceScheme ==staggered ){
		if(_entropicCorrection)//To do: study entropic correction for staggered
		{
			*_runLogFile<<"NavierStokes::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
			_runLogFile->close();
			throw CdmathException("NavierStokes::convectionMatrices: entropy scheme not available for staggered scheme");
		}
		/* Calcul de Aprim */
		convectionMatrixPrimitiveVariables(_Uroe[0], u_n, _Uroe[2],vitesse);//Ici on calcule Aprim et on le stocke dans _AroeImplicit
		/* Calcul du décentrement staggered */
		staggeredRoeUpwindingMatrixPrimitiveVariables( u_n);//Ici on calcule le décentrement staggered et on le stocke dans _absAroeImplicit
		
	}
	else
	{
		*_runLogFile<<"NavierStokes::convectionMatrices: scheme not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("NavierStokes::convectionMatrices: scheme not treated");
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

void NavierStokes::computeScaling(double maxvp)
{
	_blockDiag[0]=1;
	_invBlockDiag[0]=1;
	for(int q=1; q<_nVar-1; q++)
	{
		_blockDiag[q]=1./maxvp;//
		_invBlockDiag[q]= maxvp;//1.;//
	}
	_blockDiag[_nVar - 1]=(_compressibleFluid->constante("gamma")-1)/(maxvp*maxvp);//1
	_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
}


void IsothermalSinglePhase::addConvectionToSecondMember(const int &i, const int &j, bool isBord, string groupname)
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

void NavierStokes::addDiffusionToSecondMember(const int &i, const int &j, bool isBord)
{
	double rho = _fluides[0]->getDensityFromEnthalpy( _Vdiff[0], _Vdiff[ _nVar-1 ] );
	double T = _fluides[0]->getTemperatureFromEnthalpy( _Vdiff[ _nVar -1], rho );
	
	double lambda = _fluides[0]->getConductivity( T, _Vdiff[0] );
	double mu     = _fluides[0]-getViscosity( T, rho  );

	if(isBord )
		lambda=max(lambda,_heatTransfertCoeff);//wall nucleate boing -> larger heat transfer

	if(lambda==0 && mu ==0 && _heatTransfertCoeff==0)
		return;
	if(isBord)
	{
		for(int k=0; k<_nVar; k++)
			_Vj[k] = _Vextdiff[k];
		_inv_dxj=_inv_dxi;
	}

	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_nVar-1; k++)
		_phi[k] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*mu*(_porosityj*_Vj[k] - _porosityi*_Vi[k]);
	//contribution visqueuse sur l'energie
	_phi[_nVar-1] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*lambda*(_porosityj*_Vj[_nVar-1] - _porosityi*_Vi[_nVar-1]);
	
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

void NavierStokes::addSourceTermToSecondMember(const int i, int nbVoisinsi,const int j, int nbVoisinsj, bool isBord, int ij, double mesureFace)//To do : generalise to unstructured meshes
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"NavierStokes::addSourceTerm cell i= "<<i<< " cell j= "<< j<< " isbord "<<isBord<<endl;

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

void NavierStokes::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{
	double P = Vi[0];
	double h = Vi[ _nVar -1];
	double phirho=Ui[0], T =_fluides[0]->getTemperatureFromEnthalpy(P, h);
	double norm_u=0;
	for(int k=0; k<_Ndim; k++)
		norm_u+=Vi[1+k]*Vi[1+k];
	norm_u=sqrt(norm_u);
	if(T>_Tsat)
		Si[0]=_heatPowerField(i)/_latentHeat;
	else
		Si[0]=0;
	for(int k=1; k<_nVar-1; k++)
		Si[k]  =(_gravite[k]-_dragCoeffs[0]*norm_u*Vi[1+k])*phirho;

	Si[_nVar-1]=_heatPowerField(i);

	for(int k=0; k<_Ndim; k++)
		Si[_nVar-1] +=(_GravityField3d[k]-_dragCoeffs[0]*norm_u*Vi[1+k])*Vi[1+k]*phirho;

	if(_timeScheme==Implicit)
	{
		for(int k=0; k<_nVar*_nVar;k++)
			_GravityImplicitationMatrix[k] = 0;
		if(!_usePrimitiveVarsInNewton)
			for(int k=0; k<_nVar;k++)
				_GravityImplicitationMatrix[k*_nVar]=-_gravite[k];
		else
		{
			double P = Vi[ 0 ];
			double h = Vi[_nVar -1];
			double Temperature = _fluides[0]->getTemperatureFromEnthalpy(p,h);
			double drho_sur_dp = _fluides[0]->getDrhoDP_e(pression, Temperature);
			double drho_sur_dT = _fluides[0]->getDrhoDT_P(pression, Temperature)
			for(int k=0; k<_nVar;k++)
			{
				_GravityImplicitationMatrix[k*_nVar+0]      =-_gravite[k] * drho_sur_dp;
				_GravityImplicitationMatrix[k*_nVar+_nVar-1]=-_gravite[k] * drho_sur_dT;
			}
		}
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::sourceVector"<<endl;
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
		if(_timeScheme==Implicit)
			displayMatrix(_GravityImplicitationMatrix, _nVar, "Gravity implicitation matrix");
	}
}

void NavierStokes::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
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
		cout<<"SinglePhase::pressureLossVector K= "<<K<<endl;
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

void NavierStokes::porosityGradientSourceVector()
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
	_porosityGradientSourceVector[_nVar-1]=0;
}

void NavierStokes::jacobian(const int &j, string nameOfGroup,double * normale)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite convection bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_Jcb[k] = 0;//No implicitation at this stage

	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);
	double q_n=0;//quantité de mouvement normale à la paroi
	for(k=0; k<_Ndim; k++)
		q_n+=_externalStates[(k+1)]*normale[k];

	// loop of boundary types
	if (_limitField[nameOfGroup].bcType==Wall)
	{
		for(k=0; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
		for(k=1; k<_nVar-1;k++)
			for(int l=1; l<_nVar-1;l++)
				_Jcb[k*_nVar + l] -= 2*normale[k-1]*normale[l-1];
	}
	else if (_limitField[nameOfGroup].bcType==Inlet)
	{
		return;
		if(q_n<0){
			double v[_Ndim], ve[_Ndim], v2, ve2;

			_idm[0] = _nVar*j;
			for(k=1; k<_nVar; k++)
				_idm[k] = _idm[k-1] + 1;
			VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
			VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

			ve[0] = _limitField[nameOfGroup].v_x[0];
			v[0]=_Vj[1];
			ve2 = ve[0]*ve[0];
			v2 = v[0]*v[0];
			if (_Ndim >1){
				ve[1] = _limitField[nameOfGroup].v_y[0];
				v[1]=_Vj[2];
				ve2 += ve[1]*ve[1];
				v2 = v[1]*v[1];
			}
			if (_Ndim >2){
				ve[2] = _limitField[nameOfGroup].v_z[0];
				v[2]=_Vj[3];
				ve2 += ve[2]*ve[2];
				v2 = v[2]*v[2];
			}
			double internal_energy=_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,_Uj[0]);
			double total_energy=internal_energy+ve2/2;

			//Mass line
			_Jcb[0]=v2/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[1+k]=-v[k]/internal_energy;
			_Jcb[_nVar-1]=1/internal_energy;
			//Momentum lines
			for(int l =1;l<1+_Ndim;l++){
				_Jcb[l*_nVar]=v2*ve[l-1]/(2*internal_energy);
				for(k=0; k<_Ndim;k++)
					_Jcb[l*_nVar+1+k]=-v[k]*ve[l-1]/internal_energy;
				_Jcb[l*_nVar+_nVar-1]=ve[l-1]/internal_energy;
			}
			//Energy line
			_Jcb[(_nVar-1)*_nVar]=v2*total_energy/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[(_nVar-1)*_nVar+1+k]=-v[k]*total_energy/internal_energy;
			_Jcb[(_nVar-1)*_nVar+_nVar-1]=total_energy/internal_energy;
		}
		else
			for(k=0;k<_nVar;k++)
				_Jcb[k*_nVar+k]=1;
		//Old jacobian
		/*
		 _Jcb[0] = 1;
		_Jcb[_nVar]=_limitField[nameOfGroup].v_x[0];//Kieu
		v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
		if(_Ndim>1)
		{
			_Jcb[2*_nVar]= _limitField[nameOfGroup].v_y[0];
			v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			if(_Ndim==3){
				_Jcb[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
			}
		}
		_Jcb[(_nVar-1)*_nVar]=_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2;
		 */
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure && q_n<0){
		return;
		double v[_Ndim], v2=0;
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

		for(k=0; k<_Ndim;k++){
			v[k]=_Vj[1+k];
			v2+=v[k]*v[k];
		}

		double rho_ext=_compressibleFluid->getDensity(_limitField[nameOfGroup].p, _limitField[nameOfGroup].T);
		double rho_int = _externalStates[0];
		double density_ratio=rho_ext/rho_int;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_Jcb[l*_nVar]=-density_ratio*v[l-1];
			_Jcb[l*_nVar+l]=density_ratio;
		}
		//Energy lines
		_Jcb[(_nVar-1)*_nVar]=-v2*density_ratio;
		for(k=0; k<_Ndim;k++)
			_Jcb[(_nVar-1)*_nVar+1+k]=density_ratio*v[k];
	}
	// not wall, not inlet, not inletPressure
	else if(_limitField[nameOfGroup].bcType==Outlet || (_limitField[nameOfGroup].bcType==InletPressure && q_n>=0))
	{
		return;
		double v[_Ndim], v2=0;
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

		for(k=0; k<_Ndim;k++){
			v[k]=_Vj[1+k];
			v2+=v[k]*v[k];
		}

		double rho_ext=_compressibleFluid->getDensity(_limitField[nameOfGroup].p, _externalStates[_nVar-1]);
		double rho_int = _externalStates[0];
		double density_ratio=rho_ext/rho_int;
		double internal_energy=_compressibleFluid->getInternalEnergy(_externalStates[_nVar-1],rho_int);
		double total_energy=internal_energy+v2/2;

		//Mass line
		_Jcb[0]=density_ratio*(1-v2/(2*internal_energy));
		for(k=0; k<_Ndim;k++)
			_Jcb[1+k]=density_ratio*v[k]/internal_energy;
		_Jcb[_nVar-1]=-density_ratio/internal_energy;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_Jcb[l*_nVar]=density_ratio*v2*v[l-1]/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_Jcb[l*_nVar+1+k]=density_ratio*v[k]*v[l-1]/internal_energy;
			_Jcb[l*_nVar+1+k]-=density_ratio;
			_Jcb[l*_nVar+_nVar-1]=-density_ratio*v[l-1]/internal_energy;
		}
		//Energy line
		_Jcb[(_nVar-1)*_nVar]=density_ratio*v2*total_energy/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_Jcb[(_nVar-1)*_nVar+1+k]=density_ratio*v[k]*total_energy/internal_energy;
		_Jcb[(_nVar-1)*_nVar+_nVar-1]=density_ratio*(1-total_energy/internal_energy);
		//Old jacobian
		/*
		int idim,jdim;
		double cd = 1,cn=0,p0, gamma;
		_idm[0] = j*_nVar;// Kieu
		for(k=1; k<_nVar;k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_conservativeVars, _nVar, _idm, _phi);
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		// compute the common numerator and common denominator
		p0=_compressibleFluid->constante("p0");
		gamma =_compressibleFluid->constante("gamma");
		cn =_limitField[nameOfGroup].p +p0;
		cd = _phi[0]*_compressibleFluid->getInternalEnergy(_externalStates[_nVar-1],rho)-p0;
		cd*=cd;
		cd*=(gamma-1);
		//compute the v2
		for(k=1; k<_nVar-1;k++)
			v2+=_externalStates[k]*_externalStates[k];
		// drho_ext/dU
		_JcbDiff[0] = cn*(_phi[_nVar-1] -v2 -p0)/cd;
		for(k=1; k<_nVar-1;k++)
			_JcbDiff[k]=cn*_phi[k]/cd;
		_JcbDiff[_nVar-1]= -cn*_phi[0]/cd;
		//dq_ext/dU
		for(idim=0; idim<_Ndim;idim++)
		{
			//premiere colonne
			_JcbDiff[(1+idim)*_nVar]=-(v2*cn*_phi[idim+1])/(2*cd);
			//colonnes intermediaire
			for(jdim=0; jdim<_Ndim;jdim++)
			{
				_JcbDiff[(1+idim)*_nVar + jdim + 1] =_externalStates[idim+1]*_phi[jdim+1];
				_JcbDiff[(1+idim)*_nVar + jdim + 1]*=cn/cd;
			}
			//matrice identite*cn*(rhoe- p0)
			_JcbDiff[(1+idim)*_nVar + idim + 1] +=( cn*(_phi[0]*_compressibleFluid->getInternalEnergy(_externalStates[_nVar-1],rho)-p0))/cd;

			//derniere colonne
			_JcbDiff[(1+idim)*_nVar + _nVar-1]=-_phi[idim+1]*cn/cd;
		}
		//drhoE/dU
		_JcbDiff[_nVar*(_nVar-1)] = -(v2*_phi[_nVar -1]*cn)/(2*cd);
		for(int idim=0; idim<_Ndim;idim++)
			_JcbDiff[_nVar*(_nVar-1)+idim+1]=_externalStates[idim+1]*_phi[_nVar -1]*cn/cd;
		_JcbDiff[_nVar*_nVar -1] = -(v2/2+p0)*cn/cd;
		 */
	}
	else  if (_limitField[nameOfGroup].bcType!=Neumann)// not wall, not inlet, not outlet
	{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::jacobian: This boundary condition is not treated");
	}
}

//calcule la jacobienne pour une CL de diffusion
void  NavierStokes::jacobianDiff(const int &j, string nameOfGroup)
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"Jacobienne condition limite diffusion bord "<< nameOfGroup<<endl;

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;

	if (_limitField[nameOfGroup].bcType==Wall){
		double v[_Ndim], ve[_Ndim], v2, ve2;

		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

		ve[0] = _limitField[nameOfGroup].v_x[0];
		v[0]=_Vj[1];
		ve2 = ve[0]*ve[0];
		v2 = v[0]*v[0];
		if (_Ndim >1){
			ve[1] = _limitField[nameOfGroup].v_y[0];
			v[1]=_Vj[2];
			ve2 += ve[1]*ve[1];
			v2 = v[1]*v[1];
		}
		if (_Ndim >2){
			ve[2] = _limitField[nameOfGroup].v_z[0];
			v[2]=_Vj[3];
			ve2 += ve[2]*ve[2];
			v2 = v[2]*v[2];
		}
		double rho=_Uj[0];
		double internal_energy=_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,rho);
		double total_energy=internal_energy+ve2/2;

		//Mass line
		_JcbDiff[0]=v2/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_JcbDiff[1+k]=-v[k]/internal_energy;
		_JcbDiff[_nVar-1]=1/internal_energy;
		//Momentum lines
		for(int l =1;l<1+_Ndim;l++){
			_JcbDiff[l*_nVar]=v2*ve[l-1]/(2*internal_energy);
			for(k=0; k<_Ndim;k++)
				_JcbDiff[l*_nVar+1+k]=-v[k]*ve[l-1]/internal_energy;
			_JcbDiff[l*_nVar+_nVar-1]=ve[l-1]/internal_energy;
		}
		//Energy line
		_JcbDiff[(_nVar-1)*_nVar]=v2*total_energy/(2*internal_energy);
		for(k=0; k<_Ndim;k++)
			_JcbDiff[(_nVar-1)*_nVar+1+k]=-v[k]*total_energy/internal_energy;
		_JcbDiff[(_nVar-1)*_nVar+_nVar-1]=total_energy/internal_energy;
	}
	else if (_limitField[nameOfGroup].bcType==Outlet || _limitField[nameOfGroup].bcType==Neumann
			||_limitField[nameOfGroup].bcType==Inlet || _limitField[nameOfGroup].bcType==InletPressure)
	{
		for(k=0;k<_nVar;k++)
			_JcbDiff[k*_nVar+k]=1;
	}
	else{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		*_runLogFile<<"group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::jacobianDiff: This boundary condition is not recognised");
	}
}

void NavierStokes::primToCons(const double *P, const int &i, double *W, const int &j){
	//cout<<"SinglePhase::primToCons i="<<i<<", j="<<j<<", P[i*(_Ndim+2)]="<<P[i*(_Ndim+2)]<<", P[i*(_Ndim+2)+(_Ndim+1)]="<<P[i*(_Ndim+2)+(_Ndim+1)]<<endl;
	double rho = _fluides[0]->getDensityFromEnthalpy( P[i*(_Ndim+2)], P[i*(_Ndim+2)+(_Ndim+1)] );
	double T = _fluides[0]->getTemperatureFromEnthalpy( _Vdiff[ _nVar -1], rho );

	assert( P[i*(_Ndim+2)]>0);//Pressure should be positive
	assert( T > 0);//Temperature should be positive

	W[j*(_Ndim+2)] =  _porosityField(j)*rho;//phi*rho
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(k+1)] = W[j*(_Ndim+2)] * P[i*(_Ndim+2)+(k+1)];//phi*rho*u

	W[j*(_Ndim+2)+(_Ndim+1)] = W[j*(_Ndim+2)] * _fluides[0]->getInternalEnergy(P[i*(_Ndim+2)+ (_Ndim+1)],rho);//rho*e
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(_Ndim+1)] += W[j*(_Ndim+2)]*P[i*(_Ndim+2)+(k+1)]*P[i*(_Ndim+2)+(k+1)]*0.5;//phi*rho*e+0.5*phi*rho*u^2
}

void NavierStokes::primToConsJacobianMatrix(double *V)
{
	double P = V[0];
	double h = V[_nVar-1];
	double rho = _fluides[0]->getDensityFromEnthalpy(P, h);
	double T = _fluides[0]->getTemperatureFromEnthalpy(P,h);
	double e = h - pression/rho.0;
	double vitesse[_Ndim];

	for(int idim=0;idim<_Ndim;idim++)
		vitesse[idim]=V[1+idim];
	double v2=0;
	for(int idim=0;idim<_Ndim;idim++)
		v2+=vitesse[idim]*vitesse[idim];
	
	double E = e + v2/2.0 ;

	double drho_sur_dh = _fluides[0]->getDrhoDh_p(P, T);
	double drho_sur_dp = _fluides[0]->getDrhoDP_h(P, T); 
	double de_sur_dp = _fluides[0]->getDeDp_h(P,T);
	double de_sur_dh = _fluides[0]->getDeDh_p(P,T);
	double drhoE_sur_dp = drho_sur_dp * E + rho * (de_sur_dp + v2); 
	double drhoE_sur_dh = drho_sur_dh * E + rho * (de_sur_dh + v2)

	getDensityDerivatives(pression, h);
	for(int k=0;k<_nVar*_nVar; k++){
		_primToConsJacoMat[k]=0;
	}
	//première ligne
	_primToConsJacoMat[0] = drho_sur_dp;
	_primToConsJacoMat[_nVar -1] = drho_sur_dh;
	//Deuxième ligne à (Ndim+1) ligne
	for(int i=0; i< _Ndim; i++){
		_primToConsJacoMat[i*_nVar] = drho_sur_dh * vitesse[i];
		_primToConsJacoMat[i*_nVar + _nVar-1] = drho_sur_dh * vitesse[i];
		_primToConsJacoMat[i*_nVar + 1+i] = rho;
	}
	//Dernière ligne

	_primToConsJacoMat[(_nVar-1) * _nVar ] = drhoE_sur_dp;
	_primToConsJacoMat[(_nVar-1) * _nVar + _nVar -1] = drhoE_sur_dh;
	for(int j=0; j< _Ndim; j++){
		_primToConsJacoMat[(_nVar-1) * _nVar + 1+j] = rho * vitesse[i];

	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" NavierStokes::primToConsJacobianMatrix" << endl;
		displayVector(_Vi,_nVar," _Vi " );
		cout<<" Jacobienne primToCons: " << endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobienne primToCons: ");
	}
}

Vector NavierStokes::staggeredVFFCFlux()
{
	*_runLogFile<< "NavierStokes::staggeredVFFCFlux is not yet available "<<  endl;
	_runLogFile->close();
	throw CdmathException("NavierStokes::staggeredVFFCFlux is not yet available");
}
void NavierStokes::entropicShift(double* n)
{
		*_runLogFile<< "NavierStokes::entropicShift is not yet available "<<  endl;
	_runLogFile->close();
	throw CdmathException("NavierStokes::entropicShift is not yet available");
}

Vector NavierStokes::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"NavierStokes::convectionFlux start"<<endl;
		cout<<"Ucons"<<endl;
		cout<<U<<endl;
		cout<<"Vprim"<<endl;
		cout<<V<<endl;
	}
	double h = V[_nVar-1];
	double phirho=U[0];//phi rho
	Vector phiq[_Ndim];//phi rho u
	for(int i=0;i<_Ndim;i++)
		phiq[i]=U[1+i];

	double pression=V[0];
	Vector vitesse[_Ndim];
	for(int i=0;i<_Ndim;i++)
		vitesse[i]=V[1+i];
		v2 += vitesse[i] * vitesse[i];

	double vitessen=vitesse*normale;
	double rho=phirho/porosity;

	Vector F[_nVar];
	F[0] = phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F[1+i] = phirho*vitessen*vitesse[i]+pression*normale[i]*porosity;
	F[1+_Ndim] = phirho*( h - pression/rho.0 + v2/2.0)*vitessen;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void NavierStokes::ConvectionMatrixConservativeVariables(double u_n, double H,Vector velocity, double k, double K)
{
	/******** Construction de la matrice de Roe *********/
	//premiere ligne (masse)
	_Aroe[0]=0;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[1+idim]=_vec_normal[idim];
	_Aroe[_nVar-1]=0;

	//lignes intermadiaires(qdm)
	for(int idim=0; idim<_Ndim;idim++)
	{
		//premiere colonne
		_Aroe[(1+idim)*_nVar]=K*_vec_normal[idim] - u_n*_Uroe[1+idim];
		//colonnes intermediaires
		for(int jdim=0; jdim<_Ndim;jdim++)
			_Aroe[(1+idim)*_nVar + jdim + 1] = _Uroe[1+idim]*_vec_normal[jdim]-k*_vec_normal[idim]*_Uroe[1+jdim];
		//matrice identite
		_Aroe[(1+idim)*_nVar + idim + 1] += u_n;
		//derniere colonne
		_Aroe[(1+idim)*_nVar + _nVar-1]=k*_vec_normal[idim];
	}
	//derniere ligne (energie)
	_Aroe[_nVar*(_nVar-1)] = (K - H)*u_n;
	for(int idim=0; idim<_Ndim;idim++)
		_Aroe[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] - k*u_n*_Uroe[idim+1];
	_Aroe[_nVar*_nVar -1] = (1 + k)*u_n;
}

void NavierStokes::convectionMatrixPrimitiveVariables( double rho, double u_n, double H, Vector vitesse)
{
	//G(V_L)-G(V_R)=Aroe_prim (V_L-V_R)

	//TODO bien défini ?
	double P = _V[0];
	double h = _V[_nVar-1];
	double T = _fluides[0]->getTemperatureFromEnthalpy(P,h);
	double vitesse[_Ndim];

	for(int idim=0;idim<_Ndim;idim++)
		vitesse[idim]=V[1+idim];
	double v2=0;
	for(int idim=0;idim<_Ndim;idim++)
		v2+=vitesse[idim]*vitesse[idim];
	
	//TODO quelles valeurs de P et T ?
	double drho_sur_dh = _fluides[0]->getDrhoDh_p(P, T);
	double drho_sur_dp = _fluides[0]->getDrhoDP_h(P, T); 

	//Première ligne 
	_AroeImplicit[0*_nVar+0] = u_n * drho_sur_dp;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i] = rho * _vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim] = u_n * drho_sur_dh;
	//Deuxième ligne à (1 + _Ndim) ligne
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0] = drho_sur_dp * u_n * vitesse[i] + _vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j] = rho * vitesse[i] *_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i] += rho * u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim] = drho_sur_dh * u_n * vitesse[i];
	}
	//Dernière ligne
	_AroeImplicit[(1+_Ndim)*_nVar+0] = drho_sur_dp * u_n * H;
	for(int i=0;i<_Ndim;i++)1
		_AroeImplicit[(1+_Ndim)*_nVar+1+i] = rho * (H * _vec_normal[i] + u_n * vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim] = u_n * (H * drho_sur_dh + rho);
}

void NavierStokes::staggeredRoeUpwindingMatrixPrimitiveVariables(double rho, double u_n, double H, Vector vitesse)
{
	//Calcul de décentrement de type décalé en variables primitives
	//TODO à définir
	getDensityDerivatives(p,h);
	_AroeImplicit[0*_nVar+0] =_drho_sur_dp * u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i] =rho * _vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim] =_drho_sur_dT * u_n;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0] =_drho_sur_dp * u_n*vitesse[i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j] = rho * vitesse[i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i] += rho * u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim] =_drho_sur_dT * u_n * vitesse[i];
	}
	_AroeImplicit[(1+_Ndim)*_nVar+0]= (_drhoE_sur_dp+1) * u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(1+_Ndim)*_nVar+1+i] = rho * ( H *_vec_normal[i]+u_n * vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim] = _drhoE_sur_dT * u_n;
}


void NavierStokes::testConservation()
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

void NavierStokes::save(){
    PetscPrintf(PETSC_COMM_WORLD,"Saving numerical results at time step number %d \n\n", _nbTimeStep);
    *_runLogFile<< "Saving numerical results at time step number "<< _nbTimeStep << endl<<endl;

	string prim(_path+"/SinglePhasePrim_");///Results
	string cons(_path+"/SinglePhaseCons_");
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
			Ii = i*_nVar +_nVar-1;
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

			h   = _compressibleFluid->getEnthalpy(T,rho);

			_Enthalpy(i)=h;
			_Density(i)=rho;
			_Pressure(i)=p;
			_Temperature(i)=T;
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
			_MachNumber(i)=sqrt(v2)/_compressibleFluid->vitesseSonEnthalpie(h);
		}
		_Enthalpy.setTime(_time,_nbTimeStep);
		_Density.setTime(_time,_nbTimeStep);
		_Pressure.setTime(_time,_nbTimeStep);
		_Temperature.setTime(_time,_nbTimeStep);
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
				_Enthalpy.writeVTK(allFields+"_Enthalpy");
				_Density.writeVTK(allFields+"_Density");
				_Pressure.writeVTK(allFields+"_Pressure");
				_Temperature.writeVTK(allFields+"_Temperature");
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
				_Enthalpy.writeMED(allFields+"_Enthalpy");
				_Density.writeMED(allFields+"_Density");
				_Pressure.writeMED(allFields+"_Pressure");
				_Temperature.writeMED(allFields+"_Temperature");
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
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
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
				_Enthalpy.writeVTK(allFields+"_Enthalpy",false);
				_Density.writeVTK(allFields+"_Density",false);
				_Pressure.writeVTK(allFields+"_Pressure",false);
				_Temperature.writeVTK(allFields+"_Temperature",false);
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
				_Enthalpy.writeMED(allFields+"_Enthalpy",false);
				_Density.writeMED(allFields+"_Density",false);
				_Pressure.writeMED(allFields+"_Pressure",false);
				_Temperature.writeMED(allFields+"_Temperature",false);
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
				_Enthalpy.writeCSV(allFields+"_Enthalpy");
				_Density.writeCSV(allFields+"_Density");
				_Pressure.writeCSV(allFields+"_Pressure");
				_Temperature.writeCSV(allFields+"_Temperature");
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

Field& NavierStokes::getPressureField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getPressureField, Call initialize first");

	if(!_saveAllFields)
	{
		_Pressure=Field("Pressure",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&_Pressure(i));
		}
		_Pressure.setTime(_time,_nbTimeStep);
	}
	return _Pressure;
}

Field& NavierStokes::getTemperatureField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getTemperatureField, Call initialize first");

	if(!_saveAllFields)
	{
		_Temperature=Field("Temperature",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&_Temperature(i));
		}
		_Temperature.setTime(_time,_nbTimeStep);
	}
	return _Temperature;
}

Field& NavierStokes::getVelocityField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getVelocityField, Call initialize first");

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

Field& NavierStokes::getMachNumberField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getMachNumberField, Call initialize first");

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
	
			rho=_compressibleFluid->getDensity(p,T);
			h  =_compressibleFluid->getEnthalpy(T,rho);
			_MachNumber[i]  =sqrt(u2)/_compressibleFluid->vitesseSonEnthalpie(h);
			//cout<<"u="<<sqrt(u2)<<", c= "<<_compressibleFluid->vitesseSonEnthalpie(h)<<", MachNumberField[i] = "<<MachNumberField[i] <<endl;
		}
		_MachNumber.setTime(_time,_nbTimeStep);
	}
	//cout<<", MachNumberField = "<<MachNumberField <<endl;

	return _MachNumber;
}

Field& NavierStokes::getVelocityXField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getVelocityXField, Call initialize first");

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

Field& NavierStokes::getVelocityYField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getVelocityYField, Call initialize first");

	if(_Ndim<2)
        throw CdmathException("SinglePhase::getVelocityYField() error : dimension should be at least 2");	
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

Field& NavierStokes::getVelocityZField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getVelocityZField, Call initialize first");

	if(_Ndim<3)
        throw CdmathException("SinglePhase::getvelocityZField() error : dimension should be 3");	
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

Field& NavierStokes::getDensityField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getDensityField, Call initialize first");
		
	if(!_saveAllFields )
	{
		_Density=Field("Density",CELLS,_mesh,1);
		int Ii;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_conservativeVars,1,&Ii,&_Density(i));
		}
		_Density.setTime(_time,_nbTimeStep);
	}
	return _Density;
}

Field& NavierStokes::getMomentumField()//not yet managed by parameter _saveAllFields
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getMomentumField, Call initialize first");

	_Momentum=Field("Momentum",CELLS,_mesh,_Ndim);
	int Ii;
	for (long i = 0; i < _Nmailles; i++)
		for (int j = 0; j < _Ndim; j++)//On récupère les composantes de qdm
		{
			int Ii = i*_nVar +1+j;
			VecGetValues(_conservativeVars,1,&Ii,&_Momentum(i,j));
		}
	_Momentum.setTime(_time,_nbTimeStep);

	return _Momentum;
}

Field& NavierStokes::getTotalEnergyField()//not yet managed by parameter _saveAllFields
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getTotalEnergyField, Call initialize first");

	_TotalEnergy=Field("TotalEnergy",CELLS,_mesh,1);
	int Ii;
	for (long i = 0; i < _Nmailles; i++){
		Ii = i*_nVar +_nVar-1;
		VecGetValues(_conservativeVars,1,&Ii,&_TotalEnergy(i));
	}
	_TotalEnergy.setTime(_time,_nbTimeStep);

	return _TotalEnergy;
}

Field& NavierStokes::getEnthalpyField()
{
	if(!_initializedMemory)
		throw CdmathException("SinglePhase::getEnthalpyField, Call initialize first");

	if(!_saveAllFields )
	{
		_Enthalpy=Field("Enthalpy",CELLS,_mesh,1);
		int Ii;
		double p,T,rho;
		for (long i = 0; i < _Nmailles; i++){
			Ii = i*_nVar;
			VecGetValues(_primitiveVars,1,&Ii,&p);
			Ii = i*_nVar +_nVar-1;
			VecGetValues(_primitiveVars,1,&Ii,&T);
			
			rho=_compressibleFluid->getDensity(p,T);
			_Enthalpy(i)=_compressibleFluid->getEnthalpy(T,rho);
		}
		_Enthalpy.setTime(_time,_nbTimeStep);
	}

	return _Enthalpy;
}

vector<string> NavierStokes::getOutputFieldsNames()
{
	vector<string> result(8);
	
	result[0]="Pressure";
	result[1]="Velocity";
	result[2]="Temperature";
	result[3]="Density";
	result[4]="Momentum";
	result[5]="TotalEnergy";
	result[6]="Enthalpy";
	result[7]="VelocityX";
	
	return result;
}

Field& NavierStokes::getOutputField(const string& nameField )
{
	if(nameField=="pressure" || nameField=="Pressure" || nameField=="PRESSURE" || nameField=="PRESSION" || nameField=="Pression"  || nameField=="pression" )
		return getPressureField();
	else if(nameField=="velocity" || nameField=="Velocity" || nameField=="VELOCITY" || nameField=="Vitesse" || nameField=="VITESSE" || nameField=="vitesse" )
		return getVelocityField();
	else if(nameField=="velocityX" || nameField=="VelocityX" || nameField=="VELOCITYX" || nameField=="VitesseX" || nameField=="VITESSEX" || nameField=="vitesseX" )
		return getVelocityXField();
	else if(nameField=="velocityY" || nameField=="VelocityY" || nameField=="VELOCITYY" || nameField=="VitesseY" || nameField=="VITESSEY" || nameField=="vitesseY" )
		return getVelocityYField();
	else if(nameField=="velocityZ" || nameField=="VelocityZ" || nameField=="VELOCITYZ" || nameField=="VitesseZ" || nameField=="VITESSEZ" || nameField=="vitesseZ" )
		return getVelocityZField();
	else if(nameField=="temperature" || nameField=="Temperature" || nameField=="TEMPERATURE" || nameField=="temperature" )
		return getTemperatureField();
	else if(nameField=="density" || nameField=="Density" || nameField=="DENSITY" || nameField=="Densite" || nameField=="DENSITE" || nameField=="densite" )
		return getDensityField();
	else if(nameField=="momentum" || nameField=="Momentum" || nameField=="MOMENTUM" || nameField=="Qdm" || nameField=="QDM" || nameField=="qdm" )
		return getMomentumField();
	else if(nameField=="enthalpy" || nameField=="Enthalpy" || nameField=="ENTHALPY" || nameField=="Enthalpie" || nameField=="ENTHALPIE" || nameField=="enthalpie" )
		return getEnthalpyField();
	else if(nameField=="totalenergy" || nameField=="TotalEnergy" || nameField=="TOTALENERGY" || nameField=="ENERGIETOTALE" || nameField=="EnergieTotale" || nameField=="energietotale" )
		return getTotalEnergyField();
    else
    {
        cout<<"Error : Field name "<< nameField << " does not exist, call getOutputFieldsNames first to check" << endl;
        throw CdmathException("SinglePhase::getOutputField error : Unknown Field name");
    }
}
