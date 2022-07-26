/*
 * IsothermalSinglePhase.cxx
 *
 *  Created on: July 7 2022
 */

#include "IsothermalSinglePhase.hxx"
#include "StiffenedGas.hxx"

using namespace std;

IsothermalSinglePhase::IsothermalSinglePhase(phaseType fluid, pressureEstimate pEstimate, int dim){
	_Ndim=dim;
	_nVar=_Ndim+1;
	_nbPhases = 1;
	_dragCoeffs=vector<double>(1,0);
	_fluides.resize(1);
	if (pEstimate==around1bar300K)//EOS at 1 bar and 300K
	{
		_Temperature=300;//Constant temperature of the model
		if(fluid==Gas){
			cout<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is air around 1 bar and 300 K (27°C)"<<endl;
			_internalEnergy=2.22e5;//nitrogen internal energy at 1bar, 300K
			_fluides[0] = new StiffenedGas(1.4,743,_Temperature,_internalEnergy);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, c_v=743
		}
		else{
			cout<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			*_runLogFile<<"Fluid is water around 1 bar and 300 K (27°C)"<<endl;
			_internalEnergy=1.12e5;//water internal energy at 1 bar, 300K
			_fluides[0] = new StiffenedGas(996,1e5,_Temperature,_internalEnergy,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C
		}
	}
	else
	{
		if(fluid==Gas){//EOS at 155 bars and 618K
			cout<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			*_runLogFile<<"Fluid is Gas around saturation point 155 bars and 618 K (345°C)"<<endl;
			_Temperature=618;//Constant temperature of the model
			_internalEnergy=2.44e6;//Gas internal energy at saturation at 155 bar
			_fluides[0] = new StiffenedGas(102,1.55e7,_Temperature,_internalEnergy, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C:
		}
		else{//EOS at 155 bars and 573K
			cout<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			*_runLogFile<<"Fluid is water around saturation point 155 bars and 573 K (300°C)"<<endl;
			_Temperature=573;//Constant temperature of the model
			_internalEnergy=1.3e6;//water internal energy at saturation at 155 bar
			_fluides[0] = new StiffenedGas(726.82,1.55e7,_Temperature,_internalEnergy, 971.,5454.);  //stiffened gas law for water at pressure 155 bar and temperature 300°C:
		}
	}

	_fileName = "SolverlabIsothermalSinglePhase";
    PetscPrintf(PETSC_COMM_WORLD,"\n Isothermal single phase problem \n");
    
    _usePrimitiveVarsInNewton=true;//This class is designed only to solve linear system in primitive variables
     _Vdiff=NULL;
}

void IsothermalSinglePhase::initialize(){
	cout<<"\n Initialising the isothermal single phase model\n"<<endl;
	*_runLogFile<<"\n Initialising the isothermal single phase model\n"<<endl;

	_Uroe = new double[_nVar+1];//Deleted in ProblemFluid::terminate()
	_Vextdiff= new double[_nVar];
	_Vext= new double[_nVar];

	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
		_gravite[i+1]=_GravityField3d[i];

	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];//Deleted in ProblemFluid::terminate()

	_Vdiff = new double[_nVar];
	
	if(_saveVelocity)
		_Vitesse=Field("Velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant

	ProblemFluid::initialize();
}

void IsothermalSinglePhase::terminate(){
	delete[] _Vdiff,_Vextdiff,_Vext;
	ProblemFluid::terminate();
}

bool IsothermalSinglePhase::iterateTimeStep(bool &converged)
{    //The class does not allow the use of conservative variables in Newton iterations
	if(_timeScheme == Explicit || _usePrimitiveVarsInNewton)
		return ProblemFluid::iterateTimeStep(converged);
	else
		throw CdmathException("IsothermalSinglePhase can not use conservative variables in Newton scheme for implicit in time discretisation");
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
		_Vj=_Vext;

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Convection state: variables primitives maille " <<j <<endl;
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q] << endl;
		cout << endl;
	}
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
		cout << "Densité moyenne Roe  gauche " << i << ": " << ri*ri << ", droite " << j << ": " << rj*rj << "->" << _Uroe[0] << endl;
	for(int k=0;k<_Ndim;k++){
		xi = _Ui[k+1];
		xj = _Uj[k+1];
		_Uroe[1+k] = (xi/ri + xj/rj)/(ri + rj);
		//"moyenne" des vitesses
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout << "Vitesse de Roe composante "<< k<<"  gauche " << i << ": " << xi/(ri*ri) << ", droite " << j << ": " << xj/(rj*rj) << "->" << _Uroe[k+1] << endl;
	}

	//Computation of 1/c²// Todo :  add porosity in the sound speed formula
	_Uroe[_nVar] = (_Ui[0]-_Uj[0])/(_Vi[0]-_Vj[0]);//_Uroe has size _nVar+1 !!!
	
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar+1;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
	
}

void IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//entree: indices des cellules Vi et Vj.
	//Attention : Vj peut être un état fantôme
	//sortie: matrices et etat de diffusion (rho, q) ou (p, v) ?
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

	if(IsBord)
		_Vj=_Vextdiff;
	else
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices cellule gauche" << i << endl;
		cout << "Vi = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q]  << "\t";
		cout << endl;
		cout << "IsothermalSinglePhase::diffusionPrimitiveStateAndMatrices cellule droite" << j << endl;
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
		for(int i=1;i<(_nVar-1);i++)
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

	//Todo : treat correctly the computation of eigenvalues for incompressible flows
	double c;
	if(_Uroe[_nVar]==0.)//infinite sound speed
		if(_timeScheme==Explicit)
			throw CdmathException("Explicit scheme cannot be used for incompressible fluids");
		else
			c=0.;
	else
		c=1./sqrt(_Uroe[_nVar]);
		
	vector<std::complex<double>>vp_dist(3);
	vp_dist[0]=u_n-c;vp_dist[1]=u_n;vp_dist[2]=u_n+c;
	
	_maxvploc=fabs(u_n)+c;
	if(_maxvploc>_maxvp)
		_maxvp=_maxvploc;

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"IsothermalSinglePhase::convectionMatrices Eigenvalues "<<u_n-c<<" , "<<u_n<<" , "<<u_n+c<<endl;

	convectionMatrixPrimitiveVariables(u_n);

	if(_entropicCorrection)
	{
		*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: entropy scheme not available for IsothermalSinglePhase"<<endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::convectionMatrices: entropy scheme not available for IsothermalSinglePhase");
	}

	/******** Construction des matrices de decentrement ********/
	if( _spaceScheme ==centered){
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] = 0;
	}
	else if(_spaceScheme == upwind )
	{
		vector< complex< double > > y (3,0);
		for( int i=0 ; i<3 ; i++)
			y[i] = Polynoms::abs_generalise(vp_dist[i]);
		Polynoms::abs_par_interp_directe(3,vp_dist, _Aroe, _nVar,_precision, _absAroe,y);

	}
	else if( _spaceScheme ==staggered ){
		if(_entropicCorrection)//To do: study entropic correction for staggered
		{
			*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
			_runLogFile->close();
			throw CdmathException("IsothermalSinglePhase::convectionMatrices: entropy scheme not available for staggered scheme");
		}

		staggeredRoeUpwindingMatrixPrimitiveVariables( u_n);
	}
	else
	{
		*_runLogFile<<"IsothermalSinglePhase::convectionMatrices: scheme not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::convectionMatrices: scheme not treated");
	}

	for(int i=0; i<_nVar*_nVar;i++)
	{
		_AroeMinus[i] = (_Aroe[i]-_absAroe[i])/2;
		_AroePlus[i]  = (_Aroe[i]+_absAroe[i])/2;
	}
	if(_timeScheme==Implicit)
	{
		for(int i=0; i<_nVar*_nVar;i++)
		{
			_AroeMinusImplicit[i] = _AroeMinus[i];
			_AroePlusImplicit[i]  = _AroePlus[i];
		}
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		displayMatrix(_Aroe, _nVar,"Matrice de Roe");
		displayMatrix(_absAroe, _nVar,"Valeur absolue matrice de Roe");
		displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus");
		displayMatrix(_AroePlus, _nVar,"Matrice _AroePlus");
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
		double u_n=0;//q_n=quantité de mouvement normale à la face frontière;
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
			//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
			Cell Cj=_mesh.getCell(j);
			double hydroPress=Cj.x()*_GravityField3d[0];
			if(_Ndim>1){
				hydroPress+=Cj.y()*_GravityField3d[1];
				if(_Ndim>2)
					hydroPress+=Cj.z()*_GravityField3d[2];
			}
			hydroPress*= _fluides[0]->getDensity(_limitField[nameOfGroup].p, _Temperature) ;//multiplication by rho the total density
			_externalStates[0] = hydroPress + _limitField[nameOfGroup].p;
	
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout<<"Cond lim outlet pressure= "<<_externalStates[0]<<" gravite= "<<_GravityField3d[0]<<" Cj.x()= "<<Cj.x()<<endl;
				cout<<"Cond lim outlet reference pressure= "<<_limitField[nameOfGroup].p<<" pression hydro= "<<hydroPress<<" total= "<<_limitField[nameOfGroup].p+hydroPress<<endl;
			}
		}
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<< " _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, InletPressure and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, InletPressure and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}

	_idm[0] = 0;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	_Vext = _externalStates;
	_Vextdiff = _externalStates;
}

void IsothermalSinglePhase::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)

{
	double mu     = _fluides[0]->getViscosity(_Temperature);

	if(isBord)
	{
		_Vj = _Vextdiff;
		_inv_dxj=_inv_dxi;
	}

	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_nVar-1; k++)
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

void IsothermalSinglePhase::sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i)
{
	double phirho=Ui[0], T=Vi[_nVar-1];
	double norm_u=0;
	for(int k=0; k<_Ndim; k++)
		norm_u+=Vi[1+k]*Vi[1+k];
	norm_u=sqrt(norm_u);

	Si[0]=0;
	for(int k=1; k<_nVar-1; k++)
		Si[k]  =(_gravite[k]-_dragCoeffs[0]*norm_u*Vi[1+k])*phirho;

	if(_timeScheme==Implicit)
	{
		for(int k=0; k<_nVar*_nVar;k++)
			_GravityImplicitationMatrix[k] = 0;
		if(!_usePrimitiveVarsInNewton)
			for(int k=0; k<_nVar;k++)
				_GravityImplicitationMatrix[k*_nVar]=-_gravite[k];
		else
		{
			double pression=Vi[0];
			getDensityDerivatives( pression);
			for(int k=0; k<_nVar;k++)
				_GravityImplicitationMatrix[k*_nVar+0]      =-_gravite[k]*_drho_sur_dp;
		}
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
		if(_timeScheme==Implicit)
			displayMatrix(_GravityImplicitationMatrix, _nVar, "Gravity implicitation matrix");
	}
}

void IsothermalSinglePhase::getDensityDerivatives( double pressure)
{
	double rho=_fluides[0]->getDensity(pressure,_Temperature);
	double gamma=_fluides[0]->constante("gamma");
	double q=_fluides[0]->constante("q");

	StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_fluides[0]);
	double e = fluide0->getInternalEnergy(_Temperature);

	_drho_sur_dp=1/((gamma-1)*(e-q));

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"_drho_sur_dp= "<<_drho_sur_dp<<endl;	
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
	_blockDiag[1+_Ndim]=1;
	_invBlockDiag[1+_Ndim]=1.0;//_blockDiag[1+_Ndim];
	for(int q=1; q<_Ndim+1; q++)
	{
		_blockDiag[q]=1/(maxvp*maxvp);
		_invBlockDiag[q]=1;//_blockDiag[q];
		_blockDiag[q+1+_Ndim]=1/(maxvp*maxvp);
		_invBlockDiag[q+1+_Ndim]=1;//_blockDiag[q+1+_Ndim];
	}
}

void IsothermalSinglePhase::jacobianConvGhostState(const int &j, string nameOfGroup,double * normale)
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
		for(k=1; k<_nVar;k++)
			for(int l=1; l<_nVar;l++)
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

		double rho_ext=_fluides[0]->getDensity(_limitField[nameOfGroup].p, _limitField[nameOfGroup].T);
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

		double rho_ext=_fluides[0]->getDensity(_limitField[nameOfGroup].p, _Temperature);
		double rho_int = _externalStates[0];
		double density_ratio=rho_ext/rho_int;
		double internal_energy=_fluides[0]->getInternalEnergy(_Temperature,rho_int);
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
		p0=_fluides[0]->constante("p0");
		gamma =_fluides[0]->constante("gamma");
		cn =_limitField[nameOfGroup].p +p0;
		cd = _phi[0]*_fluides[0]->getInternalEnergy(_externalStates[_nVar-1],rho)-p0;
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
			_JcbDiff[(1+idim)*_nVar + idim + 1] +=( cn*(_phi[0]*_fluides[0]->getInternalEnergy(_externalStates[_nVar-1],rho)-p0))/cd;

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
		throw CdmathException("IsothermalSinglePhase::jacobianConvGhostState: This boundary condition is not treated");
	}
}

void IsothermalSinglePhase::jacobianDiffGhostState(const int &j, string nameOfGroup)
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
		double internal_energy=_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho);
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
		throw CdmathException("IsothermalSinglePhase::jacobianDiffGhostState: This boundary condition is not recognised");
	}
}

void IsothermalSinglePhase::primToCons(const double *P, const int &i, double *W, const int &j)
{   //We do not know the size of Wcons and Wprim 
	//Sometimes they have _nVar components, sometimes they have _Nmailles*_nVar
	double phi_rho =_porosityField(j)*_fluides[0]->getDensity(P[i*_nVar], _Temperature);
	W[j*(_Ndim+2)] =  phi_rho;//phi*rho
	for(int k=0; k<_Ndim; k++)
		W[j*_nVar+(k+1)] = phi_rho*P[i*_nVar+(k+1)];//phi*rho*u
}

void IsothermalSinglePhase::primToConsJacobianMatrix(double *V)
{//V vecteur primitif de taille _nVar
	double pression=V[0];
	double rho=_fluides[0]->getDensity(pression,_Temperature);
	double invSoundSpeed = _fluides[0]->getInverseSoundSpeed(pression,_Temperature);
	
	_primToConsJacoMat[0] = invSoundSpeed;

	for(int idim=0;idim<_Ndim;idim++)
	{
		_primToConsJacoMat[(idim+1)*_nVar]=V[1+idim]*invSoundSpeed;
		_primToConsJacoMat[(idim+1)*_nVar+idim+1]=rho*invSoundSpeed;
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" IsothermalSinglePhase::primToConsJacobianMatrix" << endl;
		displayVector(_Vi,_nVar," _Vi " );
		cout<<" Jacobienne primToCons: " << endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobienne primToCons: ");
	}
}

void IsothermalSinglePhase::consToPrim(const double *Wcons, double* Wprim,double porosity)//To do: treat porosity
{  //Wcons and Wprim are vectors with _nVar components
	//Wcons has been extracted from the vector _conservativeVars which has _Nmailles*_nVar components
		*_runLogFile<< "IsothermalSinglePhase::consToPrim should not be used" << endl;
		_runLogFile->close();
		throw CdmathException("IsothermalSinglePhase::consToPrim should not be used");
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
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
	}
}

void IsothermalSinglePhase::staggeredRoeUpwindingMatrixPrimitiveVariables( double u_n)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//Calcul de décentrement de type décalé pour formulation Roe
	double rho = _Uroe[ 0 ];
	_drho_sur_dp = _Uroe[_nVar];
	
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*_Uroe[1+i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*_Uroe[1+i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
	}
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
	double e_int=_fluides[0]->getInternalEnergy(Temperature,rho);

	Vector F(_nVar);
	F(0)=phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F(1+i)=phirho*vitessen*vitesse(i)+pression*normale(i)*porosity;
	F(1+_Ndim)=phirho*(e_int+0.5*vitesse*vitesse+pression/rho)*vitessen;

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
			double h=0;
			_MachNumber(i)=sqrt(v2)/_fluides[0]->vitesseSonEnthalpie(h);
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
