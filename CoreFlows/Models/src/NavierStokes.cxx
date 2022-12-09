/*
 * SinglePhase.cxx
 *
 *  Created on: Sep 15, 2014
 *      Author: tn236279
 */

#include "NavierStokes.hxx"

using namespace std;

NavierStokes::NavierStokes(phaseType fluid, pressureEstimate pEstimate, int dim, bool useDellacherieEOS){
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
		_Density=Field("Density",CELLS,_mesh,1);
		_Temperature=Field("Temperature",CELLS,_mesh,1);
		_MachNumber=Field("MachNumber",CELLS,_mesh,1);
		_VitesseX=Field("Velocity x",CELLS,_mesh,1);
		if(_Ndim>1)
		{
			_VitesseY=Field("Velocity y",CELLS,_mesh,1);
			if(_Ndim>2)
				_VitesseZ=Field("Velocity z",CELLS,_mesh,1);
		}
	}

	if(_entropicCorrection)
		_entropicShift=vector<double>(3,0);//at most 3 distinct eigenvalues

	ProblemFluid::initialize();
}
//#Todo _constantPressureVector ?
void IsothermalSinglePhase::terminate(){
	delete[] _Vdiff,_Vextdiff,_Vext;
	
	if( _isSingularSystem )
		VecDestroy(&_constantPressureVector);

	ProblemFluid::terminate();
}

void NavierStokes::convectionState( const long &i, const long &j, const bool &IsBord){
	//First conservative state then further down we will compute interface (Roe) state and then compute primitive state
	_idm[0] = _nVar*i; 
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);

	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	if(IsBord)
		VecGetValues(_Uext, _nVar, _idm, _Uj);
	else
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection CONSERVATIVE Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Ui[k]<<endl;
		cout<<"Convection CONSERVATIVE Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Uj[k]<<endl;
	}
	if(_Ui[0]<0||_Uj[0]<0)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!densité negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!densité negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("densité negative, arret de calcul");
	
	_idm[0] = _nVar*i; 
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);

	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	if(IsBord)
		VecGetValues(_Vext, _nVar, _idm, _Vj);
	else
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection PRIMITIVE Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Vi[k]<<endl;
		cout<<"Convection PRIMITIVE Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Vj[k]<<endl;
	}
	if(_Vi[0]<0||_Vj[0]<0)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!pression negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		*_runLogFile<<"!!!!!!!!!!!!!!!!!!!!!!!!pression negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		_runLogFile->close();
		throw CdmathException("pression negative, arret de calcul");
	}
	
	//#Todo Etat de Roe à passer en variables primitives
	//Computation of Roe state
	PetscScalar ri, rj, xi, xj;
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
	double H_i = _compressibleFluid->getEnthalpy(_Vi[2], _Ui[0] );
	double H_j = _compressibleFluid->getEnthalpy(_Vj[2], _Uj[0] );
	_Uroe[_nVar - 1] = (ri*H_i + H_j*xj)/(ri + rj);
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout << "Enthalpie totale de Roe H  gauche " << i << ": " << xi << ", droite " << j << ": " << xj << "->" << _Uroe[_nVar-1] << endl;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Convection interfacial state"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<< _Uroe[k]<<" , "<<endl;
	}
	
}

void NavierStokes::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat de diffusion (p, u, T)
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	for(int k=0; k<_nVar; k++)
		_idn[k] = k;

	if(IsBord)
		VecGetValues(_Vextdiff, _nVar, _idn, _Vj);
	else
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "NavierStokes::diffusionStateAndMatrices primitives cellule gauche" << i << endl;
		cout << "Vi = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vi[q]  << "\t";
		cout << endl;
		cout << "NavierStokes::diffusionStateAndMatrices primitives cellule droite" << j << endl;
		cout << "Vj = ";
		for(int q=0; q<_nVar; q++)
			cout << _Vj[q]  << "\t";
		cout << endl;
	}

	// #Todo Faut-il modifier ces valeurs en variables primitives ?
	for(int k=0; k<_nVar; k++)
		_Vdiff[k] = (_Vi[k]/_porosityi+_Vj[k]/_porosityj)/2;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "SinglePhase::diffusionStateAndMatrices conservative diffusion state" << endl;
		cout << "_Udiff = ";
		for(int q=0; q<_nVar; q++)
			cout << _Udiff[q]  << "\t";
		cout << endl;
		cout << "porosite gauche= "<<_porosityi<< ", porosite droite= "<<_porosityj<<endl;
	}
	consToPrim(_Udiff,_phi,1);
	_Udiff[_nVar-1]=_phi[_nVar-1];

	if(_timeScheme==Implicit)
	{
		double q_2=0;
		for (int i = 0; i<_Ndim;i++)
			q_2+=_Udiff[i+1]*_Udiff[i+1];
		double mu = _compressibleFluid->getViscosity(_Udiff[_nVar-1]);
		double lambda = _compressibleFluid->getConductivity(_Udiff[_nVar-1]);
		double Cv= _compressibleFluid->constante("Cv");
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		for(int i=1;i<(_nVar-1);i++)
		{
			_Diffusion[i*_nVar] =  mu*_Udiff[i]/(_Udiff[0]*_Udiff[0]);
			_Diffusion[i*_nVar+i] = -mu/_Udiff[0];
		}
		int i = (_nVar-1)*_nVar;
		_Diffusion[i]=lambda*(_Udiff[_nVar-1]/_Udiff[0]-q_2/(2*Cv*_Udiff[0]*_Udiff[0]*_Udiff[0]));
		for(int k=1;k<(_nVar-1);k++)
		{
			_Diffusion[i+k]= lambda*_Udiff[k]/(_Udiff[0]*_Udiff[0]*Cv);
		}
		_Diffusion[_nVar*_nVar-1]=-lambda/(_Udiff[0]*Cv);
	}
}
void NavierStokes::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	double porosityj=_porosityField(j);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "setBoundaryState for group "<< nameOfGroup << ", inner cell j= "<<j<< " face unit normal vector "<<endl;
		for(int k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}

	if (_limitField[nameOfGroup].bcType==Wall){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne conservatif
		double q_n=0;//q_n=quantité de mouvement normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_n+=_externalStates[(k+1)]*normale[k];
			
		//Pour la convection, inversion du sens de la vitesse normale
		for(int k=0; k<_Ndim; k++)
			_externalStates[(k+1)]-= 2*q_n*normale[k];

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		//Pour la diffusion, paroi à vitesse et temperature imposees
		_idm[0] = _nVar*j;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//L'état fantome contient à présent les variables primitives internes
		double pression=_externalStates[0];
		double T=_limitField[nameOfGroup].T;
		double rho=_compressibleFluid->getDensity(pression,T);

		_externalStates[0]=porosityj*rho;
		_externalStates[1]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
		double v2=0;
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
		_externalStates[_nVar-1] = _externalStates[0]*(_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On prend l'état fantôme égal à l'état interne (conditions limites de Neumann)

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne (conditions limites de Neumann)
		double q_int_n=0;//q_int_n=composante normale de la quantité de mouvement à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_int_n+=_externalStates[(k+1)]*normale[k];//On calcule la vitesse normale sortante

		double q_ext_n=_limitField[nameOfGroup].v_x[0]*normale[0];
		if(_Ndim>1)
			{
				q_ext_n+=_limitField[nameOfGroup].v_y[0]*normale[1];
				if(_Ndim>2)
						q_ext_n+=_limitField[nameOfGroup].v_z[0]*normale[2];
			}

		if(q_int_n+q_ext_n<=0){//Interfacial velocity goes inward
			VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);//On met à jour l'état fantome avec les variables primitives internes
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_compressibleFluid->getDensity(pression,T);

			_externalStates[0]=porosityj*rho;//Composante fantome de masse
			_externalStates[1]=_externalStates[0]*(_limitField[nameOfGroup].v_x[0]);//Composante fantome de qdm x
			double v2=0;
			v2 +=(_limitField[nameOfGroup].v_x[0])*(_limitField[nameOfGroup].v_x[0]);
			if(_Ndim>1)
			{
				v2 +=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
				_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];//Composante fantome de qdm y
				if(_Ndim==3)
				{
					_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];//Composante fantome de qdm z
					v2 +=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[0]*(_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);//Composante fantome de l'nrj
		}
		else if(_nbTimeStep%_freqSave ==0)
			cout<< "Warning : fluid possibly going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletRotationVelocity){
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double u_int_n=0;//u_int_n=composante normale de la vitesse intérieure à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_int_n+=_externalStates[(k+1)]*normale[k];

		double u_ext_n=0;
        //ghost velocity
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
			u_ext_n=-0.01*tangent_vel.norm();
			//Changing external state velocity
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]=u_ext_n*normale[k] + tangent_vel[k];
        }

		if(u_ext_n + u_int_n <= 0){
			double pression=_externalStates[0];
			double T=_limitField[nameOfGroup].T;
			double rho=_compressibleFluid->getDensity(pression,T);

			double v2=0;
			v2 +=_externalStates[1]*_externalStates[1];//v_x*v_x
			_externalStates[0]=porosityj*rho;
			_externalStates[1]*=_externalStates[0];
			if(_Ndim>1)
			{
				v2 +=_externalStates[2]*_externalStates[2];//+v_y*v_y
				_externalStates[2]*=_externalStates[0];
				if(_Ndim==3)
				{
					v2 +=_externalStates[3]*_externalStates[3];//+v_z*v_z
					_externalStates[3]*=_externalStates[0];
				}
			}
			_externalStates[_nVar-1] = _externalStates[0]*(_compressibleFluid->getInternalEnergy(_limitField[nameOfGroup].T,rho) + v2/2);
		}
		else if(_nbTimeStep%_freqSave ==0)
		{
			/*
			 * cout<< "Warning : fluid going out through inlet boundary "<<nameOfGroup<<". Applying Neumann boundary condition"<<endl;
			 */ 
			VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On définit l'état fantôme avec l'état interne
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Wall boundary condition."<<endl;
			
			//Changing external state momentum
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]-=2*_externalStates[0]*u_int_n*normale[k];
		}

		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density

		//Building the primitive external state
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		double u_n=0;//u_n=vitesse normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			u_n+=_externalStates[(k+1)]*normale[k];
        
		if(u_n<=0)
		{
			_externalStates[0]=porosityj*_compressibleFluid->getDensity(_limitField[nameOfGroup].p+hydroPress,_limitField[nameOfGroup].T);
			
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
	                _externalStates[(k+1)]=u_n*normale[k] + abs(u_n)*tangent_vel[k];
	        }
		}
		else{
			/*
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Neumann boundary condition for velocity and temperature (only pressure value is imposed as in outlet BC)."<<endl;
			_externalStates[0]=porosityj*_compressibleFluid->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
			*/
			if(_nbTimeStep%_freqSave ==0)
				cout<< "Warning : fluid going out through inletPressure boundary "<<nameOfGroup<<". Applying Wall boundary condition."<<endl;
			_externalStates[0]=porosityj*_compressibleFluid->getDensity(_externalStates[0]+hydroPress, _externalStates[_nVar-1]);
			//Changing external state velocity
            for(int k=0; k<_Ndim; k++)
                _externalStates[(k+1)]-=2*u_n*normale[k];
		}

		double v2=0;
		for(int k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;//qdm component
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_compressibleFluid->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);//nrj component


		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne conservatif
		double q_n=0;//q_n=quantité de mouvement normale à la face frontière;
		for(int k=0; k<_Ndim; k++)
			q_n+=_externalStates[(k+1)]*normale[k];

		if(q_n < -_precision &&  _nbTimeStep%_freqSave ==0)
        {
			cout<< "Warning : fluid going in through outlet boundary "<<nameOfGroup<<" with flow rate "<< q_n<<endl;
            cout<< "Applying Neumann boundary condition for velocity and temperature"<<endl;
        }

		double hydroPress=0;
		if( _VV.getSpaceDimension()>1)
		{
			//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
			Cell Cj=_mesh.getCell(j);
			hydroPress=(Cj.x()-_gravityReferencePoint[0])*_GravityField3d[0];
			if(_Ndim>1){
				hydroPress+=(Cj.y()-_gravityReferencePoint[1])*_GravityField3d[1];
				if(_Ndim>2)
					hydroPress+(Cj.z()-_gravityReferencePoint[2])*_GravityField3d[2];
			}
			hydroPress*=_externalStates[0]/porosityj;//multiplication by rho the total density
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout<<"Cond lim outlet densite= "<<_externalStates[0]<<" gravite= "<<_GravityField3d[0]<<" Cj.x()= "<<Cj.x()<<endl;
				cout<<"Cond lim outlet pression ref= "<<_limitField[nameOfGroup].p<<" pression hydro= "<<hydroPress<<" total= "<<_limitField[nameOfGroup].p+hydroPress<<endl;
			}
		}

		//Building the external state
		_idm[0] = _nVar*j;// Kieu
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);

		_externalStates[0]=porosityj*_compressibleFluid->getDensity(_limitField[nameOfGroup].p+hydroPress, _externalStates[_nVar-1]);
		double v2=0;
		for(int k=0; k<_Ndim; k++)
		{
			v2+=_externalStates[(k+1)]*_externalStates[(k+1)];
			_externalStates[(k+1)]*=_externalStates[0] ;
		}
		_externalStates[_nVar-1] = _externalStates[0]*(_compressibleFluid->getInternalEnergy( _externalStates[_nVar-1],_externalStates[0]) + v2/2);
		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);
	}else {
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<< " _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		*_runLogFile<<"Boundary condition not set for boundary named. Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		_runLogFile->close();
		throw CdmathException("Unknown boundary condition");
	}
}

void NavierStokes::convectionMatrices()
{
	//entree: URoe = p, u, T
	//sortie: matrices Roe+  et Roe-

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"NavierStokes::convectionMatrices()"<<endl;

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
		c = _compressibleFluid->vitesseSonEnthalpie(H-u_2/2);//vitesse du son a l'interface
		k = _compressibleFluid->constante("gamma") - 1;//A generaliser pour porosite et stephane gas law
		K = u_2*k/2; //g-1/2 *|u|²

		vp_dist[0]=u_n-c;vp_dist[1]=u_n;vp_dist[2]=u_n+c;

		_maxvploc=fabs(u_n)+c;
		if(_maxvploc>_maxvp)
			_maxvp=_maxvploc;

		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<"SinglePhase::convectionMatrices Eigenvalues "<<u_n-c<<" , "<<u_n<<" , "<<u_n+c<<endl;

		RoeMatrixConservativeVariables( u_n, H,vitesse,k,K);

		/******** Construction des matrices de decentrement ********/
		if( _spaceScheme ==centered){
			if(_entropicCorrection)
			{
				*_runLogFile<<"SinglePhase::convectionMatrices: entropy scheme not available for centered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhase::convectionMatrices: entropy scheme not available for centered scheme");
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
				*_runLogFile<<"SinglePhase::convectionMatrices: entropy scheme not available for staggered scheme"<<endl;
				_runLogFile->close();
				throw CdmathException("SinglePhase::convectionMatrices: entropy scheme not available for staggered scheme");
			}

			staggeredRoeUpwindingMatrixConservativeVariables( u_n, H, vitesse, k, K);
		}
		else
		{
			*_runLogFile<<"SinglePhase::convectionMatrices: scheme not treated"<<endl;
			_runLogFile->close();
			throw CdmathException("SinglePhase::convectionMatrices: scheme not treated");
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
				_Vij[0]=_compressibleFluid->getPressureFromEnthalpy(_Uroe[_nVar-1]-u_2/2, _Uroe[0]);//pressure
				_Vij[_nVar-1]=_compressibleFluid->getTemperatureFromPressure( _Vij[0], _Uroe[0]);//Temperature
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
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			displayMatrix(_Aroe, _nVar,"Matrice de Roe");
			displayMatrix(_absAroe, _nVar,"Valeur absolue matrice de Roe");
			displayMatrix(_AroeMinus, _nVar,"Matrice _AroeMinus");
			displayMatrix(_AroePlus, _nVar,"Matrice _AroePlus");
		}
	}

	if(_verbose && _nbTimeStep%_freqSave ==0 && _timeScheme==Implicit)
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
		*_runLogFile<<"SinglePhase::convectionMatrices: well balanced option not treated"<<endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::convectionMatrices: well balanced option not treated");
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

void NavierStokes::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)

{
	double lambda = _compressibleFluid->getConductivity(_Udiff[_nVar-1]);
	double mu     = _compressibleFluid->getViscosity(_Udiff[_nVar-1]);

	if(isBord )
		lambda=max(lambda,_heatTransfertCoeff);//wall nucleate boing -> larger heat transfer

	if(lambda==0 && mu ==0 && _heatTransfertCoeff==0)
		return;

	_idm[0] = 0;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	if(isBord)
	{
		VecGetValues(_Uextdiff, _nVar, _idm, _Uj);
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

void NavierStokes::sourceVector(PetscScalar * Si, PetscScalar * Ui, PetscScalar * Vi, int i)
{
	double phirho=Ui[0], T=Vi[_nVar-1];
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
			double pression=Vi[0];
			getDensityDerivatives( pression, T, norm_u*norm_u);
			for(int k=0; k<_nVar;k++)
			{
				_GravityImplicitationMatrix[k*_nVar+0]      =-_gravite[k]*_drho_sur_dp;
				_GravityImplicitationMatrix[k*_nVar+_nVar-1]=-_gravite[k]*_drho_sur_dT;
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

	assert( P[i*(_Ndim+2)]>0);//Pressure should be positive
	assert( P[i*(_Ndim+2)+(_Ndim+1)]>0);//Temperature should be positive

	double rho=_compressibleFluid->getDensity(P[i*(_Ndim+2)], P[i*(_Ndim+2)+(_Ndim+1)]);
	W[j*(_Ndim+2)] =  _porosityField(j)*rho;//phi*rho
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(k+1)] = W[j*(_Ndim+2)]*P[i*(_Ndim+2)+(k+1)];//phi*rho*u

	W[j*(_Ndim+2)+(_Ndim+1)] = W[j*(_Ndim+2)]*_compressibleFluid->getInternalEnergy(P[i*(_Ndim+2)+ (_Ndim+1)],rho);//rho*e
	for(int k=0; k<_Ndim; k++)
		W[j*(_Ndim+2)+(_Ndim+1)] += W[j*(_Ndim+2)]*P[i*(_Ndim+2)+(k+1)]*P[i*(_Ndim+2)+(k+1)]*0.5;//phi*rho*e+0.5*phi*rho*u^2
}

void NavierStokes::primToConsJacobianMatrix(double *V)
{
	double pression=V[0];
	double temperature=V[_nVar-1];
	double vitesse[_Ndim];
	for(int idim=0;idim<_Ndim;idim++)
		vitesse[idim]=V[1+idim];
	double v2=0;
	for(int idim=0;idim<_Ndim;idim++)
		v2+=vitesse[idim]*vitesse[idim];

	double rho=_compressibleFluid->getDensity(pression,temperature);
	double gamma=_compressibleFluid->constante("gamma");
	double cv=_compressibleFluid->constante("cv");

	for(int k=0;k<_nVar*_nVar; k++)
		_primToConsJacoMat[k]=0;

	if(		!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_compressibleFluid);
		double e=fluide0->getInternalEnergy(temperature);
		double E=e+0.5*v2;
		/* To do : replace the formulas usind p0 and q by calls to sound speed */
		double Pinf = fluide0->constante("p0");
		double    q = fluide0->constante("q");

		_primToConsJacoMat[0]=1/((gamma-1)*(e-q));
		_primToConsJacoMat[_nVar-1]=-rho*cv/(e-q);

		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[_nVar+idim*_nVar]=vitesse[idim]/((gamma-1)*(e-q));
			_primToConsJacoMat[_nVar+idim*_nVar+1+idim]=rho;
			_primToConsJacoMat[_nVar+idim*_nVar+_nVar-1]=-rho*vitesse[idim]*cv/(e-q);
		}
		_primToConsJacoMat[(_nVar-1)*_nVar]=E/((gamma-1)*(e-q));
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+1+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*cv*(1-E/(e-q));
	}
	else if(	_useDellacherieEOS)
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_compressibleFluid);
		double h=fluide0->getEnthalpy(temperature);
		double H=h+0.5*v2;
		double cp=_compressibleFluid->constante("cp");
		/* To do : replace the formulas usind p0 and q by calls to sound speed */
		double Pinf = fluide0->constante("p0");
		double    q = fluide0->constante("q");

		_primToConsJacoMat[0]=gamma/((gamma-1)*(h-q));
		_primToConsJacoMat[_nVar-1]=-rho*cp/(h-q);

		for(int idim=0;idim<_Ndim;idim++)
		{
			_primToConsJacoMat[_nVar+idim*_nVar]=gamma*vitesse[idim]/((gamma-1)*(h-q));
			_primToConsJacoMat[_nVar+idim*_nVar+1+idim]=rho;
			_primToConsJacoMat[_nVar+idim*_nVar+_nVar-1]=-rho*vitesse[idim]*cp/(h-q);
		}
		_primToConsJacoMat[(_nVar-1)*_nVar]=gamma*H/((gamma-1)*(h-q))-1;
		for(int idim=0;idim<_Ndim;idim++)
			_primToConsJacoMat[(_nVar-1)*_nVar+1+idim]=rho*vitesse[idim];
		_primToConsJacoMat[(_nVar-1)*_nVar+_nVar-1]=rho*cp*(1-H/(h-q));
	}
	else
	{
		*_runLogFile<<"SinglePhase::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie"<<endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::primToConsJacobianMatrix: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" SinglePhase::primToConsJacobianMatrix" << endl;
		displayVector(_Vi,_nVar," _Vi " );
		cout<<" Jacobienne primToCons: " << endl;
		displayMatrix(_primToConsJacoMat,_nVar," Jacobienne primToCons: ");
	}
}

void NavierStokes::consToPrim(const double *Wcons, double* Wprim,double porosity)
{
	assert( Wcons[0]>0);//Density should be positive
	assert( Wcons[_nVar-1]>0);//Total energy should be positive

	double q_2 = 0;
	for(int k=1;k<=_Ndim;k++)
		q_2 += Wcons[k]*Wcons[k];
	q_2 /= Wcons[0];	//phi rho u²
	double rhoe=(Wcons[(_Ndim+2)-1]-q_2/2)/porosity;
	double rho=Wcons[0]/porosity;
	Wprim[0] =  _compressibleFluid->getPressure(rhoe,rho);//pressure p
	if (Wprim[0]<0){
		cout << "pressure = "<< Wprim[0] << " < 0 " << endl;
		*_runLogFile<< "pressure = "<< Wprim[0] << " < 0 " << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::consToPrim: negative pressure");
	}
	for(int k=1;k<=_Ndim;k++)
		Wprim[k] = Wcons[k]/Wcons[0];//velocity u
	Wprim[(_Ndim+2)-1] =  _compressibleFluid->getTemperatureFromPressure(Wprim[0],Wcons[0]/porosity);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"ConsToPrim Vecteur conservatif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Wcons[k]<<endl;
		cout<<"ConsToPrim Vecteur primitif"<<endl;
		for(int k=0;k<_nVar;k++)
			cout<<Wprim[k]<<endl;
	}
}

void NavierStokes::entropicShift(double* n)//TO do: make sure _Vi and _Vj are well set
{
	/*Left and right values */
	double ul_n = 0, ul_2=0, ur_n=0,	ur_2 = 0; //valeurs de vitesse normale et |u|² a droite et a gauche
	for(int i=0;i<_Ndim;i++)
	{
		ul_n += _Vi[1+i]*n[i];
		ul_2 += _Vi[1+i]*_Vi[1+i];
		ur_n += _Vj[1+i]*n[i];
		ur_2 += _Vj[1+i]*_Vj[1+i];
	}

	double cl = _compressibleFluid->vitesseSonEnthalpie(_Vi[_Ndim+1]-ul_2/2);//vitesse du son a l'interface
	double cr = _compressibleFluid->vitesseSonEnthalpie(_Vj[_Ndim+1]-ur_2/2);//vitesse du son a l'interface

	_entropicShift[0]=fabs(ul_n-cl - (ur_n-cr));
	_entropicShift[1]=fabs(ul_n     - ur_n);
	_entropicShift[2]=fabs(ul_n+cl - (ur_n+cr));

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Entropic shift left eigenvalues: "<<endl;
		cout<<"("<< ul_n-cl<< ", " <<ul_n<<", "<<ul_n+cl << ")";
		cout<<endl;
		cout<<"Entropic shift right eigenvalues: "<<endl;
		cout<<"("<< ur_n-cr<< ", " <<ur_n<<", "<<ur_n+cr << ")";
		cout<< endl;
		cout<<"eigenvalue jumps "<<endl;
		cout<< _entropicShift[0] << ", " << _entropicShift[1] << ", "<< _entropicShift[2] <<endl;
	}
}

Vector NavierStokes::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::convectionFlux start"<<endl;
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
	double e_int=_compressibleFluid->getInternalEnergy(Temperature,rho);

	Vector F(_nVar);
	F(0)=phirho*vitessen;
	for(int i=0;i<_Ndim;i++)
		F(1+i)=phirho*vitessen*vitesse(i)+pression*normale(i)*porosity;
	F(1+_Ndim)=phirho*(e_int+0.5*vitesse*vitesse+pression/rho)*vitessen;

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"SinglePhase::convectionFlux end"<<endl;
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void NavierStokes::RoeMatrixConservativeVariables(double u_n, double H,Vector velocity, double k, double K)
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
void NavierStokes::convectionMatrixPrimitiveVariables( double rho, double u_n, double H,Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//On remplit la matrice de Roe en variables primitives : F(V_L)-F(V_R)=Aroe (V_L-V_R)
	//EOS is more involved with primitive variables
	// call to getDensityDerivatives(double concentration, double pression, double temperature,double v2) needed
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim]=_drho_sur_dT*u_n;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*vitesse[i]+_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*vitesse[i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim]=_drho_sur_dT*u_n*vitesse[i];
	}
	_AroeImplicit[(1+_Ndim)*_nVar+0]=(_drhoE_sur_dp+1)*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(1+_Ndim)*_nVar+1+i]=rho*(H*_vec_normal[i]+u_n*vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim]=_drhoE_sur_dT*u_n;
}
void NavierStokes::staggeredRoeUpwindingMatrixConservativeVariables( double u_n, double H,Vector velocity, double k, double K)
{
	//Calcul de décentrement de type décalé pour formulation de Roe
	if(fabs(u_n)>_precision)
	{
		//premiere ligne (masse)
		_absAroe[0]=0;
		for(int idim=0; idim<_Ndim;idim++)
			_absAroe[1+idim]=_vec_normal[idim];
		_absAroe[_nVar-1]=0;

		//lignes intermadiaires(qdm)
		for(int idim=0; idim<_Ndim;idim++)
		{
			//premiere colonne
			_absAroe[(1+idim)*_nVar]=-K*_vec_normal[idim] - u_n*_Uroe[1+idim];
			//colonnes intermediaires
			for(int jdim=0; jdim<_Ndim;jdim++)
				_absAroe[(1+idim)*_nVar + jdim + 1] = _Uroe[1+idim]*_vec_normal[jdim]+k*_vec_normal[idim]*_Uroe[1+jdim];
			//matrice identite
			_absAroe[(1+idim)*_nVar + idim + 1] += u_n;
			//derniere colonne
			_absAroe[(1+idim)*_nVar + _nVar-1]=-k*_vec_normal[idim];
		}

		//derniere ligne (energie)
		_absAroe[_nVar*(_nVar-1)] = (-K - H)*u_n;
		for(int idim=0; idim<_Ndim;idim++)
			_absAroe[_nVar*(_nVar-1)+idim+1]=H*_vec_normal[idim] + k*u_n*_Uroe[idim+1];
		_absAroe[_nVar*_nVar -1] = (1 - k)*u_n;

		double signu;
		if(u_n>0)
			signu=1;
		else if (u_n<0)
			signu=-1;

		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] *= signu;
	}
	else//umn=0 ->centered scheme
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i] =0;
	}
}
void NavierStokes::staggeredRoeUpwindingMatrixPrimitiveVariables(double rho, double u_n,double H, Vector vitesse)
{
	//Not used. Suppress or use in alternative implicitation in primitive variable of the staggered-roe scheme
	//Calcul de décentrement de type décalé pour formulation Roe
	_AroeImplicit[0*_nVar+0]=_drho_sur_dp*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[0*_nVar+1+i]=rho*_vec_normal[i];
	_AroeImplicit[0*_nVar+1+_Ndim]=_drho_sur_dT*u_n;
	for(int i=0;i<_Ndim;i++)
	{
		_AroeImplicit[(1+i)*_nVar+0]=_drho_sur_dp *u_n*vitesse[i]-_vec_normal[i];
		for(int j=0;j<_Ndim;j++)
			_AroeImplicit[(1+i)*_nVar+1+j]=rho*vitesse[i]*_vec_normal[j];
		_AroeImplicit[(1+i)*_nVar+1+i]+=rho*u_n;
		_AroeImplicit[(1+i)*_nVar+1+_Ndim]=_drho_sur_dT*u_n*vitesse[i];
	}
	_AroeImplicit[(1+_Ndim)*_nVar+0]=(_drhoE_sur_dp+1)*u_n;
	for(int i=0;i<_Ndim;i++)
		_AroeImplicit[(1+_Ndim)*_nVar+1+i]=rho*(H*_vec_normal[i]+u_n*vitesse[i]);
	_AroeImplicit[(1+_Ndim)*_nVar+1+_Ndim]=_drhoE_sur_dT*u_n;
}

Vector NavierStokes::staggeredVFFCFlux()
{
	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout<<"SinglePhase::staggeredVFFCFlux() start"<<endl;

	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
	{
		*_runLogFile<< "SinglePhase::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding, pressure = "<<  endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding");
	}
	else//_spaceScheme==staggered && _nonLinearFormulation==VFFC
	{
		Vector Fij(_nVar);

		double uijn=0, phiqn=0, uin=0, ujn=0;
		for(int idim=0;idim<_Ndim;idim++)
		{
			uijn+=_vec_normal[idim]*_Uroe[1+idim];//URoe = rho, u, H
			uin +=_vec_normal[idim]*_Ui[1+idim];
			ujn +=_vec_normal[idim]*_Uj[1+idim];
		}
		
		if( (uin>0 && ujn >0) || (uin>=0 && ujn <=0 && uijn>0) ) // formerly (uijn>_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Ui[1+idim];//phi rho u n
			Fij(0)=phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phiqn*_Vi[1+idim]+_Vj[0]*_vec_normal[idim]*_porosityj;
			Fij(_nVar-1)=phiqn/_Ui[0]*(_Ui[_nVar-1]+_Vj[0]*sqrt(_porosityj/_porosityi));
		}
		else if( (uin<0 && ujn <0) || (uin>=0 && ujn <=0 && uijn<0) ) // formerly (uijn<-_precision)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phiqn+=_vec_normal[idim]*_Uj[1+idim];//phi rho u n
			Fij(0)=phiqn;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phiqn*_Vj[1+idim]+_Vi[0]*_vec_normal[idim]*_porosityi;
			Fij(_nVar-1)=phiqn/_Uj[0]*(_Uj[_nVar-1]+_Vi[0]*sqrt(_porosityi/_porosityj));
		}
		else//case (uin<=0 && ujn >=0) or (uin>=0 && ujn <=0 && uijn==0), apply centered scheme
		{
			Vector Ui(_nVar), Uj(_nVar), Vi(_nVar), Vj(_nVar), Fi(_nVar), Fj(_nVar);
			Vector normale(_Ndim);
			for(int i1=0;i1<_Ndim;i1++)
				normale(i1)=_vec_normal[i1];
			for(int i1=0;i1<_nVar;i1++)
			{
				Ui(i1)=_Ui[i1];
				Uj(i1)=_Uj[i1];
				Vi(i1)=_Vi[i1];
				Vj(i1)=_Vj[i1];
			}
			Fi=convectionFlux(Ui,Vi,normale,_porosityi);
			Fj=convectionFlux(Uj,Vj,normale,_porosityj);
			Fij=(Fi+Fj)/2;//+_maxvploc*(Ui-Uj)/2;
		}
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"SinglePhase::staggeredVFFCFlux() endf uijn="<<uijn<<endl;
			cout<<Fij<<endl;
		}
		return Fij;
	}
}



void SinglePhase::testConservation()
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

void SinglePhase::getDensityDerivatives( double pressure, double temperature, double v2)
{
	double rho=_compressibleFluid->getDensity(pressure,temperature);
	double gamma=_compressibleFluid->constante("gamma");

	if(	!_useDellacherieEOS)
	{
		StiffenedGas* fluide0=dynamic_cast<StiffenedGas*>(_compressibleFluid);
		double e = fluide0->getInternalEnergy(temperature);
		double cv=fluide0->constante("cv");
		double E=e+0.5*v2;
		/* To do : replace the formula using q by calls to sound speed */
		double q = fluide0->constante("q");

		_drho_sur_dp=1/((gamma-1)*(e-q));
		_drho_sur_dT=-rho*cv/(e-q);
		_drhoE_sur_dp=E/((gamma-1)*(e-q));
		_drhoE_sur_dT=rho*cv*(1-E/(e-q));
	}
	else if(_useDellacherieEOS )
	{
		StiffenedGasDellacherie* fluide0=dynamic_cast<StiffenedGasDellacherie*>(_compressibleFluid);
		double h=fluide0->getEnthalpy(temperature);
		double H=h+0.5*v2;
		double cp=fluide0->constante("cp");
		/* To do : replace the formula using q by calls to sound speed */
		double q = fluide0->constante("q");

		_drho_sur_dp=gamma/((gamma-1)*(h-q));
		_drho_sur_dT=-rho*cp/(h-q);
		_drhoE_sur_dp=gamma*H/((gamma-1)*(h-q))-1;
		_drhoE_sur_dT=rho*cp*(1-H/(h-q));
	}
	else
	{
		*_runLogFile<< "SinglePhase::getDensityDerivatives: eos should be StiffenedGas or StiffenedGasDellacherie" << endl;
		_runLogFile->close();
		throw CdmathException("SinglePhase::getDensityDerivatives: eos should be StiffenedGas or StiffenedGasDellacherie");
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"_drho_sur_dp= "<<_drho_sur_dp<<", _drho_sur_dT= "<<_drho_sur_dT<<endl;
		cout<<"_drhoE_sur_dp= "<<_drhoE_sur_dp<<", _drhoE_sur_dT= "<<_drhoE_sur_dT<<endl;
	}
}
void SinglePhase::save(){
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

Field& SinglePhase::getPressureField()
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

Field& SinglePhase::getTemperatureField()
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

Field& SinglePhase::getVelocityField()
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

Field& SinglePhase::getMachNumberField()
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

Field& SinglePhase::getVelocityXField()
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

Field& SinglePhase::getVelocityYField()
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

Field& SinglePhase::getVelocityZField()
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

Field& SinglePhase::getDensityField()
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

Field& SinglePhase::getMomentumField()//not yet managed by parameter _saveAllFields
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

Field& SinglePhase::getTotalEnergyField()//not yet managed by parameter _saveAllFields
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

Field& SinglePhase::getEnthalpyField()
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

vector<string> SinglePhase::getOutputFieldsNames()
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

Field& SinglePhase::getOutputField(const string& nameField )
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
