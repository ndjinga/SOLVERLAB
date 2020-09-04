/*
 * IsothermalTwoFluid.cxx
 *
 *  Created on: Sep 16, 2014
 *      Author: tn236279
 */

#include "IsothermalTwoFluid.hxx"

using namespace std;

IsothermalTwoFluid::IsothermalTwoFluid(pressureEstimate pEstimate, int dim){
	_Ndim=dim;
	_nVar=2*(_Ndim+1);
	_nbPhases = 2;
	_dragCoeffs=vector<double>(2,0);
	_fluides.resize(2);
	if (pEstimate==around1bar300K)//EOS at 1 bar and 300K
	{
		cout<<"Fluid is air-water mixture around 1 bar and 300 K (27°C)"<<endl;
		*_runLogFile<<"Fluid is air-water mixture around 1 bar and 300 K (27°C)"<<endl;
		_Temperature=300;//Constant temperature of the model
		_internalEnergy1=2.22e5;//nitrogen internal energy at 1bar, 300K
		_internalEnergy2=1.12e5;//water internal energy at 1 bar, 300K
		_fluides[0] = new StiffenedGas(1.4,743,_Temperature,_internalEnergy1);  //ideal gas law for nitrogen at pressure 1 bar and temperature 27°C, c_v=743
		_fluides[1] = new StiffenedGas(996,1e5,_Temperature,_internalEnergy2,1501,4130);  //stiffened gas law for water at pressure 1 bar and temperature 27°C
	}
	else//EOS at 155 bars and 618K
	{
		cout<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		*_runLogFile<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		_Temperature=618;//Constant temperature of the model
		_internalEnergy1=2.44e6;//Gas internal energy at saturation at 155 bar
		_internalEnergy2=1.6e6;//water internal energy at saturation at 155 bar
		_fluides[0] = new StiffenedGas(102,1.55e7,_Temperature,_internalEnergy1, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C:
		_fluides[1] = new StiffenedGas(594,1.55e7,_Temperature,_internalEnergy2, 621,3100);  //stiffened gas law for water at pressure 155 bar and temperature 345°C:
	}
	_intPressCoeff=1.5;
}

void IsothermalTwoFluid::initialize(){
	cout<<"Initialising the isothermal two-fluid model"<<endl;
	*_runLogFile<<"Initialising the isothermal two-fluid model"<<endl;

	_Uroe = new double[_nVar+1];

	_guessalpha = _VV(0,0);

	_gravite = vector<double>(_nVar,0);//Not to be confused with _GravityField3d (size _Ndim). _gravite (size _Nvar) is usefull for dealing with source term and implicitation of gravity vector
	for(int i=0; i<_Ndim; i++)
	{
		_gravite[i+1]=_GravityField3d[i];
		_gravite[i+1 +_Ndim+1]=_GravityField3d[i];
	}
	_GravityImplicitationMatrix = new PetscScalar[_nVar*_nVar];

	if(_saveVelocity){
		_Vitesse1=Field("Gas velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant
		_Vitesse2=Field("Liquid velocity",CELLS,_mesh,3);//Forcement en dimension 3 pour le posttraitement des lignes de courant
	}
	if(_entropicCorrection)
		_entropicShift=vector<double>(3,0);

	ProblemFluid::initialize();
}

void IsothermalTwoFluid::convectionState( const long &i, const long &j, const bool &IsBord){
	//sortie: WRoe en (alpha, p, u1, u2, dm1,dm2,dalpha1,dp)z
	//entree: _conservativeVars en (alpha1 rho1, alpha1 rho1 u1, alpha2 rho2, alpha2 rho2 u2)

	// _l always inside the domain (index i)
	// _r is maybe the boundary cell (negative index)
	// _l and _r are primative vectors (alp, P, u1, u2)

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
		cout<<"Convection Left state cell " << i<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Ui[k]<<endl;
		cout<<"Convection Right state cell " << j<< ": "<<endl;
		for(int k =0; k<_nVar; k++)
			cout<< _Uj[k]<<endl;
	}
	//       if(_Ui[0]<-(_precision) || _Uj[0]<-(_precision) || _Ui[_Ndim+1]<-(_precision) || _Uj[_Ndim+1]<-(_precision))
	// 	{
	// 	  cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Masse partielle negative, arret de calcul!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	// 	  cout<< "valeurs a gauche: "<<_Ui[0]<<", "<<_Ui[_Ndim+1]<<"valeurs a droite: "<<_Uj[0]<<", "<<_Uj[_Ndim+1]<<endl;
	// 	  throw CdmathException(" Masse partielle negative, arret de calcul");
	// 	}
	//       else
	{
		_Ui[0]=max(0.,_Ui[0]);// mass1 a gauche
		_Uj[0]=max(0.,_Uj[0]);// mass1 a droite
		_Ui[_Ndim+1]=max(0.,_Ui[_Ndim+1]);// mass2 a gauche
		_Uj[_Ndim+1]=max(0.,_Uj[_Ndim+1]);// mass2 a droite
	}

	PetscScalar ri1, ri2, rj1, rj2, xi, xj;
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_primitiveVars, _nVar, _idm, _l);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Etat de Roe, etat primitif gauche: "<<endl;
		for(int i =0; i<_nVar; i++)
			cout<< _l[i]<<endl;
	}

	// boundary : compute the _r
	if(IsBord)
	{
		_guessalpha=_l[0];
		consToPrim(_Uj, _r);
	}
	// inside the domain : extract from the primative vector
	else
	{
		_idm[0] = _nVar*j;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _r);
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Etat de Roe, etat primitif droite: "<<endl;
		for(int i =0; i<_nVar; i++)
			cout<< _r[i]<<endl;
	}

	// Using Toumi linearisation (read in the article of Toumi) : alpha^{Roe}_l
	if(2-_l[0]-_r[0] > _precision)
		_Uroe[0] = 1- 2*(1-_l[0])*(1-_r[0])/(2-_l[0]-_r[0]);
	// Using an average : (alp_l + alp_r)/2 : suggestion of Michael (no theory)
	else
		_Uroe[0] = (_l[0]+_r[0])/2;
	// Pressure is computed as function of alp and P_l, P_r (Toumi article)
	if(_l[0]+_r[0] > _precision)
		_Uroe[1] = (_l[1]*_l[0]+_r[1]*_r[0])/(_l[0]+_r[0]);
	else
		_Uroe[1] = (_l[1]*(1-_l[0])+_r[1]*(1-_r[0]))/(2-_l[0]-_r[0]);
	// i :left, j : right (U is normally conservative variable)
	ri1 = sqrt(_Ui[0]); ri2 = sqrt(_Ui[_Ndim+1]);
	rj1 = sqrt(_Uj[0]); rj2 = sqrt(_Uj[_Ndim+1]);
	// Roe average formula of the velocities
	for(int k=0;k<_Ndim;k++)
	{
		xi = _Ui[k+1];
		xj = _Uj[k+1];
		// avoid dividing by zero, if mass is zero, do not consider the distribution of such a phase
		if(ri1>_precision && rj1>_precision)
			_Uroe[2+k] = (xi/ri1 + xj/rj1)/(ri1 + rj1);
		else if(ri1<_precision && rj1>_precision)
			_Uroe[2+k] =  xj/_Uj[0];
		else if(ri1>_precision && rj1<_precision)
			_Uroe[2+k] =  xi/_Ui[0];
		else
			_Uroe[2+k] =(_Ui[k+1+_Ndim+1]/ri2 + _Uj[k+1+_Ndim+1]/rj2)/(ri2 + rj2);

		xi = _Ui[k+1+_Ndim+1];
		xj = _Uj[k+1+_Ndim+1];
		if(ri2>_precision && rj2>_precision)
			_Uroe[2+k+_Ndim] = (xi/ri2 + xj/rj2)/(ri2 + rj2);
		else if(ri2<_precision && rj2>_precision)
			_Uroe[2+k+_Ndim] = xj/_Uj[_Ndim+1];
		else if(ri2>_precision && rj2<_precision)
			_Uroe[2+k+_Ndim] = xi/_Ui[_Ndim+1];
		else
			_Uroe[2+k+_Ndim] = (xi/ri1 + xj/rj1)/(ri1 + rj1);
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Etat de Roe calcule: "<<endl;
		for(int k=0;k<_nVar; k++)//At this point _Uroe[_nVar] is not yet set
			cout<< _Uroe[k]<<endl;
	}
}

void IsothermalTwoFluid::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat de diffusion (m_v, q_v, m_l, q_l)
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _Ui);
	_idm[0] = _nVar*j;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	if(IsBord)
		VecGetValues(_Uextdiff, _nVar, _idm, _Uj);
	else
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

	for(int k=0; k<_nVar; k++)
		_Udiff[k] = (_Ui[k]+_Uj[k])/2;
	double q_2=0;
	for (int i = 0; i<_Ndim;i++)
		q_2+=_Udiff[i+1]*_Udiff[i+1];

	if(_timeScheme==Implicit)
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		double mu1 = _fluides[0]->getViscosity(_Temperature);
		double mu2 = _fluides[1]->getViscosity(_Temperature);
		for(int idim=1;idim<_Ndim+1;idim++)
		{
			_Diffusion[idim*_nVar] =  mu1*_Udiff[idim]/(_Udiff[0]*_Udiff[0]);
			_Diffusion[idim*_nVar+idim] = -mu1/_Udiff[0];
			_Diffusion[(idim+_Ndim+1)*_nVar] =  mu2*_Udiff[idim+_Ndim+1]/(_Udiff[_Ndim+1]*_Udiff[_Ndim+1]);
			_Diffusion[(idim+_Ndim+1)*_nVar+idim+_Ndim+1] = -mu2/_Udiff[_Ndim+1];
		}
	}
}

void IsothermalTwoFluid::convectionMatrices()
{
	//entree: URoe = alpha, p, u1, u2 + ajout dpi pour calcul flux ultérieur
	//sortie: matrices Roe+  et Roe- +Roe si scheme centre

	if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)
		throw CdmathException("Implicitation with primitive variables not yet available for IsothermalTwoFluid model");

	/*Definitions */
	complex< double > tmp;
	double u1_n, u1_2, u2_n, u2_2, u_r2;
	u1_2 = 0; u2_2=0;
	u1_n = 0; u2_n=0;
	// relative velocity
	u_r2=0;
	for(int i=0;i<_Ndim;i++)
	{
		u1_2 += _Uroe[2+i]*_Uroe[2+i];
		u1_n += _Uroe[2+i]*_vec_normal[i];
		u2_2 += _Uroe[2+i+_Ndim]*_Uroe[2+i+_Ndim];
		u2_n += _Uroe[2+i+_Ndim]*_vec_normal[i];
		u_r2 += (_Uroe[2+i]-_Uroe[2+i+_Ndim])*(_Uroe[2+i]-_Uroe[2+i+_Ndim]);
	}
	//Ancienne construction Mat Roe (Dm1,Dm2,Dalp,Dp)
	double alpha = _Uroe[0];
	double p = _Uroe[1];
	double rho1 = _fluides[0]->getDensity(p, _Temperature);
	double rho2 = _fluides[1]->getDensity(p, _Temperature);
	double dpi1 = intPressDef(alpha, u_r2, rho1, rho2);
	double dpi2 = dpi1;
	double invcsq1 = 1/_fluides[0]->vitesseSonTemperature(_Temperature,rho1);
	invcsq1*=invcsq1;
	double invcsq2 = 1/_fluides[1]->vitesseSonTemperature(_Temperature,rho2);
	invcsq2*=invcsq2;
	double g2 = 1/(alpha*rho2*invcsq1+(1-alpha)*rho1*invcsq2);
	double g2press=g2, g2alpha=g2;
	//saving dpi value for flux calculation later
	_Uroe[_nVar]=dpi1 ;

	/***********Calcul des valeurs propres ********/
	Polynoms Poly;
	vector< double > pol_car= Poly.polynome_caracteristique(alpha,  1-alpha, u1_n, u2_n, rho1, rho2,invcsq1,invcsq2, dpi1, dpi2, g2press, g2alpha, g2,_precision);
	for(int ct=0;ct<4;ct++){
		if (abs(pol_car[ct])<_precision*_precision)
			pol_car[ct]=0;
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"pol caract= "<<endl;
		for(int i =0; i<5; i++)
			cout<<pol_car[i]<<"  ";
		cout<< endl;
		cout<<"alpha= "<<alpha<<", p= " << p << ", rho1= " << rho1<< ", rho2= " << rho2<< ", c1= " <<sqrt(1/invcsq1)<<
				", c2= " <<sqrt(1/invcsq2)<<endl;
		cout<< "u1_n= "<<u1_n<<", u2_n= "<<u2_n<< ", dpi1= "<<dpi1<< ", dpi2= "<<dpi2<< endl;
	}
	vector< complex<double> > valeurs_propres = getRacines(pol_car);

	//On ajoute les valeurs propres triviales
	if(_Ndim>1)
	{
		if( !Poly.belongTo(u1_n,valeurs_propres, _precision) )
			valeurs_propres.push_back(u1_n);//vp vapor energy
		if( !Poly.belongTo(u2_n,valeurs_propres, _precision) )
			valeurs_propres.push_back(u2_n);//vp liquid energy
	}
	bool doubleeigenval = norm(valeurs_propres[0] - valeurs_propres[1])<_precision;//norm= suqare of the magnitude
	if(doubleeigenval)
	{
		valeurs_propres[0] = valeurs_propres[valeurs_propres.size()-1];
		valeurs_propres.pop_back();
	}

	int taille_vp = Poly.new_tri_selectif(valeurs_propres,valeurs_propres.size(),_precision);//valeurs_propres.size();//

	_maxvploc=0;
	for(int i =0; i<taille_vp; i++)
		if(fabs(valeurs_propres[i].real())>_maxvploc)
			_maxvploc=fabs(valeurs_propres[i].real());
	if(_maxvploc>_maxvp)
		_maxvp=_maxvploc;

	int existVpCplx = 0,pos_conj;
	double vp_imag_iter;
	for (int ct=0; ct<taille_vp; ct++) {
		vp_imag_iter = valeurs_propres[ct].imag();
		if ( fabs(vp_imag_iter) > _precision ) {
			existVpCplx +=1;
			if ( _part_imag_max < fabs(vp_imag_iter))
				_part_imag_max = fabs(vp_imag_iter);
			//On cherhe le conjugue
			pos_conj = ct+1;
			while(pos_conj<taille_vp && fabs(valeurs_propres[pos_conj].imag()+vp_imag_iter)>_precision)
				pos_conj++;
			if(pos_conj!=ct+1 && pos_conj<taille_vp )
			{
				tmp=valeurs_propres[ct+1];
				valeurs_propres[ct+1]=valeurs_propres[pos_conj];
				valeurs_propres[pos_conj] = tmp;
				ct++;
			}
		}
	}
	if (existVpCplx >0)
		_nbVpCplx +=1;

	//on ordonne les deux premieres valeurs
	/*
	if(valeurs_propres[1].real()<valeurs_propres[0].real())
	{
		tmp=valeurs_propres[0];
		valeurs_propres[0]=valeurs_propres[1];
		valeurs_propres[1]=tmp;
	}
	 */
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" Vp apres tri " << valeurs_propres.size()<<endl;
		for(int ct =0; ct<taille_vp; ct++)
			cout<< "("<<valeurs_propres[ct].real()<< ", " <<valeurs_propres[ct].imag() <<")  ";
		cout<< endl;
	}

	/******** Construction de la matrice de Roe *********/
	//lignes de masse
	for(int i=0; i<_nVar*_nVar;i++)
		_Aroe[i]=0;

	for(int idim=0; idim<_Ndim;idim++)
	{
		_Aroe[1+idim]=_vec_normal[idim];
		_Aroe[1+idim+_Ndim+1]=0;
		_Aroe[(_Ndim+1)*_nVar+1+idim]=0;
		_Aroe[(_Ndim+1)*_nVar+1+idim+_Ndim+1]=_vec_normal[idim];
	}
	//lignes de qdm
	for(int idim=0; idim<_Ndim;idim++)
	{
		//premiere colonne (masse gaz)
		_Aroe[                 (1+idim)*_nVar]=   (alpha *rho2*g2press+dpi1*(1-alpha)*invcsq2*g2alpha)*_vec_normal[idim] - u1_n*_Uroe[2+idim];
		_Aroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar]=((1-alpha)*rho2*g2press-dpi2*(1-alpha)*invcsq2*g2alpha)*_vec_normal[idim];
		//colonnes intermediaires
		for(int jdim=0; jdim<_Ndim;jdim++)
		{
			_Aroe[                 (1+idim)*_nVar + jdim + 1]         = _Uroe[      2+idim]*_vec_normal[jdim];
			_Aroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar + jdim + 1+_Ndim+1] = _Uroe[_Ndim+2+idim]*_vec_normal[jdim];
		}
		//matrice identite
		_Aroe[                 (1+idim)*_nVar +          idim + 1] += u1_n;
		_Aroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar + _Ndim+1+ idim + 1] += u2_n;
		//troisieme colonne (masse liquide)
		_Aroe[                (1+idim)*_nVar + _Ndim+1]=   (alpha *rho1*g2press-dpi1*alpha*invcsq1*g2alpha)*_vec_normal[idim];
		_Aroe[(_Ndim+1)*_nVar+(1+idim)*_nVar + _Ndim+1]=((1-alpha)*rho1*g2press+dpi2*alpha*invcsq1*g2alpha)*_vec_normal[idim] - u2_n*_Uroe[1+idim+_Ndim+1];
	}

	/******* Construction des matrices de decentrement *****/
	if(_spaceScheme == centered){
		if(_entropicCorrection)
			throw CdmathException("IsothermalTwoFluid::convectionMatrices: entropic scheme not available for centered scheme");
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i]=0;
	}
	if( _spaceScheme ==staggered){
		if(_entropicCorrection)//To do: study entropic correction for staggered
			throw CdmathException("IsothermalTwoFluid::convectionMatrices: entropic scheme not yet available for staggered scheme");
		/******** Construction du decentrement du type decale *********/
		//lignes de masse
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i]=0;

		for(int idim=0; idim<_Ndim;idim++)
		{
			_absAroe[1+idim]=_vec_normal[idim];
			_absAroe[1+idim+_Ndim+1]=0;
			_absAroe[(_Ndim+1)*_nVar+1+idim]=0;
			_absAroe[(_Ndim+1)*_nVar+1+idim+_Ndim+1]=_vec_normal[idim];
		}
		//lignes de qdm
		for(int idim=0; idim<_Ndim;idim++)
		{
			//premiere colonne (masse gaz)
			_absAroe[                 (1+idim)*_nVar]=   (-alpha *rho2*g2press+dpi1*(1-alpha)*invcsq2*g2alpha)*_vec_normal[idim] - u1_n*_Uroe[2+idim];
			_absAroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar]=(-(1-alpha)*rho2*g2press-dpi2*(1-alpha)*invcsq2*g2alpha)*_vec_normal[idim];
			//colonnes intermediaires
			for(int jdim=0; jdim<_Ndim;jdim++)
			{
				_absAroe[                 (1+idim)*_nVar + jdim + 1]         = _Uroe[      2+idim]*_vec_normal[jdim];
				_absAroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar + jdim + 1+_Ndim+1] = _Uroe[_Ndim+2+idim]*_vec_normal[jdim];
			}
			//matrice identite
			_absAroe[                 (1+idim)*_nVar +          idim + 1] += u1_n;
			_absAroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar + _Ndim+1+ idim + 1] += u2_n;
			//troisieme colonne (masse liquide)
			_absAroe[                (1+idim)*_nVar + _Ndim+1]=   (-alpha *rho1*g2press-dpi1*alpha*invcsq1*g2alpha)*_vec_normal[idim];
			_absAroe[(_Ndim+1)*_nVar+(1+idim)*_nVar + _Ndim+1]=(-(1-alpha)*rho1*g2press+dpi2*alpha*invcsq1*g2alpha)*_vec_normal[idim] - u2_n*_Uroe[1+idim+_Ndim+1];
		}

		double signu1,signu2;
		if(u1_n>0)
			signu1=1;
		else if (u1_n<0)
			signu1=-1;
		else
			signu1=0;
		if(u2_n>0)
			signu2=1;
		else if (u2_n<0)
			signu2=-1;
		else
			signu2=0;

		for(int i=0; i<(1+_Ndim)*_nVar;i++){
			_absAroe[i] *= signu1;
			_absAroe[i+(1+_Ndim)*_nVar] *= signu2;
		}
	}
	if(_spaceScheme==upwind || _spaceScheme ==lowMach)
	{
		vector< complex< double > > y (taille_vp,0);
		Polynoms Poly;
		for( int i=0 ; i<taille_vp ; i++)
			y[i] = Poly.abs_generalise(valeurs_propres[i]);

		if(_entropicCorrection)
		{
			entropicShift(_vec_normal);
			y[0] +=_entropicShift[0];
			y[taille_vp-1] +=_entropicShift[2];
			for( int i=1 ; i<taille_vp-1 ; i++)
				y[i] +=_entropicShift[1];
		}

		Poly.abs_par_interp_directe(taille_vp,valeurs_propres, _Aroe, _nVar,_precision, _absAroe,y);
		if( _spaceScheme ==pressureCorrection){
			for( int i=0 ; i<_Ndim ; i++)
				for( int j=0 ; j<_Ndim ; j++){
					_absAroe[(1+i)*_nVar+1+j]-=alpha*(valeurs_propres[1].real()-valeurs_propres[0].real())/2*_vec_normal[i]*_vec_normal[j];
					_absAroe[(2+_Ndim+i)*_nVar+2+_Ndim+j]-=(1-alpha)*(valeurs_propres[1].real()-valeurs_propres[0].real())/2*_vec_normal[i]*_vec_normal[j];
				}
		}
		else if( _spaceScheme ==lowMach){
			double M=max(fabs(u1_2),fabs(u2_2))/_maxvploc;
			for( int i=0 ; i<_Ndim ; i++)
				for( int j=0 ; j<_Ndim ; j++){
					_absAroe[(1+i)*_nVar+1+j]-=(1-M)*alpha*(valeurs_propres[1].real()-valeurs_propres[0].real())/2*_vec_normal[i]*_vec_normal[j];
					_absAroe[(2+_Ndim+i)*_nVar+2+_Ndim+j]-=(1-M)*(1-alpha)*(valeurs_propres[1].real()-valeurs_propres[0].real())/2*_vec_normal[i]*_vec_normal[j];
				}
		}

	}
	//Calcul de la matrice signe pour VFFC, VFRoe et décentrement des termes source
	vector< complex< double > > valeurs_propres_dist(taille_vp,0);
	for( int i=0 ; i<taille_vp ; i++)
		valeurs_propres_dist[i] = valeurs_propres[i];

	if(_entropicCorrection || _spaceScheme ==pressureCorrection || _spaceScheme ==lowMach){
		InvMatriceRoe( valeurs_propres_dist);
		Poly.matrixProduct(_absAroe, _nVar, _nVar, _invAroe, _nVar, _nVar, _signAroe);
	}
	else if (_spaceScheme==upwind)//upwind sans entropic
		SigneMatriceRoe( valeurs_propres_dist);
	else if (_spaceScheme==centered)//centre  sans entropic
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
	}
	else  if(_spaceScheme ==staggered )
	{
		double signu1,signu2;
		if(u1_n>0)
			signu1=1;
		else if (u1_n<0)
			signu1=-1;
		else
			signu1=0;
		if(u2_n>0)
			signu2=1;
		else if (u2_n<0)
			signu2=-1;
		else
			signu2=0;
		for(int i=0; i<_nVar*_nVar;i++)
			_signAroe[i] = 0;
		_signAroe[0] = signu1;
		_signAroe[(1+_Ndim)*_nVar +1+_Ndim] = signu2;
		for(int i=2; i<_nVar-1;i++){
			_signAroe[i*_nVar+i] = -signu1;
			_signAroe[(i+1+_Ndim)*_nVar+i+1+_Ndim] = -signu2;
		}
	}
	else
		throw CdmathException("IsothermalTwoFluid::convectionMatrices: well balanced option not treated");

	for(int i=0; i<_nVar*_nVar;i++)
	{
		_AroeMinus[i] = (_Aroe[i]-_absAroe[i])/2;
		_AroePlus[i]  = (_Aroe[i]+_absAroe[i])/2;
	}
	if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)//Implicitation using primitive variables
		for(int i=0; i<_nVar*_nVar;i++)
			_AroeMinusImplicit[i] = (_AroeImplicit[i]-_absAroeImplicit[i])/2;
	else
		for(int i=0; i<_nVar*_nVar;i++)
			_AroeMinusImplicit[i] = _AroeMinus[i];

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<endl<<"Matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++)
		{
			for(int j=0; j<_nVar;j++)
				cout << _Aroe[i*_nVar+j]<< " , ";
			cout<<endl;
		}
		cout<<"Valeur absolue matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++){
			for(int j=0; j<_nVar;j++)
				cout<<_absAroe[i*_nVar+j]<<" , ";
			cout<<endl;
		}
		cout<<endl<<"Matrice _AroeMinus"<<endl;
		for(int i=0; i<_nVar;i++)
		{
			for(int j=0; j<_nVar;j++)
				cout << _AroeMinus[i*_nVar+j]<< " , ";
			cout<<endl;
		}
		cout<<endl<<"Matrice _AroePlus"<<endl;
		for(int i=0; i<_nVar;i++)
		{
			for(int j=0; j<_nVar;j++)
				cout << _AroePlus[i*_nVar+j]<< " , ";
			cout<<endl;
		}
	}
}

void IsothermalTwoFluid::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	//To do controle signe des vitesses pour CL entree/sortie
	int k;
	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne
	double q1_n=0, q2_n=0;//quantité de mouvement normale à la face limite
	for(k=0; k<_Ndim; k++){
		q1_n+=_externalStates[(k+1)]*normale[k];
		q2_n+=_externalStates[(k+1+1+_Ndim)]*normale[k];
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Boundary conditions for group "<< nameOfGroup<< ", inner cell j= "<<j << " face unit normal vector "<<endl;
		for(k=0; k<_Ndim; k++){
			cout<<normale[k]<<", ";
		}
		cout<<endl;
	}
	if (_limitField[nameOfGroup].bcType==Wall){

		//Pour la convection, inversion du sens des vitesses
		for(k=0; k<_Ndim; k++){
			_externalStates[(k+1)]-= 2*q1_n*normale[k];
			_externalStates[(k+1+1+_Ndim)]-= 2*q2_n*normale[k];
		}

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		//Pour la diffusion, paroi à vitesses imposees
		_externalStates[1]=_externalStates[0]*_limitField[nameOfGroup].v_x[0];
		_externalStates[2+_Ndim]=_externalStates[1+_Ndim]*_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
			_externalStates[3+_Ndim]=_externalStates[1+_Ndim]*_limitField[nameOfGroup].v_y[1];
			if(_Ndim==3)
			{
				_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
				_externalStates[4+_Ndim]=_externalStates[1+_Ndim]*_limitField[nameOfGroup].v_z[1];
			}
		}
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
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double alpha=_limitField[nameOfGroup].alpha;
		double pression=_Vj[1];
		double T=_Temperature;
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		_externalStates[0]=alpha*rho_v;
		_externalStates[1]=alpha*rho_v*(_limitField[nameOfGroup].v_x[0]);
		_externalStates[1+_Ndim]=(1-alpha)*rho_l;
		_externalStates[2+_Ndim]=(1-alpha)*rho_l*(_limitField[nameOfGroup].v_x[1]);
		if(_Ndim>1)
		{
			_externalStates[2]=_externalStates[0]*_limitField[nameOfGroup].v_y[0];
			_externalStates[3+_Ndim]=_externalStates[1+_Ndim]*_limitField[nameOfGroup].v_y[1];
			if(_Ndim==3)
			{
				_externalStates[3]=_externalStates[0]*_limitField[nameOfGroup].v_z[0];
				_externalStates[4+_Ndim]=_externalStates[1+_Ndim]*_limitField[nameOfGroup].v_z[1];
			}
		}
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
		hydroPress*=_externalStates[0]+_externalStates[_Ndim];//multiplication by rho the total density

		//Building the external state
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double pression=_limitField[nameOfGroup].p+hydroPress;
		double alpha=_limitField[nameOfGroup].alpha;
		double T=_Temperature;
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		_externalStates[0]=alpha*rho_v;
		_externalStates[1+_Ndim]=(1-alpha)*rho_l;

		for(k=1;k<1+_Ndim;k++){
			_externalStates[k]=_externalStates[0]*_Vj[1+k];
			_externalStates[k+1+_Ndim]=_externalStates[1+_Ndim]*_Vj[k+1+_Ndim];
		}
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Uextdiff);

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"Etat fantôme InletPressure"<<endl;
			for(int k=0;k<_nVar;k++)
				cout<<_externalStates[k]<<endl;
		}
	}
	else if (_limitField[nameOfGroup].bcType==Outlet){
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		//Computation of the hydrostatic contribution : scalar product between gravity vector and position vector
		Cell Cj=_mesh.getCell(j);
		double hydroPress=Cj.x()*_GravityField3d[0];
		if(_Ndim>1){
			hydroPress+=Cj.y()*_GravityField3d[1];
			if(_Ndim>2)
				hydroPress+=Cj.z()*_GravityField3d[2];
		}
		hydroPress*=_externalStates[0]+_externalStates[_Ndim];//multiplication by rho the total density

		//Building the external state
		VecGetValues(_primitiveVars, _nVar, _idm,_Vj);
		double pression_int=_Vj[1];
		double pression_ext=_limitField[nameOfGroup].p+hydroPress;
		double T=_Vj[_nVar-1];
		double rho_v_int=_fluides[0]->getDensity(pression_int,T);
		double rho_l_int=_fluides[1]->getDensity(pression_int,T);
		double rho_v_ext=_fluides[0]->getDensity(pression_ext,T);
		double rho_l_ext=_fluides[1]->getDensity(pression_ext,T);

		for(k=0;k<1+_Ndim;k++){
			_externalStates[k]*=rho_v_ext/rho_v_int;
			_externalStates[k+1+_Ndim]*=rho_l_ext/rho_l_int;
		}
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
	else{
		cout<<"!!!!!!!!!!!!!!!!! Error IsothermalTwoFluid::setBoundaryState !!!!!!!!!!"<<endl;
		cout<<"!!!!!!!!!!!!! Boundary condition not set for boundary named"<<nameOfGroup<< ", _limitField[nameOfGroup].bcType= "<<_limitField[nameOfGroup].bcType<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		throw CdmathException("Unknown boundary condition");
	}
}

void IsothermalTwoFluid::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)

{
	double mu1 = _fluides[0]->getViscosity(_Temperature);
	double mu2 = _fluides[1]->getViscosity(_Temperature);

	if(mu1==0 && mu2==0)
		return;

	//extraction des valeurs
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_primitiveVars, _nVar, _idm, _Vi);
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion: variables primitives maille " << i<<endl;
		for(int q=0; q<_nVar; q++)
		{
			cout << _Vi[q] << endl;
		}
		cout << endl;
	}

	if(!isBord ){
		for(int k=0; k<_nVar; k++)
			_idm[k] = _nVar*j + k;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
	}
	else
	{
		for(int k=0; k<_nVar; k++)
			_idm[k] = k;

		VecGetValues(_Uextdiff, _nVar, _idm, _phi);
		consToPrim(_phi,_Vj);
	}

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion: variables primitives maille " <<j <<endl;
		for(int q=0; q<_nVar; q++)
		{
			cout << _Vj[q] << endl;
		}
		cout << endl;
	}


	double alpha=(_Vj[0]+_Vi[0])/2;
	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	_phi[_Ndim+1]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_Ndim+1; k++)
	{
		_phi[k]         = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*alpha*mu1*(_Vj[k] - _Vi[k]);//attention car primitif=alpha p u1 u2
		_phi[k+_Ndim+1] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*(1-alpha)*mu2*(_Vj[1+k+_Ndim] - _Vi[1+k+_Ndim]);
	}

	_idm[0] = i;
	VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion au 2nd membre pour la maille " << i << ": "<<endl;
		for(int q=0; q<_nVar; q++)
		{
			cout << _phi[q] << endl;
		}
		cout << endl;
	}

	if(!isBord)
	{
		//On change de signe pour l'autre contribution
		for(int k=0; k<_nVar; k++)
			_phi[k] *= -_inv_dxj/_inv_dxi;
		_idm[0] = j;

		VecSetValuesBlocked(_b, 1, _idm, _phi, ADD_VALUES);
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout << "Contribution diffusion au 2nd membre pour la maille  " << j << ": "<<endl;
		for(int q=0; q<_nVar; q++)
		{
			cout << _phi[q] << endl;
		}
		cout << endl;

		if(_timeScheme==Implicit)
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
}

void IsothermalTwoFluid::sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i)
{
	double m1=Ui[0],m2=Ui[1+_Ndim],alpha=Vi[0], P=Vi[1];
	double norm_ur=0, Gamma;
	for(int k=0; k<_Ndim; k++)
		norm_ur+=(Vi[2+k]-Vi[2+k+_Ndim])*(Vi[2+k]-Vi[2+k+_Ndim]);
	norm_ur=sqrt(norm_ur);

	if(i>=0 &&_Temperature>_Tsat && alpha<1-_precision)
		Gamma=_heatPowerField(i)/_latentHeat;
	else
		Gamma=0;
	for(int k=1; k<_Ndim+1; k++)
	{
		//cout<<"Vi[1+"<<k+_Ndim<<"]="<<Vi[1+k+_Ndim]<<endl;
		Si[k] =_gravite[k]*m1 -_dragCoeffs[0]*norm_ur*(Vi[1+k]-Vi[1+k+_Ndim])+ Gamma*Vi[1+k+_Ndim];//interfacial velocity= ul
		Si[k+_Ndim+1] =_gravite[k+_Ndim+1]*m2 + _dragCoeffs[0]*norm_ur*(Vi[1+k]-Vi[1+k+_Ndim])- Gamma*Vi[1+k+_Ndim];
	}
	if(true){//heated boiling
		Si[0]=Gamma;
		Si[1+_Ndim]=-Gamma;
	}
	else if (P<_Psat && alpha<1-_precision){//flash boiling
		Si[0]=-_dHsatl_over_dp*_dp_over_dt(i)/_latentHeat;
		Si[1+_Ndim]=_dHsatl_over_dp*_dp_over_dt(i)/_latentHeat;
	}
	else{
		Si[0]=0;
		Si[1+_Ndim]=0;
	}

	if(_timeScheme==Implicit)
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_GravityImplicitationMatrix[i] = 0;
		for(int i=0; i<_nVar/2;i++)
			_GravityImplicitationMatrix[i*_nVar]=-_gravite[i];
		for(int i=_nVar/2; i<_nVar;i++)
			_GravityImplicitationMatrix[i*_nVar+_nVar/2]=-_gravite[i];
	}

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalTwoFluid::sourceVector"<<endl;
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
void IsothermalTwoFluid::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{
	double norm_u1=0, u1_n=0, norm_u2=0, u2_n=0, m1, m2;
	for(int i=0;i<_Ndim;i++){
		u1_n += _Uroe[2+i]      *_vec_normal[i];
		u2_n += _Uroe[2+i+_Ndim]*_vec_normal[i];
	}
	pressureLoss[0]=0;
	pressureLoss[1+_Ndim]=0;
	if(u1_n>0){
		for(int i=0;i<_Ndim;i++)
			norm_u1 += Vi[1+i]*Vi[1+i];
		norm_u1=sqrt(norm_u1);
		m1=Ui[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-K*m1*norm_u1*Vi[1+i];
	}
	else{
		for(int i=0;i<_Ndim;i++)
			norm_u1 += Vj[1+i]*Vj[1+i];
		norm_u1=sqrt(norm_u1);
		m1=Uj[0];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[1+i]=-K*m1*norm_u1*Vj[1+i];
	}
	if(u2_n>0){
		for(int i=0;i<_Ndim;i++)
			norm_u2 += Vi[2+i+_Ndim]*Vi[2+i+_Ndim];
		norm_u2=sqrt(norm_u2);
		m2=Ui[1+_Ndim];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[2+i+_Ndim]=-K*m2*norm_u2*Vi[2+i+_Ndim];
	}
	else{
		for(int i=0;i<_Ndim;i++)
			norm_u2 += Vj[2+i+_Ndim]*Vj[2+i+_Ndim];
		norm_u2=sqrt(norm_u2);
		m2=Uj[1+_Ndim];
		for(int i=0;i<_Ndim;i++)
			pressureLoss[2+i+_Ndim]=-K*m2*norm_u2*Vj[2+i+_Ndim];
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"IsothermalTwoFluid::pressureLossVector K= "<<K<<endl;
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

void IsothermalTwoFluid::porosityGradientSourceVector()
{
	double u1_ni=0, u1_nj=0, u2_ni=0, u2_nj=0, rho1i, rho2i, rho1j, rho2j, pi=_Vi[1], pj=_Vj[1], pij1, pij2, alphaij=_Uroe[0];
	for(int i=0;i<_Ndim;i++) {
		u1_ni += _Vi[2+i]*_vec_normal[i];
		u2_ni += _Vi[2+_Ndim+i]*_vec_normal[i];
		u1_nj += _Vj[2+i]*_vec_normal[i];
		u2_nj += _Vj[2+_Ndim+i]*_vec_normal[i];
	}
	_porosityGradientSourceVector[0]=0;
	_porosityGradientSourceVector[1+_Ndim]=0;
	rho1i = _fluides[0]->getDensity(pi, _Temperature);
	rho2i = _fluides[1]->getDensity(pi, _Temperature);
	rho1j = _fluides[0]->getDensity(pj, _Temperature);
	rho2j = _fluides[1]->getDensity(pj, _Temperature);
	pij1=(pi+pj)/2+rho1i*rho1j/2/(rho1i+rho1j)*(u1_ni-u1_nj)*(u1_ni-u1_nj);
	pij2=(pi+pj)/2+rho2i*rho2j/2/(rho2i+rho2j)*(u2_ni-u2_nj)*(u2_ni-u2_nj);
	for(int i=0;i<_Ndim;i++){
		_porosityGradientSourceVector[1+i]      =alphaij*pij1*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
		_porosityGradientSourceVector[2+_Ndim+i]=alphaij*pij2*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
	}
}

/* Funtion of equations of states */
double IsothermalTwoFluid::ecartPression(double m1,double m2, double alpha, double e1, double e2){
	if(alpha>_precision*_precision&& alpha<1-_precision*_precision)
		return _fluides[0]->getPressure((m1/alpha)*e1,m1/alpha) - _fluides[1]->getPressure((m2/(1-alpha))*e2,m2/(1-alpha));
	else if(alpha<=_precision*_precision)
	{
		//cout<<"Warning ecartPression, alpha close to 0"<<endl;
		//       if(m1<_precision*_precision*_precision)
		// 	return 0;
		//       else
		return alpha/_precision*_precision;
	}
	else
	{
		//cout<<"Warning ecartPression, alpha close to 1"<<endl;
		//       if(m2<_precision*_precision*_precision)
		// 	return 0;
		//       else
		return -(1-alpha)/_precision*_precision;
	}
}

double IsothermalTwoFluid::ecartPressionDerivee(double m1,double m2, double alpha, double e1, double e2){
	if(alpha>_precision*_precision && alpha<1-_precision*_precision )
		return -(m1/alpha)/alpha*e1*_fluides[0]->getPressureDerivativeRhoE() - (m2/(1-alpha))/(1-alpha)*e2*_fluides[1]->getPressureDerivativeRhoE();
	else if (alpha<_precision*_precision)
		return -1/_precision*_precision;
	else
		return 1/_precision*_precision;
}

double IsothermalTwoFluid::intPressDef(double alpha, double u_r2, double rho1, double rho2)
{
	return  _intPressCoeff*alpha*(1-alpha)*rho1*rho2*u_r2/( alpha*rho2+(1-alpha)*rho1)
			+alpha*(1-alpha)*rho1*rho2*u_r2/((alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1))*u_r2
			*(alpha*alpha*rho2-(1-alpha)*(1-alpha)*rho1)
			*(alpha*alpha*rho2*rho2/(_fluides[0]->vitesseSonTemperature(_Temperature,rho1)*_fluides[0]->vitesseSonTemperature(_Temperature,rho1))
					-(1-alpha)*(1-alpha)*rho1*rho1/(_fluides[1]->vitesseSonTemperature(_Temperature,rho2)*_fluides[1]->vitesseSonTemperature(_Temperature,rho2)));
}

void IsothermalTwoFluid::entropicShift(double* n)
{
	vector< double > pol_car;

	/*Left values */
	double u1_n = 0, u1_2 = 0, u2_n = 0, u2_2 = 0, u_r2 = 0;
	for(int i=0;i<_Ndim;i++)
	{
		u1_2 += _l[2+i]*_l[2+i];
		u1_n += _l[2+i]*n[i];
		u2_2 += _l[1+i+1+_Ndim]*_l[1+i+1+_Ndim];
		u2_n += _l[1+i+1+_Ndim]*n[i];
		u_r2 += (_l[2+i]-_l[1+i+1+_Ndim])*(_l[2+i]-_l[1+i+1+_Ndim]);
	}
	double alpha = _l[0];
	double p = _l[1];
	double rho1 = _fluides[0]->getDensity(p, _Temperature);
	double rho2 = _fluides[1]->getDensity(p, _Temperature);
	double dpi1 = intPressDef(alpha, u_r2, rho1, rho2);
	double dpi2 = dpi1;
	double invcsq1 = 1/_fluides[0]->vitesseSonTemperature(_Temperature,rho1);
	invcsq1*=invcsq1;
	double invcsq2 = 1/_fluides[1]->vitesseSonTemperature(_Temperature,rho2);
	invcsq2*=invcsq2;

	Polynoms Poly;
	pol_car= Poly.polynome_caracteristique(alpha,  1-alpha, u1_n, u2_n, rho1, rho2, invcsq1, invcsq2, dpi1, dpi2);
	for(int ct=0;ct<4;ct++){
		if (abs(pol_car[ct])<_precision*_precision)
			pol_car[ct]=0;
	}
	vector< complex<double> > vp_left = getRacines(pol_car);
	int taille_vp_left = Poly.new_tri_selectif(vp_left,vp_left.size(),_precision);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Entropic shift left eigenvalues: "<<endl;
		for(unsigned int ct =0; ct<vp_left.size(); ct++)
			cout<<"("<< vp_left[ct].real()<< ", " <<vp_left[ct].imag() << ")";
		cout<<endl;
	}

	/*right values */
	u1_2 = 0; u2_2=0;
	u1_n = 0; u2_n=0;
	u_r2=0;
	for(int i=0;i<_Ndim;i++)
	{
		u1_2 += _r[2+i]*_r[2+i];
		u1_n += _r[2+i]*n[i];
		u2_2 += _r[1+i+1+_Ndim]*_r[1+i+1+_Ndim];
		u2_n += _r[1+i+1+_Ndim]*n[i];
		u_r2 += (_r[2+i]-_r[1+i+1+_Ndim])*(_r[2+i]-_r[1+i+1+_Ndim]);
	}
	alpha = _r[0];
	p = _r[1];
	rho1 = _fluides[0]->getDensity(p, _Temperature);
	rho2 = _fluides[1]->getDensity(p, _Temperature);
	dpi1 = intPressDef(alpha, u_r2, rho1, rho2);
	dpi2 = dpi1;
	invcsq1 = 1/_fluides[0]->vitesseSonTemperature(_Temperature,rho1);
	invcsq1*=invcsq1;
	invcsq2 = 1/_fluides[1]->vitesseSonTemperature(_Temperature,rho2);
	invcsq2*=invcsq2;

	pol_car= Poly.polynome_caracteristique(alpha,  1-alpha, u1_n, u2_n, rho1, rho2, invcsq1, invcsq2, dpi1, dpi2);
	for(int ct=0;ct<4;ct++){
		if (abs(pol_car[ct])<_precision*_precision)
			pol_car[ct]=0;
	}
	vector< complex<double> > vp_right = getRacines(pol_car);
	int taille_vp_right = Poly.new_tri_selectif(vp_right,vp_right.size(),_precision);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Entropic shift right eigenvalues: "<<endl;
		for(unsigned int ct =0; ct<vp_right.size(); ct++)
			cout<<"("<<vp_right[ct].real()<< ", " <<vp_right[ct].imag() <<")";
		cout<< endl;
	}
	_entropicShift[0] = abs(vp_left[0]-vp_right[0]);
	_entropicShift[2] = abs(vp_left[taille_vp_left-1]-vp_right[taille_vp_right-1]);
	_entropicShift[1]=0;
	for(int i=1;i<min(taille_vp_right-1,taille_vp_left-1);i++)
		_entropicShift[1] = max(_entropicShift[1],abs(vp_left[i]-vp_right[i]));
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"eigenvalue jumps "<<endl;
		cout<< _entropicShift[0] << ", " << _entropicShift[1] << ", "<< _entropicShift[2] <<endl;
	}
}

void IsothermalTwoFluid::computeScaling(double maxvp)
{
	//	double alphaScaling;
	//	if(_guessalpha>_precision && _guessalpha<1-_precision)
	//		alphaScaling=_guessalpha;
	//	else
	//		alphaScaling=0.5;

	_blockDiag[0]=1;//alphaScaling;
	_invBlockDiag[0]=1;//_blockDiag[0];
	_blockDiag[1+_Ndim]=1;//-alphaScaling;
	_invBlockDiag[1+_Ndim]=1.0;//_blockDiag[1+_Ndim];
	for(int q=1; q<_Ndim+1; q++)
	{
		_blockDiag[q]=1/(maxvp*maxvp);
		_invBlockDiag[q]=1;//_blockDiag[q];
		_blockDiag[q+1+_Ndim]=1/(maxvp*maxvp);
		_invBlockDiag[q+1+_Ndim]=1;//_blockDiag[q+1+_Ndim];
	}
}

void IsothermalTwoFluid::jacobian(const int &j, string nameOfGroup,double * normale)
{
	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_Jcb[k] = 0;//No implicitation at this stage

	// loop on boundaries
	if (_limitField[nameOfGroup].bcType==Wall)
	{
		for(k=0; k<_nVar;k++)
			_Jcb[k*_nVar + k] = 1;
		for(k=1; k<1+_Ndim;k++)
			for(int l=1; l<1+_Ndim;l++){
				_Jcb[k*_nVar + l] -= 2*normale[k-1]*normale[l-1];
				_Jcb[(k+1+_Ndim)*_nVar + l+1+_Ndim] -= 2*normale[k-1]*normale[l-1];			
			}
	}
	//not wall
	else if (_limitField[nameOfGroup].bcType==Inlet)
	{
		_Jcb[0] = 1;
		_Jcb[(1+_Ndim)*_nVar +1+_Ndim] = 1;
		_Jcb[_nVar]=_limitField[nameOfGroup].v_x[0];
		_Jcb[(2+_Ndim)*_nVar +1+_Ndim] =_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_Jcb[2*_nVar]=_limitField[nameOfGroup].v_y[0];
			_Jcb[(3+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_y[1];
			if(_Ndim==3)
			{
				_Jcb[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				_Jcb[(4+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_z[1];

			}
		}
	}
	// not wall, not inlet
	else if (_limitField[nameOfGroup].bcType==Outlet){
		_idm[0] = j*_nVar;
		for(k=1; k<_nVar;k++)
		{_idm[k] = _idm[k-1] + 1;}
		VecGetValues(_conservativeVars, _nVar, _idm, _phi);
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
	}
	else  if (_limitField[nameOfGroup].bcType!=Neumann && _limitField[nameOfGroup].bcType!=InletPressure)// not wall, not inlet, not outlet
	{
		cout << "group named "<<nameOfGroup << " : unknown boundary condition" << endl;
		throw CdmathException("IsothermalTwoFluid::jacobianDiff: This boundary condition is not treated");
	}
}

void IsothermalTwoFluid::jacobianDiff(const int &j, string nameOfGroup)
{

	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;
	if (_limitField[nameOfGroup].bcType==Wall){
		_JcbDiff[0] = 1;
		_JcbDiff[(1+_Ndim)*_nVar +1+_Ndim] = 1;
		_JcbDiff[_nVar]=_limitField[nameOfGroup].v_x[0];
		_JcbDiff[(2+_Ndim)*_nVar +1+_Ndim] =_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_JcbDiff[2*_nVar]=_limitField[nameOfGroup].v_y[0];
			_JcbDiff[(3+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_y[1];

			if(_Ndim==3)
			{
				_JcbDiff[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				_JcbDiff[(4+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_z[1];
			}
		}
	} else if (_limitField[nameOfGroup].bcType==Inlet){
		/*
		_JcbDiff[0] = 1;
		_JcbDiff[(1+_Ndim)*_nVar +1+_Ndim] = 1;
		_JcbDiff[_nVar]=_limitField[nameOfGroup].v_x[0];
		_JcbDiff[(2+_Ndim)*_nVar +1+_Ndim] =_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_JcbDiff[2*_nVar]=_limitField[nameOfGroup].v_y[0];
			_JcbDiff[(3+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_y[1];

			if(_Ndim==3)
			{
				_JcbDiff[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				_JcbDiff[(4+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_z[1];
			}
		}
		 */
	} else if (_limitField[nameOfGroup].bcType==Outlet){
		/*
		//extraction de l etat courant et primitives
		_idm[0] = j*_nVar;
		for(k=1; k<_nVar;k++)
		{_idm[k] = _idm[k-1] + 1;}
		VecGetValues(_conservativeVars, _nVar, _idm, _phi);
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		 */
	}
	else if (_limitField[nameOfGroup].bcType!=Neumann && _limitField[nameOfGroup].bcType!=InletPressure){
		cout<<"Condition  limite non traitee pour le bord "<<nameOfGroup<< endl;
		throw CdmathException("IsothermalTwoFluid::jacobianDiff: Condition  limite non traitee");
	}
}

void IsothermalTwoFluid::primToCons(const double *P, const int &i, double *W, const int &j){
	//P=alpha,p,u1,u2
	//W=m1,q1,m2,q2
	W[j*_nVar] = P[i*_nVar]*_fluides[0]->getDensity(P[i*_nVar+1],_Temperature);
	W[j*_nVar+1+_Ndim] = (1-P[i*_nVar])*_fluides[1]->getDensity(P[i*_nVar+1],_Temperature);
	for(int k=0; k<_Ndim; k++)
	{
		W[j*_nVar+(k+1)] = W[j*_nVar]*P[i*_nVar+(k+2)];
		W[j*_nVar+(k+1)+1+_Ndim] = W[j*_nVar+1+_Ndim]*P[i*_nVar+(k+2)+_Ndim];
	}

}

void IsothermalTwoFluid::consToPrim(const double *Wcons, double* Wprim,double porosity)//To do: treat porosity
{
	//P=alpha,p,u1,u2
	//W=m1,q1,m2,q2
	double m1=Wcons[0];
	double m2=Wcons[_Ndim+1];
	double e1 = _internalEnergy1;
	double e2 = _internalEnergy2;

	_minm1=min(m1,_minm1);
	_minm2=min(m2,_minm2);
	if(m1<-_precision || m2<-_precision)
		_nbMaillesNeg+=1;
	if(fabs(m1)<_precision*_precision)
	{
		Wprim[0]=0;
		Wprim[1]=_fluides[1]->getPressure(m2*e2,m2);

		for(int idim=0; idim<_Ndim; idim++)
		{
			Wprim[2+_Ndim+idim] = Wcons[2+_Ndim+idim]/Wcons[1+_Ndim];
			Wprim[2+idim] = Wprim[2+_Ndim+idim];
		}
	}
	else if(fabs(m2)<_precision*_precision)
	{
		Wprim[0]=1;
		Wprim[1]=_fluides[0]->getPressure(m1*e1,m1);
		for(int idim=0; idim<_Ndim; idim++)
		{
			Wprim[2+idim] = Wcons[1+idim]/Wcons[0];
			Wprim[2+_Ndim+idim] = Wprim[2+idim];
		}
	}
	else
	{
		for(int idim=0; idim<_Ndim; idim++)
		{
			Wprim[2+idim] = Wcons[1+idim]/Wcons[0];
			Wprim[2+_Ndim+idim] = Wcons[2+_Ndim+idim]/Wcons[1+_Ndim];
		}
		double alphanewton, alphainf=0, alphasup=1;
		int iterMax=50, iter=0;

		double   dp=ecartPression( m1, m2, _guessalpha, e1, e2);
		double   dpprim=ecartPressionDerivee( m1, m2, _guessalpha, e1, e2);

		while(fabs(dp)>1e3*_precision && iter<iterMax)
		{

			//  if (dp>0)
			// 	     	alphainf = _guessalpha;
			// 	      else
			// 	     	alphansup = _guessalpha;

			// 	      _guessalpha = (alphainf+alphansup)/2;

			if (dp>0)
				alphainf = _guessalpha;
			else
				alphasup = _guessalpha;

			alphanewton = _guessalpha-dp/dpprim;

			if(alphanewton<=alphainf)
			{
				_guessalpha = (9*alphainf+alphasup)/10;
				//cout<< "dichotomie"<<endl;
			}
			else if(alphanewton>=alphasup)
			{
				_guessalpha = (alphainf+9*alphasup)/10;
				//cout<< "dichotomie"<<endl;
			}
			else
				_guessalpha=alphanewton;


			if(_verbose && _nbTimeStep%_freqSave ==0)
				cout<<"consToPrim diphasique iter= " <<iter<<" dp= " << dp<<" dpprim= " << dpprim<< " _guessalpha= " << _guessalpha<<endl;
			dp=ecartPression( m1, m2, _guessalpha, e1, e2);
			dpprim=ecartPressionDerivee( m1, m2, _guessalpha, e1, e2);

			iter++;
		}
		if(_verbose && _nbTimeStep%_freqSave ==0)
			cout<<endl;

		if(iter>=iterMax)
		{
			_err_press_max = max(_err_press_max,fabs(dp/1.e5));
			// 	      if(_guessalpha>0.5)
			// 		_guessalpha=1;
			// 	      else
			// 		_guessalpha=0;
		}
		if(_guessalpha>0.5)
			Wprim[1]=_fluides[0]->getPressure((m1/_guessalpha)*e1,m1/_guessalpha);
		else
			Wprim[1]=_fluides[1]->getPressure((m2/(1-_guessalpha))*e2,m2/(1-_guessalpha));
		Wprim[0]=_guessalpha;
	}
	if (Wprim[1]<0){
		cout << "pressure = "<< Wprim[1] << " < 0 " << endl;
		cout << "Conservative state = ";
		for(int k=0; k<_nVar; k++){
			cout<<Wcons[k]<<", ";
		}
		cout<<endl;
		*_runLogFile<< "IsothermalTwoFluid::consToPrim: negative pressure = "<< Wprim[1] << " < 0 " << endl;
		throw CdmathException("IsothermalTwoFluid::consToPrim: negative pressure");
	}
}

Vector IsothermalTwoFluid::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
	if(_verbose){
		cout<<"Ucons"<<endl;
		cout<<U<<endl;
		cout<<"Vprim"<<endl;
		cout<<V<<endl;
	}

	double phim1=U(0);//phi alpha1 rho1
	double phim2=U(1+_Ndim);//phi alpha2 rho2
	Vector phiq1(_Ndim),phiq2(_Ndim);//phi alpha1 rho1 u1, phi alpha2 rho2 u2
	for(int i=0;i<_Ndim;i++){
		phiq1(i)=U(1+i);
		phiq2(i)=U(2+_Ndim+i);
	}
	double alpha=V(0);
	double pression=V(1);
	Vector vitesse1(_Ndim),vitesse2(_Ndim);
	for(int i=0;i<_Ndim;i++){
		vitesse1(i)=V(2+i);
		vitesse2(i)=V(2+_Ndim+i);
	}

	double vitesse1n=vitesse1*normale;
	double vitesse2n=vitesse2*normale;

	double alpha_roe = _Uroe[0];//Toumi formula
	// interfacial pressure term (hyperbolic correction)
	double dpi = _Uroe[_nVar];

	Vector F(_nVar);
	F(0)=phim1*vitesse1n;
	F(1+_Ndim)=phim2*vitesse2n;
	for(int i=0;i<_Ndim;i++){
		F(1+i)=phim1*vitesse1n*vitesse1(i)+(alpha_roe*pression+dpi*alpha)*porosity*normale(i);
		F(2+_Ndim+i)=phim2*vitesse2n*vitesse2(i)+((1-alpha_roe)*pression-dpi*alpha)*normale(i)*porosity;
	}

	if(_verbose){
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

Vector IsothermalTwoFluid::staggeredVFFCFlux()
{
	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
		throw CdmathException("IsothermalTwoFluid::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding");
	else//_spaceScheme==staggered
	{
		Vector Fij(_nVar);
		double alpha_roe = _Uroe[0];//Toumi formula
		// interfacial pressure term (hyperbolic correction)
		double dpi = _Uroe[_nVar];

		double u1ijn=0, u2ijn=0, phialphaq1n=0, phialphaq2n=0;
		for(int idim=0;idim<_Ndim;idim++){//URoe = alpha, p, u1, u2, dpi
			u1ijn+=_vec_normal[idim]*_Uroe[2+idim];
			u2ijn+=_vec_normal[idim]*_Uroe[2+_Ndim+idim];
		}
		if(u1ijn>=0)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phialphaq1n+=_vec_normal[idim]*_Ui[1+idim];//phi alpha rho u n
			Fij(0)=phialphaq1n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phialphaq1n*_Vi[2+idim]+(alpha_roe*_Vj[1]*_porosityj+dpi*_Vi[0]*_porosityi)*_vec_normal[idim];
		}
		else
		{
			for(int idim=0;idim<_Ndim;idim++)
				phialphaq2n+=_vec_normal[idim]*_Uj[1+idim];//phi alpha rho u n
			Fij(0)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phialphaq2n*_Vj[2+idim]+(alpha_roe*_Vi[1]*_porosityi+dpi*_Vj[0]*_porosityj)*_vec_normal[idim];
		}

		if(u2ijn>=0)
		{
			for(int idim=0;idim<_Ndim;idim++)
				phialphaq2n+=_vec_normal[idim]*_Ui[2+_Ndim+idim];//phi alpha rho u n
			Fij(1+_Ndim)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+_Ndim+idim)=phialphaq2n*_Vi[2+_Ndim+idim]+((1-alpha_roe)*_Vj[1]*_porosityj-dpi*_Vi[0]*_porosityi)*_vec_normal[idim];
		}
		else
		{
			for(int idim=0;idim<_Ndim;idim++)
				phialphaq2n+=_vec_normal[idim]*_Uj[2+_Ndim+idim];//phi alpha rho u n
			Fij(1+_Ndim)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+_Ndim+idim)=phialphaq2n*_Vj[2+_Ndim+idim]+((1-alpha_roe)*_Vi[1]*_porosityi-dpi*_Vj[0]*_porosityj)*_vec_normal[idim];
		}
		return Fij;
	}
}

void IsothermalTwoFluid::applyVFRoeLowMachCorrections(bool isBord, string groupname)
{
	if(_nonLinearFormulation!=VFRoe)
		throw CdmathException("IsothermalTwoFluid::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation");
	else//_nonLinearFormulation==VFRoe
	{
		if(_spaceScheme==lowMach){
			double u1_2=0, u2_2=0;
			for(int i=0;i<_Ndim;i++){
				u1_2 += _Uroe[2+i]*_Uroe[2+i];
				u2_2 += _Uroe[2+i+_Ndim]*_Uroe[2+i+_Ndim];
			}

			double 	c = _maxvploc;//mixture sound speed
			double M=max(sqrt(u1_2),sqrt(u2_2))/c;//Mach number
			_Vij[1]=M*_Vij[1]+(1-M)*(_Vi[1]+_Vj[1])/2;
			primToCons(_Vij,0,_Uij,0);
		}
		else if(_spaceScheme==pressureCorrection)
		{
			if(_pressureCorrectionOrder>2)
				throw CdmathException("IsothermalTwoFluid::applyVFRoeLowMachCorrections pressure correction order can be only 1 or 2 for Isothermal two-fluid model");

			double norm_uij=0, uij_n=0, ui_n=0, uj_n=0;//mean velocities
			double rho1 =  _fluides[0]->getDensity(_Uroe[1],_Temperature);
			double rho2 =  _fluides[1]->getDensity(_Uroe[1],_Temperature);
			double m1=_Uroe[0]*rho1,  m2=(1-_Uroe[0])*rho2;
			double rhom=m1+m2;
			for(int i=0;i<_Ndim;i++)
			{
				norm_uij += (m1*_Uroe[2+i]+m2*_Uroe[2+i+_Ndim])*(m1*_Uroe[2+i]+m2*_Uroe[2+i+_Ndim]);
				uij_n += (m1*_Uroe[2+i]+m2*_Uroe[2+i+_Ndim])*_vec_normal[i];
				ui_n += _Vi[2+i]*_vec_normal[i];
				uj_n += _Vj[2+i]*_vec_normal[i];
			}
			norm_uij=sqrt(norm_uij)/rhom;
			uij_n/=rhom;
			if(norm_uij>_precision)//avoid division by zero
				_Vij[1]=(_Vi[1]+_Vj[1])/2 + uij_n/norm_uij*(_Vj[1]-_Vi[1])/4 - rhom*norm_uij*(uj_n-ui_n)/4;
			else
				_Vij[1]=(_Vi[1]+_Vj[1])/2                                    - rhom*norm_uij*(uj_n-ui_n)/4;

		}
		else if(_spaceScheme==staggered)
		{
			double qij_n=0;
			double rho1 =  _fluides[0]->getDensity(_Uroe[1],_Uroe[_nVar-1]);
			double rho2 =  _fluides[1]->getDensity(_Uroe[1],_Uroe[_nVar-1]);
			double m1=_Uroe[0]*rho1,  m2=(1-_Uroe[0])*rho2;
			double rhom=m1+m2;
			for(int i=0;i<_Ndim;i++)
				qij_n += (m1*_Uroe[2+i]+m2*_Uroe[2+i+_Ndim])*_vec_normal[i];

			if(qij_n>=0){
				_Vij[0]=_Vi[0];
				_Vij[1]=_Vj[1];
				for(int i=0;i<_Ndim;i++)
				{
					_Vij[2+i]		=_Vi[2+i];
					_Vij[2+i+_Ndim] =_Vi[2+i+_Ndim];
				}
			}
			else{
				_Vij[0]=_Vj[0];
				_Vij[1]=_Vi[1];
				for(int i=0;i<_Ndim;i++)
				{
					_Vij[2+i]		=_Vj[2+i];
					_Vij[2+i+_Ndim] =_Vj[2+i+_Ndim];
				}
			}
			primToCons(_Vij,0,_Uij,0);
		}
	}
}

void IsothermalTwoFluid::testConservation()
{
	double SUM, DELTA, x;
	int I;
	for(int i=0; i<_nVar; i++)
	{
		{
			if(i == 0)
				cout << "Masse totale phase " << 0 <<" (kg): ";
			else if( i == 1+_Ndim)
				cout << "Masse totale phase " << 1 <<" (kg): ";
			else
			{
				if(i < 1+_Ndim)
					cout << "Quantite de mouvement totale phase 0 (kg.m.s^-1): ";
				else
					cout << "Quantite de mouvement totale phase 1 (kg.m.s^-1): ";
			}
		}
		SUM = 0;
		DELTA = 0;
		I = i;
		for(int j=0; j<_Nmailles; j++)
		{
			VecGetValues(_old, 1, &I, &x);//on recupere la valeur du champ
			SUM += x*_mesh.getCell(j).getMeasure();
			VecGetValues(_newtonVariation, 1, &I, &x);//on recupere la variation du champ
			DELTA += x*_mesh.getCell(j).getMeasure();
			I += _nVar;
		}
		if(fabs(SUM)>_precision)
			cout << SUM<< ", variation relative: " << fabs(DELTA /SUM)  << endl;
		else
			cout << " a une somme nulle,  variation absolue: " << fabs(DELTA) << endl;
	}
}

void IsothermalTwoFluid::save(){
	string prim(_path+"/IsothermalTwoFluidPrim_");
	string cons(_path+"/IsothermalTwoFluidCons_");
	prim+=_fileName;
	cons+=_fileName;

	PetscInt Ii;
	for (long i = 0; i < _Nmailles; i++){
		/* j = 0 : void fraction
			   j = 1 : pressure
			   j = 2, 3, 4: velocity phase 1
			   j = 5, 6, 7: velocity phase 2 */
		for (int j = 0; j < _nVar; j++){
			Ii = i*_nVar +j;
			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
		}
	}
	if(_saveConservativeField){
		for (long i = 0; i < _Nmailles; i++){
			for (int j = 0; j < _nVar; j++){
				Ii = i*_nVar +j;
				//				cout<<"i= "<<i<< " j= "<<j<<" UU(i,j)= "<<_UU(i,j)<<endl;
				VecGetValues(_conservativeVars,1,&Ii,&_UU(i,j));
			}
		}
		_UU.setTime(_time,_nbTimeStep);
	}
	_VV.setTime(_time,_nbTimeStep);
	if (_nbTimeStep ==0 || _restartWithNewFileName){
		string prim_suppress ="rm -rf "+prim+"_*";
		string cons_suppress ="rm -rf "+cons+"_*";
		system(prim_suppress.c_str());//Nettoyage des précédents calculs identiques
		system(cons_suppress.c_str());//Nettoyage des précédents calculs identiques
		_VV.setInfoOnComponent(0,"Void_fraction");
		_VV.setInfoOnComponent(1,"Pressure_(Pa)");
		_VV.setInfoOnComponent(2,"Velocity1_x_m/s");
		if (_Ndim>1)
			_VV.setInfoOnComponent(3,"Velocity1_y_m/s");
		if (_Ndim>2)
			_VV.setInfoOnComponent(4,"Velocity1_z_m/s");

		_VV.setInfoOnComponent(2+_Ndim,"Velocity2_x_m/s");
		if (_Ndim>1)
			_VV.setInfoOnComponent(3+_Ndim,"Velocity2_y_m/s");
		if (_Ndim>2)
			_VV.setInfoOnComponent(4+_Ndim,"Velocity2_z_m/s");

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
			_UU.setInfoOnComponent(0,"Partial_density1");// (kg/m^3)
			_UU.setInfoOnComponent(1,"Momentum1_x");// phase1  (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(2,"Momentum1_y");// phase1 (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(3,"Momentum1_z");// phase1 (kg/m^2/s)

			_UU.setInfoOnComponent(1+_Ndim,"Partial_density2");// phase2 (kg/m^3)
			_UU.setInfoOnComponent(2+_Ndim,"Momentum2_x");// phase2 (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(3+_Ndim,"Momentum2_y");// phase2 (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(4+_Ndim,"Momentum2_z");// phase2 (kg/m^2/s)

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
	}
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
			// j = 0 : concentration, j=1 : pressure; j = _nVar - 1: temperature; j = 2,..,_nVar-2: velocity
			for (int j = 0; j < _Ndim; j++)//On récupère les composantes de vitesse
			{
				Ii = i*_nVar +2+j;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse1(i,j));
				Ii=i*_nVar +2+j+_Ndim;
				VecGetValues(_primitiveVars,1,&Ii,&_Vitesse2(i,j));
			}
			for (int j = _Ndim; j < 3; j++){//On met à zero les composantes de vitesse si la dimension est <3
				_Vitesse1(i,j)=0;
				_Vitesse2(i,j)=0;
			}
		}
		_Vitesse1.setTime(_time,_nbTimeStep);
		_Vitesse2.setTime(_time,_nbTimeStep);
		if (_nbTimeStep ==0 || _restartWithNewFileName){		
			_Vitesse1.setInfoOnComponent(0,"Velocity_x_(m/s)");
			_Vitesse1.setInfoOnComponent(1,"Velocity_y_(m/s)");
			_Vitesse1.setInfoOnComponent(2,"Velocity_z_(m/s)");

			_Vitesse2.setInfoOnComponent(0,"Velocity_x_(m/s)");
			_Vitesse2.setInfoOnComponent(1,"Velocity_y_(m/s)");
			_Vitesse2.setInfoOnComponent(2,"Velocity_z_(m/s)");

			switch(_saveFormat)
			{
			case VTK :
				_Vitesse1.writeVTK(prim+"_GasVelocity");
				_Vitesse2.writeVTK(prim+"_LiquidVelocity");
				break;
			case MED :
				_Vitesse1.writeMED(prim+"_GasVelocity");
				_Vitesse2.writeMED(prim+"_LiquidVelocity");
				break;
			case CSV :
				_Vitesse1.writeCSV(prim+"_GasVelocity");
				_Vitesse2.writeCSV(prim+"_LiquidVelocity");
				break;
			}
		}
		else{
			switch(_saveFormat)
			{
			case VTK :
				_Vitesse1.writeVTK(prim+"_GasVelocity",false);
				_Vitesse2.writeVTK(prim+"_LiquidVelocity",false);
				break;
			case MED :
				_Vitesse1.writeMED(prim+"_GasVelocity",false);
				_Vitesse2.writeMED(prim+"_LiquidVelocity",false);
				break;
			case CSV :
				_Vitesse1.writeCSV(prim+"_GasVelocity");
				_Vitesse2.writeCSV(prim+"_LiquidVelocity");
				break;
			}
		}
	}

	if (_restartWithNewFileName)
		_restartWithNewFileName=false;
}

