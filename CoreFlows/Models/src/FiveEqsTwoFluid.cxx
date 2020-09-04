/*
 * FiveEqsTwoFluid.cxx
 *
 *  Created on: Jan 28, 2015
 *      Author: M. Ndjinga
 */
#include "FiveEqsTwoFluid.hxx"
#include <cstdlib>

using namespace std;

extern "C" int dgeev_(char *jobvl, char *jobvr, int *n, double *
		a, int *lda, double *wr, double *wi, double *vl,
		int *ldvl, double *vr, int *ldvr, double *work,
		int *lwork, int *info);


FiveEqsTwoFluid::FiveEqsTwoFluid(pressureEstimate pEstimate, int dim){
	_Ndim=dim;
	_nVar=2*(_Ndim+1)+1;
	_nbPhases = 2;
	_dragCoeffs=vector<double>(2,0);
	_fluides.resize(2);
	//Ne pas utiliser la loi de Stephane Dellacherie mais la stiffened gas standard
	if (pEstimate==around1bar300K)//EOS at 1 bar and 373K
	{
		cout<<"Fluid is water-Gas mixture around saturation point 1 bar and 373 K (100°C)"<<endl;
		*_runLogFile<<"Fluid is water-Gas mixture around saturation point 1 bar and 373 K (100°C)"<<endl;
		_fluides[0] = new StiffenedGas(1.34,1555,373,2.5e6);  //ideal gas law for Gas at pressure 1 bar and temperature 100°C: eref1=2.5e6
		_fluides[1] = new StiffenedGas(958,1e5,373,4.2e5,1543,3769);  //stiffened gas law for water at pressure 1 bar and temperature 100°C: eref2=5e5
		_Tsat=373;//saturation temperature at 1 bar
		_hsatl=4.2e5;//water enthalpy at saturation at 1 bar
		_hsatv=2.5e6;//Gas enthalpy at saturation at 1 bar
	}
	else//EOS at 155 bars and 618K
	{
		cout<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		*_runLogFile<<"Fluid is water-Gas mixture around saturation point 155 bars and 618 K (345°C)"<<endl;
		_fluides[0] = new StiffenedGas(102,1.55e7,618,2.44e6, 433,3633);  //stiffened gas law for Gas at pressure 155 bar and temperature 345°C: eref1=2.4e6
		_fluides[1] = new StiffenedGas(594,1.55e7,618,1.6e6, 621,3100);  //stiffened gas law for water at pressure 155 bar and temperature 345°C: eref2=1.6e6
		_Tsat=618;//saturation temperature at 155 bars
		_hsatl=1.63e6;//water enthalpy at saturation at 155 bars
		_hsatv=2.6e6;//Gas enthalpy at saturation at 155 bars
	}
	_latentHeat=_hsatv-_hsatl;
	_intPressCoeff=1.5;
}

void FiveEqsTwoFluid::initialize()
{
	cout<<"Initialising the five equation two fluid model"<<endl;
	*_runLogFile<<"Initialising the five equation two fluid model"<<endl;

	if(static_cast<StiffenedGas*>(_fluides[0])==NULL || static_cast<StiffenedGas*>(_fluides[1])==NULL)
		throw CdmathException("FiveEqsTwoFluid::initialize: both phase must have stiffened gas EOS");

	_Uroe = new double[_nVar+1];

	_lCon = new PetscScalar[_nVar];//should be deleted in ::terminate
	_rCon = new PetscScalar[_nVar];//should be deleted in ::terminate
	_JacoMat = new PetscScalar[_nVar*_nVar];//should be deleted in ::terminate

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
		_entropicShift=vector<double>(_nVar);

	ProblemFluid::initialize();
}

void FiveEqsTwoFluid::convectionState( const long &i, const long &j, const bool &IsBord){
	//sortie: WRoe en (alpha, p, u1, u2, T, dm1,dm2,dalpha1,dp)
	//entree: _conservativeVars en (rho1, rho1 u1, rho2, rho2 u2)

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
	/*
	if(_Ui[0]<-(_precision) || _Uj[0]<-(_precision) || _Ui[_Ndim+1]<-(_precision) || _Uj[_Ndim+1]<-(_precision))
	{
		cout<<"Warning: masse partielle negative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cout<< "valeurs a gauche: "<<_Ui[0]<<", "<<_Ui[_Ndim+1]<<", valeurs a droite: "<<_Uj[0]<<", "<<_Uj[_Ndim+1]<<endl;
		// 	  throw CdmathException(" Masse partielle negative, arret de calcul");
	}
	else
		_Ui[0]=max(0.,_Ui[0]);
		_Uj[0]=max(0.,_Uj[0]);
		_Ui[_Ndim+1]=max(0.,_Ui[_Ndim+1]);
		_Uj[_Ndim+1]=max(0.,_Uj[_Ndim+1]);
	 */

	PetscScalar ri1, ri2, rj1, rj2, xi, xj;
	//get _l and _r the primitive  states left and right of the interface
	_idm[0] = _nVar*i;
	for(int k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;
	VecGetValues(_primitiveVars, _nVar, _idm, _l);

	if(IsBord)
	{
		//cout<<"_r is border"<<endl;
		//consToPrim(_Uj, _r);
		_idm[0] = 0;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_Vext, _nVar, _idm, _r);
	}
	else
	{
		_idm[0] = _nVar*j;
		for(int k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _r);
	}
	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"_l: "<<endl;
		for(int k=0;k<_nVar; k++)
			cout<< _l[k]<<endl;
		cout<<"_r: "<<endl;
		for(int k=0;k<_nVar; k++)
			cout<< _r[k]<<endl;
	}
	// _Uroe[0] = \tilde{\alpha_v} = 1 - \tilde{\alpha_l} (formula of Toumi)
	if(2-_l[0]-_r[0] > _precision)
		_Uroe[0] = 1- 2*(1-_l[0])*(1-_r[0])/(2-_l[0]-_r[0]);
	else
		_Uroe[0] = (_l[0]+_r[0])/2;

	if(_l[0]+_r[0] > _precision)
		_Uroe[1] = (_l[1]*_l[0]+_r[1]*_r[0])/(_l[0]+_r[0]);
	else
		_Uroe[1] = (_l[1]*(1-_l[0])+_r[1]*(1-_r[0]))/(2-_l[0]-_r[0]);

	ri1 = sqrt(_Ui[0]); ri2 = sqrt(_Ui[_Ndim+1]);
	rj1 = sqrt(_Uj[0]); rj2 = sqrt(_Uj[_Ndim+1]);
	for(int k=0;k<_Ndim;k++)
	{
		xi = _Ui[k+1];
		xj = _Uj[k+1];
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
			_Uroe[1+k+_Ndim+1] = (xi/ri2 + xj/rj2)/(ri2 + rj2);
		else if(ri2<_precision && rj2>_precision)
			_Uroe[1+k+_Ndim+1] = xj/_Uj[_Ndim+1];
		else if(ri2>_precision && rj2<_precision)
			_Uroe[1+k+_Ndim+1] = xi/_Ui[_Ndim+1];
		else
			_Uroe[1+k+_Ndim+1] = (xi/ri1 + xj/rj1)/(ri1 + rj1);
	}
	_Uroe[_nVar-1]=.5*(_l[_nVar-1]+_r[_nVar-1]);

	//Fin du remplissage dans la fonction convectionMatrices

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"Etat de Roe calcule: "<<endl;
		for(int k=0;k<_nVar; k++)
			cout<< _Uroe[k]<<endl;
	}
}

void FiveEqsTwoFluid::diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord){
	//sortie: matrices et etat Diffusion (alpha1 rho1, q1, alpha2 rho2, q2,T)
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
	double q1_2,q2_2=0;
	for (int i = 0; i<_Ndim;i++){
		q1_2+=_Udiff[        i+1]*_Udiff[        i+1];
		q2_2+=_Udiff[1+_Ndim+i+1]*_Udiff[1+_Ndim+i+1];
	}
	consToPrim(_Udiff,_phi);
	_Udiff[_nVar-1]=_phi[_nVar-1];
	double alpha=_phi[0];
	double Tm=_phi[_nVar-1];
	double mu1 = _fluides[0]->getViscosity(Tm);
	double mu2 = _fluides[1]->getViscosity(Tm);
	double lambda = alpha*_fluides[0]->getConductivity(Tm)+(1-alpha)*_fluides[1]->getConductivity(Tm);
	double Cv1= _fluides[0]->constante("Cv");
	double Cv2= _fluides[1]->constante("Cv");

	if(_timeScheme==Implicit)
	{
		for(int i=0; i<_nVar*_nVar;i++)
			_Diffusion[i] = 0;
		for(int idim=1;idim<_Ndim+1;idim++)
		{
			if(alpha>_precision){
				_Diffusion[idim*_nVar] = alpha* mu1*_Udiff[idim]/(_Udiff[0]*_Udiff[0]);
				_Diffusion[idim*_nVar+idim] = -alpha* mu1/_Udiff[0];
			}
			if(1-alpha>_precision){
				_Diffusion[(idim+_Ndim+1)*_nVar] =  (1-alpha)* mu2*_Udiff[idim+_Ndim+1]/(_Udiff[_Ndim+1]*_Udiff[_Ndim+1]);
				_Diffusion[(idim+_Ndim+1)*_nVar+idim+_Ndim+1] = -(1-alpha)* mu2/_Udiff[_Ndim+1];
			}
		}
		/*//Should correct the formula before using
		int i = (_nVar-1)*_nVar;
		_Diffusion[i]=lambda*(Tm/_Udiff[0]-q1_2/(2*Cv1*_Udiff[0]*_Udiff[0]*_Udiff[0]));
		_Diffusion[i+1+_Ndim]=lambda*(Tm/_Udiff[1+_Ndim]-q2_2/(2*Cv2*_Udiff[1+_Ndim]*_Udiff[1+_Ndim]*_Udiff[1+_Ndim]));
		for(int k=1;k<1+_Ndim;k++)
		{
			_Diffusion[i+k]= lambda*_Udiff[k]/(_Udiff[0]*_Udiff[0]*Cv1);
			_Diffusion[i+k+1+_Ndim]= lambda*_Udiff[k+1+_Ndim]/(_Udiff[1+_Ndim]*_Udiff[+1+_Ndim]*Cv2);
		}
		_Diffusion[_nVar*_nVar-1]=-lambda/(_Udiff[0]*Cv1+_Udiff[1+_Ndim]*Cv2);
		 */
	}
}

void FiveEqsTwoFluid::sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i)
{
	double m1=Ui[0],m2=Ui[1+_Ndim],rho=m1+m2, rhoE=Ui[_nVar-1], T=Vi[_nVar-1],alpha=Vi[0], P=Vi[1];
	double norm_ur=0,norm_u1sq=0,norm_u2sq=0, Gamma;

	for(int k=0; k<_Ndim; k++){
		norm_ur+=(Vi[2+k]-Vi[2+k+_Ndim])*(Vi[2+k]-Vi[2+k+_Ndim]);
		norm_u1sq+=Vi[2+k]*Vi[2+k];
		norm_u2sq+=Vi[2+k+_Ndim]*Vi[2+k+_Ndim];
	}
	norm_ur=sqrt(norm_ur);
	double h=(rhoE-0.5*m1*norm_u1sq-0.5*m2*norm_u2sq+P)/rho;
	double u_int[_Ndim];
	for(int k=0; k<_Ndim; k++)
		u_int[k] = 0.5*(Vi[2+k]+Vi[2+k+_Ndim]);
	//		u_int[k] = Vi[0]*Vi[2+k+_Ndim] + (1-Vi[0])*Vi[2+k];
	if(i>=0 && T>_Tsat && alpha<1-_precision)//if(i>=0 && _hsatv>h  && h>_hsatl && alpha<1-_precision)
		Gamma=_heatPowerField(i)/_latentHeat;
	else//boundary cell, no phase change
		Gamma=0;

	for(int k=1; k<_Ndim+1; k++)
	{
		Si[k] =_gravite[k]*m1-_dragCoeffs[0]*norm_ur*(Vi[1+k]-Vi[1+k+_Ndim]) + Gamma*u_int[k-1];//interfacial velocity= ul
		Si[k+_Ndim+1] =_gravite[k+_Ndim+1]*m2+ _dragCoeffs[0]*norm_ur*(Vi[1+k]-Vi[1+k+_Ndim]) - Gamma*u_int[k-1];
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
	if(i>=0)
		Si[_nVar-1]=_heatPowerField(i);
	else//boundary cell, no heating
		Si[_nVar-1]=0;
	for(int k=0; k<_Ndim; k++)
		Si[_nVar-1] +=_GravityField3d[k]*(Ui[1+k]+Ui[2+k+_Ndim])-_dragCoeffs[0]*norm_ur*(Vi[2+k]-Vi[2+k+_Ndim])*(Vi[2+k]-Vi[2+k+_Ndim]);

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
		cout<<"FiveEqsTwoFluid::sourceVector"<<endl;
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

void FiveEqsTwoFluid::pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj)
{
	double norm_u1=0, u1_n=0, norm_u2=0, u2_n=0, m1, m2;
	for(int i=0;i<_Ndim;i++){
		u1_n += _Uroe[1+i]      *_vec_normal[i];
		u2_n += _Uroe[1+i+_Ndim]*_vec_normal[i];
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
	pressureLoss[_nVar-1]=-K*(m1*norm_u1*norm_u1*norm_u1+m2*norm_u2*norm_u2*norm_u2);

	if(_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<"FiveEqsTwoFluid::pressureLossVector K= "<<K<<endl;
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

void FiveEqsTwoFluid::porosityGradientSourceVector()
{
	double u1_ni=0, u1_nj=0, u2_ni=0, u2_nj=0, rho1i, rho2i, rho1j, rho2j, pi=_Vi[1], pj=_Vj[1], Ti=_Vi[_nVar-1], Tj=_Vj[_nVar-1];
	double pij1, pij2, alphaij=_Uroe[0];
	for(int i=0;i<_Ndim;i++) {
		u1_ni += _Vi[2+i]*_vec_normal[i];
		u2_ni += _Vi[2+_Ndim+i]*_vec_normal[i];
		u1_nj += _Vj[2+i]*_vec_normal[i];
		u2_nj += _Vj[2+_Ndim+i]*_vec_normal[i];
	}
	_porosityGradientSourceVector[0]=0;
	_porosityGradientSourceVector[1+_Ndim]=0;
	rho1i = _fluides[0]->getDensity(pi, Ti);
	rho2i = _fluides[1]->getDensity(pi, Ti);
	rho1j = _fluides[0]->getDensity(pj, Tj);
	rho2j = _fluides[1]->getDensity(pj, Tj);
	pij1=(pi+pj)/2+rho1i*rho1j/2/(rho1i+rho1j)*(u1_ni-u1_nj)*(u1_ni-u1_nj);
	pij2=(pi+pj)/2+rho2i*rho2j/2/(rho2i+rho2j)*(u2_ni-u2_nj)*(u2_ni-u2_nj);
	for(int i=0;i<_Ndim;i++){
		_porosityGradientSourceVector[1+i]      =alphaij*pij1*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
		_porosityGradientSourceVector[2+_Ndim+i]=alphaij*pij2*(_porosityi-_porosityj)*2/(1/_inv_dxi+1/_inv_dxj);
	}
	_porosityGradientSourceVector[_nVar-1]=0;
}

double FiveEqsTwoFluid::intPressDef(double alpha, double u_r2, double rho1, double rho2, double Temperature)
{
	return  _intPressCoeff*alpha*(1-alpha)*rho1*rho2*u_r2/( alpha*rho2+(1-alpha)*rho1);
	+alpha*(1-alpha)*rho1*rho2*u_r2/((alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1)*(alpha*rho2+(1-alpha)*rho1))*u_r2
			*(alpha*alpha*rho2-(1-alpha)*(1-alpha)*rho1)
			*(alpha*alpha*rho2*rho2/(_fluides[0]->vitesseSonTemperature(Temperature,rho1)*_fluides[0]->vitesseSonTemperature(Temperature,rho1))
					-(1-alpha)*(1-alpha)*rho1*rho1/(_fluides[1]->vitesseSonTemperature(Temperature,rho2)*_fluides[1]->vitesseSonTemperature(Temperature,rho2)));
}

Vector FiveEqsTwoFluid::convectionFlux(Vector U,Vector V, Vector normale, double porosity){
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
	double Temperature= V(_nVar-1);

	double vitesse1n=vitesse1*normale;
	double vitesse2n=vitesse2*normale;
	double rho1=_fluides[0]->getDensity(pression,Temperature);
	double rho2=_fluides[1]->getDensity(pression,Temperature);
	double e1_int=_fluides[0]->getInternalEnergy(Temperature,rho1);
	double e2_int=_fluides[1]->getInternalEnergy(Temperature,rho2);

	double alpha_roe = _Uroe[0];//Toumi formula
	// interfacial pressure term (hyperbolic correction)
	double dpi = _Uroe[_nVar];

	Vector F(_nVar);
	F(0)=phim1*vitesse1n;
	F(1+_Ndim)=phim2*vitesse2n;
	for(int i=0;i<_Ndim;i++){
		F(1+i)=phim1*vitesse1n*vitesse1(i)+(alpha_roe*pression+dpi*alpha)*porosity*normale(i);
		F(2+_Ndim+i)=phim2*vitesse2n*vitesse2(i)+((1-alpha_roe)*pression+dpi*(1-alpha))*normale(i)*porosity;
	}
	F(_nVar-1)=phim1*(e1_int+0.5*vitesse1*vitesse1+pression/rho1)*vitesse1n+phim2*(e2_int+0.5*vitesse2*vitesse2+pression/rho2)*vitesse2n;

	if(_verbose){
		cout<<"Flux F(U,V)"<<endl;
		cout<<F<<endl;
	}

	return F;
}

void FiveEqsTwoFluid::convectionJacobianMatrix(double *V, double *n)
{
	complex< double > tmp;
	// enter : V(nVar) : primitive variables
	//            n    : normal vector
	double alp = V[0];
	double p = V[1];
	double Tm = V[_nVar-1];
	double rho1 = _fluides[0]->getDensity(p, Tm);
	double rho2 = _fluides[1]->getDensity(p, Tm);
	double ur_2 = 0;
	for (int idim=0; idim<_Ndim; idim++){
		ur_2 += (V[2+idim]-V[2+idim+_Ndim])*(V[2+idim]-V[2+idim+_Ndim]);
	}
	// interfacial pressure term (hyperbolic correction)
	double dpi1 = intPressDef(alp,ur_2, rho1, rho2,Tm);
	double dpi2 = dpi1;

	/********Prepare the parameters to compute the Jacobian Matrix********/
	/**** coefficients a, b, c ****/
	double inv_a1_2,inv_a2_2,b1,c1,a2,b2,c2;
	double e1,e2;
	e1 = _fluides[0]->getInternalEnergy(V[_nVar-1],rho1);// primitive variable _l[_nVar-1]=Tm
	e2 = _fluides[1]->getInternalEnergy(V[_nVar-1],rho2);
	inv_a1_2 = static_cast<StiffenedGas*>(_fluides[0])->getDiffDensPress(e1);
	inv_a2_2 = static_cast<StiffenedGas*>(_fluides[1])->getDiffDensPress(e2);
	//double getJumpDensInternalEnergy(const double p_l,const double p_r,const double e_l,const double e_r);
	b1 = static_cast<StiffenedGas*>(_fluides[0])->getDiffDensInternalEnergy(p,e1);
	b2 = static_cast<StiffenedGas*>(_fluides[1])->getDiffDensInternalEnergy(p,e2);
	//double getJumpInternalEnergyTemperature();
	c1 = static_cast<StiffenedGas*>(_fluides[0])->getDiffInternalEnergyTemperature();
	c2 = static_cast<StiffenedGas*>(_fluides[1])->getDiffInternalEnergyTemperature();
	/**** coefficients eta,  varrho_2 ****/
	double eta[_Ndim], varrho_2;
	// prefix m is arithmetic mean
	double m1,m2, eta_n;
	double u1[_Ndim],u2[_Ndim], alp_u1[_Ndim],alp_u2[_Ndim];
	m1 = alp*rho1;
	m2 = (1-alp)*rho2;
	varrho_2 =1/((alp*rho2)*inv_a1_2+((1-alp)*rho1)*inv_a2_2);
	eta_n = 0.;
	for (int idim=0; idim<_Ndim; idim++){
		u1[idim] = V[idim+2];
		u2[idim] = V[_Ndim+idim+2];
		alp_u1[idim] = alp*u1[idim];
		alp_u2[idim] = (1-alp)*u2[idim];
		eta_n += (alp_u1[idim]*(1-p/rho1*inv_a1_2)+alp_u2[idim]*(1-p/rho2*inv_a2_2))*n[idim];
	}
	double eta_varrho_2n = eta_n*varrho_2;
	/**** compute jump of Delta T, Delta e1, Delta e2 ****/
	double inv_cm = 1/(c1*m1+c2*m2);
	double DeltaT [_nVar], Delta_e1[_nVar], Delta_e2[_nVar];
	// initialize DeltaT
	DeltaT[0] =-e1;
	DeltaT[1+_Ndim] =-e2;
	DeltaT[_nVar-1] = 1 ;
	for (int idim=0; idim<_Ndim; idim++){
		DeltaT[idim+1] = 0.;
		DeltaT[1+_Ndim+idim+1] = 0;
	}
	for (int idim=0; idim<_Ndim; idim++){
		// wrt mass gas
		DeltaT[0] += 0.5*u1[idim]*u1[idim];
		// wrt mass liquid
		DeltaT[_Ndim+1] += 0.5*u2[idim]*u2[idim];
		// wrt momentum gass
		DeltaT[idim+1] += - u1[idim];//*n[idim]
		// wrt momentum liquid
		DeltaT[_Ndim+idim+2] += - u2[idim];//*n[idim]
	}
	// finalize  DeltaT, Delta_e1 and Delta_e2
	for (int i =0; i< _nVar; ++i){
		DeltaT[i] = inv_cm*DeltaT[i];
		Delta_e1[i] = c1*DeltaT[i];
		Delta_e2[i] = c2*DeltaT[i];
	}
	/**** compute differential flux (energy equation) A5 ****/

	double dF5[_nVar];
	// initialize A5
	for (int i=0; i<_nVar; i++){
		dF5[i]=0;
	}
	dF5[0] = eta_varrho_2n*rho2; // mass gas
	dF5[_Ndim+1] = eta_varrho_2n*rho1; // mass liquid
	for (int idim=0; idim<_Ndim; idim++){
		// momentum gas
		dF5[idim+1]= (e1+p/rho1)*n[idim];
		// momentum liquid
		dF5[_Ndim+idim+2]=(e2+p/rho2)*n[idim];
	}
	// assign the value of A5 (last row of the Roe matrix)
	for (int idim=0; idim<_Ndim; idim++){
		for (int jdim=0; jdim<_Ndim;jdim++){
			dF5[0] -= u1[idim]*u1[jdim]*u1[jdim]*n[idim];// -uin * ujn^2
			dF5[_Ndim+1] -= u2[idim]*u2[jdim]*u2[jdim]*n[idim];
			//momentum gas
			dF5[idim+1] += u1[idim]*u1[jdim]*n[jdim]+0.5*(u1[jdim]*u1[jdim])*n[idim];
			//momentum liquid
			dF5[_Ndim+idim+2] += u2[idim]*u2[jdim]*n[jdim]+0.5*(u2[jdim]*u2[jdim])*n[idim];
		}
	}
	// final dF5
	double coef_e1, coef_e2;
	coef_e1 = - eta_varrho_2n*alp*rho2*b1;
	coef_e2 = - eta_varrho_2n*(1-alp)*rho1*b2;
	for (int idim=0; idim<_Ndim; idim++){
		coef_e1 += (alp*rho1 - alp*p*b1/rho1)*u1[idim]*n[idim];
		coef_e2 += ((1-alp)*rho2 - (1-alp)*p*b2/rho2)*u2[idim]*n[idim];
	}
	for (int i =0; i< _nVar; i++){
		dF5[i] += coef_e1*Delta_e1[i] + coef_e2*Delta_e2[i];
	}
	/******** Construction de la matrice J *********/
	//lignes de masse
	for(int i=0; i<_nVar*_nVar;i++)
		_JacoMat[i]=0.;

	for(int idim=0; idim<_Ndim;idim++)
	{
		_JacoMat[1+idim]=n[idim];
		_JacoMat[1+idim+_Ndim+1]=0.;
		_JacoMat[(_Ndim+1)*_nVar+1+idim]=0.;
		_JacoMat[(_Ndim+1)*_nVar+1+idim+_Ndim+1]=n[idim];
	}
	// Roe Matrix new version
	for(int idim=0; idim<_Ndim;idim++)
		for (int jdim=0; jdim<_Ndim;jdim++){
			// momentum gas (neglect delta alpha and delta P)
			_JacoMat[ (1+idim)*_nVar] += -u1[idim]*u1[jdim]*n[jdim];
			_JacoMat[(1+idim)*_nVar+jdim+1] += u1[idim]*n[jdim];
			_JacoMat[(1+idim)*_nVar+idim+1] += u1[jdim]*n[jdim];
			// momentum liquid (neglect delta alpha and delta P)
			_JacoMat[(2+_Ndim+idim)*_nVar+_Ndim+1] += -u2[idim]*u2[jdim]*n[jdim];
			_JacoMat[(2+_Ndim+idim)*_nVar+_Ndim+1+jdim+1] += u2[idim]*n[jdim];
			_JacoMat[(2+_Ndim+idim)*_nVar+_Ndim+1+idim+1] += u2[jdim]*n[jdim];
		}
	// update \Delta alpha
	/*
	 * (alpha *rho2*varrho_2+dpi1*(1-alpha)*inv_a2_2*varrho_2)*n[idim]
	 */
	for (int idim=0; idim<_Ndim; idim++){
		_JacoMat[ (1+idim)*_nVar] += dpi1*varrho_2*(1-alp)*inv_a2_2*n[idim];
		_JacoMat[ (1+idim)*_nVar+_Ndim+1] += -dpi1*varrho_2*alp*inv_a1_2*n[idim];
		_JacoMat[(2+_Ndim+idim)*_nVar] += - dpi2*varrho_2*(1-alp)*inv_a2_2*n[idim];
		_JacoMat[(2+_Ndim+idim)*_nVar+_Ndim+1] += dpi2*varrho_2*alp*inv_a1_2*n[idim];
		for  (int i=0; i<_nVar; i++){
			_JacoMat[ (1+idim)*_nVar+i]+=dpi1*varrho_2*(-alp*(1-alp)*inv_a2_2*b1*Delta_e1[i]+alp*(1-alp)*inv_a1_2*b2*Delta_e2[i])*n[idim];
			_JacoMat[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i]+=-dpi2*varrho_2*(-alp*(1-alp)*inv_a2_2*b1*Delta_e1[i]+alp*(1-alp)*inv_a1_2*b2*Delta_e2[i])*n[idim];
		}
	}
	// update \Delta P
	for (int idim=0; idim<_Ndim; idim++){
		_JacoMat[ (1+idim)*_nVar] += alp*varrho_2*rho2*n[idim];
		_JacoMat[ (1+idim)*_nVar+_Ndim+1] +=alp* varrho_2*rho1*n[idim];
		_JacoMat[(2+_Ndim+idim)*_nVar] +=  (1-alp)*varrho_2*rho2*n[idim];
		_JacoMat[(2+_Ndim+idim)*_nVar+_Ndim+1] += (1-alp)* varrho_2*rho1*n[idim];
		for  (int i=0; i<_nVar; i++){
			_JacoMat[ (1+idim)*_nVar+i]+=alp*varrho_2*(-alp*rho2*b1*Delta_e1[i] -(1-alp)*rho1*b2*Delta_e2[i])*n[idim];
			_JacoMat[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i]+=(1-alp)*varrho_2*(-alp*rho2*b1*Delta_e1[i] -(1-alp)*rho1*b2*Delta_e2[i])*n[idim];
		}
	}
	// last row (total energy)
	for (int i=0; i<_nVar; i++){
		_JacoMat[(2*_Ndim+2)*_nVar +i] = dF5[i];
	}
}

void FiveEqsTwoFluid::convectionMatrices()
{
	//entree: URoe = alpha, p, u1, u2, T + ajout dpi
	//sortie: matrices Roe+  et Roe- +Roe si schéma centre

	if(_timeScheme==Implicit && _usePrimitiveVarsInNewton)
		throw CdmathException("Implicitation with primitive variables not yet available for FiveEqsTwoFluid model");

	/*Definitions */
	complex< double > tmp;
	double  u_r2 = 0, u1_n=0, u2_n=0;
	// u1_n = u1.n;u2_n = u2.n (scalar product)
	for(int i=0;i<_Ndim;i++)
	{
		//u_r2 += (_Uroe[2+i]-_Uroe[1+i+1+_Ndim])*(_Uroe[2+i]-_Uroe[1+i+1+_Ndim]);
		u_r2 += 0.5*((_l[2+i]-_l[2+i+_Ndim])*(_l[2+i]-_l[2+i+_Ndim])+(_r[2+i]-_r[2+i+_Ndim])*(_r[2+i]-_r[2+i+_Ndim])); //Kieu
		u1_n += _Uroe[2+i]*_vec_normal[i];
		u2_n += _Uroe[1+i+1+_Ndim]*_vec_normal[i];
	}

	//double alpha = (_l[0]+_r[0])*0.5;//Kieu formula
	double alpha = _Uroe[0];//Toumi formula
	//double p = _Uroe[1];//Toumi pressure

	/***********Calcul des valeurs propres ********/

	// ********Prepare the parameters to compute the Roe Matrix******** //
	// **** coefficients eta,  varrho_2 **** //
	double eta[_Ndim], varrho_2;
	double rho1_l =  _fluides[0]->getDensity(_l[1], _l[_Ndim*2+2]);//(p,T)_Ndim*2+2
	double rho2_l =  _fluides[1]->getDensity(_l[1], _l[_Ndim*2+2]);
	double rho1_r =  _fluides[0]->getDensity(_r[1], _r[_Ndim*2+2]);
	double rho2_r =  _fluides[1]->getDensity(_r[1], _r[_Ndim*2+2]);
	// **** coefficients a, b, c **** //
	double inv_a1_2,inv_a2_2,b1,c1,a2,b2,c2;
	double e1_l,e1_r,e2_l,e2_r;
	e1_l = _fluides[0]->getInternalEnergy(_l[_nVar-1],rho1_l);// primitive variable _l[_nVar-1]=Tm
	e2_l = _fluides[1]->getInternalEnergy(_l[_nVar-1],rho2_l);
	e1_r = _fluides[0]->getInternalEnergy(_r[_nVar-1],rho1_r);
	e2_r = _fluides[1]->getInternalEnergy(_r[_nVar-1],rho2_r);
	inv_a1_2 = static_cast<StiffenedGas*>(_fluides[0])->getJumpDensPress(e1_l,e1_r);
	inv_a2_2 = static_cast<StiffenedGas*>(_fluides[1])->getJumpDensPress(e2_l,e2_r);
	//double getJumpDensInternalEnergy(const double p_l,const double p_r,const double e_l,const double e_r);
	b1 = static_cast<StiffenedGas*>(_fluides[0])->getJumpDensInternalEnergy(_l[1],_r[1],e1_l,e1_r);
	b2 = static_cast<StiffenedGas*>(_fluides[1])->getJumpDensInternalEnergy(_l[1],_r[1],e2_l,e2_r);
	//double getJumpInternalEnergyTemperature();
	c1 = static_cast<StiffenedGas*>(_fluides[0])->getJumpInternalEnergyTemperature();
	c2 = static_cast<StiffenedGas*>(_fluides[1])->getJumpInternalEnergyTemperature();

	// prefix m is arithmetic mean
	double m_alp1,m_rho1,m_rho2,m_P,m_e1,m_e2,m_m1,m_m2, eta_n;
	double m_u1[_Ndim],m_u2[_Ndim], m_alp_u1[_Ndim],m_alp_u2[_Ndim];
	double u1_l[_Ndim], u2_l[_Ndim],u1_r[_Ndim],u2_r[_Ndim];
	double u1l_2 = 0, u1r_2 = 0, u2l_2 = 0, u2r_2 = 0;
	m_alp1 =  0.5*(_l[0]+_r[0]);
	m_rho1 = 0.5*(rho1_l+rho1_r);
	m_rho2 = 0.5*(rho2_l+rho2_r);
	m_P = 0.5*(_l[1]+_r[1]);
	m_e1 = 0.5*(e1_l+e1_r);
	m_e2 = 0.5*(e2_l+e2_r);
	m_m1 = 0.5*(_l[0]*rho1_l+_r[0]*rho1_r);
	m_m2 = 0.5*((1-_l[0])*rho2_l+(1-_r[0])*rho2_r);
	varrho_2 =1/((m_alp1*m_rho2)*inv_a1_2+((1-m_alp1)*m_rho1)*inv_a2_2);
	eta_n = 0.;

	for (int idim=0; idim<_Ndim; idim++){
		u1_l[idim] = _l[idim+2];
		u1_r[idim] = _r[idim+2];
		u2_l[idim] = _l[_Ndim+idim+2];
		u2_r[idim] = _r[_Ndim+idim+2];
		m_u1[idim] = 0.5*(u1_l[idim] + u1_r[idim]);
		m_u2[idim] = 0.5*(u2_l[idim] + u2_r[idim]);
		m_alp_u1[idim] = 0.5*(_l[0]*u1_l[idim]+_r[0]*u1_r[idim]);
		m_alp_u2[idim] = 0.5*((1-_l[0])*u2_l[idim]+(1-_r[0])*u2_r[idim]);
		eta_n += (m_alp_u1[idim]*(1-m_P/m_rho1*inv_a1_2)+m_alp_u2[idim]*(1-m_P/m_rho2*inv_a2_2))*_vec_normal[idim];
	}

	double eta_varrho_2n = eta_n*varrho_2;
	//  **** compute jump of Delta T, Delta e1, Delta e2 **** //
	for  (int idim=0; idim<_Ndim; idim++){
		u1l_2 += u1_l[idim]*u1_l[idim];
		u1r_2 += u1_r[idim]*u1_r[idim];
		u2l_2 += u2_l[idim]*u2_l[idim];
		u2r_2 += u2_r[idim]*u2_r[idim];
	}
	double inv_m_cm = 1/(c1*m_m1+c2*m_m2);
	double DeltaT [_nVar], Delta_e1[_nVar], Delta_e2[_nVar];
	// initialize DeltaT
	for (int i=0; i<_nVar; i++){
		DeltaT[i] = 0.;
	}
	DeltaT[0] += -m_e1;
	DeltaT[1+_Ndim] += -m_e2;
	DeltaT[_nVar-1] += 1 ;
	for (int idim=0; idim<_Ndim; idim++){
		// wrt mass gas
		DeltaT[0] += 0.5*_l[idim+2]	*_r[idim+2];//0.5*\tilde{u_g}^2
		// wrt mass liquid
		DeltaT[_Ndim+1] += 0.5*_l[_Ndim+idim+2]	*_r[_Ndim+idim+2];
		// wrt momentum gas
		DeltaT[idim+1] += - m_u1[idim];
		// wrt momentum liquid
		DeltaT[_Ndim+idim+2] += - m_u2[idim];
	}

	// finalize  DeltaT, Delta_e1 and Delta_e2
	for (int i =0; i< _nVar; ++i){
		DeltaT[i] = inv_m_cm*DeltaT[i];
		Delta_e1[i] = c1*DeltaT[i];
		Delta_e2[i] = c2*DeltaT[i];
	}

	// *** compute jump flux (energy equation) A5 *** //

	double A5[_nVar];
	// initialize A5
	for (int i=0; i<_nVar; i++){
		A5[i]=0;
	}
	A5[0] = eta_varrho_2n*m_rho2; // mass gas
	A5[_Ndim+1] = eta_varrho_2n*m_rho1; // mass liquid
	for (int idim=0; idim<_Ndim; idim++){
		// momentum gas
		A5[idim+1] = (m_e1+m_P/m_rho1)*_vec_normal[idim];
		// momentum liquid
		A5[_Ndim+idim+2] = (m_e2+m_P/m_rho2)*_vec_normal[idim];
	}
	// assign the value of A5 (last row of the Roe matrix)
	for (int idim=0; idim<_Ndim; idim++){
		for (int jdim=0; jdim<_Ndim;jdim++){
			A5[0] += 0.5*(0.5*(_l[idim+2]*_l[jdim+2]*_l[jdim+2]+_r[idim+2]*_r[jdim+2]*_r[jdim+2])*_vec_normal[idim]);// m_(uin*uj^2)
			A5[0] -= m_u1[idim]*m_u1[jdim]*m_u1[jdim]*_vec_normal[idim];// -m_uin * m_ujn^2
			A5[0] -= 0.5*(m_u1[idim]*_vec_normal[idim]*0.5*(_l[jdim+2]*_l[jdim+2]+_r[jdim+2]*_r[jdim+2]));//-0.5*m_uin*m_uj^2
			A5[_Ndim+1]+= 0.5*(0.5*(_l[_Ndim+idim+2]*_l[_Ndim+jdim+2]*_l[_Ndim+jdim+2]+_r[_Ndim+idim+2]*_r[_Ndim+jdim+2]*_r[_Ndim+jdim+2])*_vec_normal[idim]);// m_(uin*uj^2)
			A5[_Ndim+1] -= m_u2[idim]*m_u2[jdim]*m_u2[jdim]*_vec_normal[idim];// -m_uin * m_ujn^2
			A5[_Ndim+1] -= 0.5*(m_u2[idim]*0.5*(_l[_Ndim+jdim+2]*_l[_Ndim+jdim+2]+_r[_Ndim+jdim+2]*_r[_Ndim+jdim+2]))*_vec_normal[idim];//-0.5*m_uin*m_uj^2
			//momentum gas
			A5[idim+1] += m_u1[jdim]*_vec_normal[jdim]*m_u1[idim]+0.5*_vec_normal[idim]*0.5*(_l[jdim+2]*_l[jdim+2]+_r[jdim+2]*_r[jdim+2]);
			//momentum liquid
			A5[_Ndim+idim+2] += m_u2[jdim]*_vec_normal[jdim]*m_u2[idim]+0.5*_vec_normal[idim]*0.5*(_l[_Ndim+jdim+2]*_l[_Ndim+jdim+2]+_r[_Ndim+jdim+2]*_r[_Ndim+jdim+2]);
		}
	}

	// final A5
	double coef_e1, coef_e2;
	coef_e1 = - eta_varrho_2n*m_alp1*m_rho2*b1;
	coef_e2 = - eta_varrho_2n*(1-m_alp1)*m_rho1*b2;
	for (int idim=0; idim<_Ndim; idim++){
		coef_e1 += (0.5*(_l[0]*rho1_l*_l[idim+2]+_r[0]*rho1_r*_r[idim+2]) - m_alp_u1[idim]*m_P*b1/m_rho1)*_vec_normal[idim];
		coef_e2 += (0.5*((1-_l[0])*rho2_l*_l[_Ndim+idim+2]+(1-_r[0])*rho2_r*_r[_Ndim+idim+2])-m_alp_u2[idim]*m_P*b2/m_rho2)*_vec_normal[idim];
	}
	for (int i =0; i< _nVar; i++){
		A5[i] += coef_e1*Delta_e1[i] + coef_e2*Delta_e2[i];
	}
	// ******* Construction de la matrice de Roe ******** //
	// interfacial pressure correction
	double T=_Uroe[_nVar-1];
	double dpi1 = intPressDef(alpha, u_r2, m_rho1,m_rho2,T);
	double dpi2 = dpi1;
	//saving dpi value for flux calculation later
	_Uroe[_nVar]=dpi1 ;

	//lignes de masse
	for(int i=0; i<_nVar*_nVar;i++)
		_Aroe[i]=0;
	//	alpha = 0.; dpi1 = 0.; dpi2 = 0.;
	for(int idim=0; idim<_Ndim;idim++)
	{
		_Aroe[1+idim]=_vec_normal[idim];
		_Aroe[1+idim+_Ndim+1]=0;
		_Aroe[(_Ndim+1)*_nVar+1+idim]=0;
		_Aroe[(_Ndim+1)*_nVar+1+idim+_Ndim+1]=_vec_normal[idim];
	}
	// Take into account the convection term in the momentum eqts
	for(int idim=0; idim<_Ndim;idim++)
		for (int jdim=0; jdim<_Ndim;jdim++){
			// momentum gas (neglect delta alpha and delta P)
			_Aroe[ (1+idim)*_nVar] += (0.5*(_l[2+idim]*_l[2+jdim]+_r[2+idim]*_r[2+jdim])-2*m_u1[idim]*m_u1[jdim])*_vec_normal[jdim];
			_Aroe[(1+idim)*_nVar+jdim+1] += m_u1[idim]*_vec_normal[jdim];
			_Aroe[(1+idim)*_nVar+idim+1] += m_u1[jdim]*_vec_normal[jdim];
			// momentum liquid (neglect delta alpha and delta P)
			_Aroe[(2+_Ndim+idim)*_nVar+_Ndim+1] += (0.5*(_l[_Ndim+2+idim]*_l[_Ndim+2+jdim]+_r[_Ndim+2+idim]*_r[_Ndim+2+jdim])-2*m_u2[idim]*m_u2[jdim])*_vec_normal[jdim];
			_Aroe[(2+_Ndim+idim)*_nVar+_Ndim+1+jdim+1] += m_u2[idim]*_vec_normal[jdim];
			_Aroe[(2+_Ndim+idim)*_nVar+_Ndim+1+idim+1] += m_u2[jdim]*_vec_normal[jdim];
		}
	// update \Delta alpha
	for (int idim=0; idim<_Ndim; idim++){
		_Aroe[ (1+idim)*_nVar] += dpi1*varrho_2*(1-m_alp1)*inv_a2_2*_vec_normal[idim];
		_Aroe[ (1+idim)*_nVar+_Ndim+1] += -dpi1*varrho_2*(m_alp1)*inv_a1_2*_vec_normal[idim];
		_Aroe[(2+_Ndim+idim)*_nVar] += - dpi2*varrho_2*(1-m_alp1)*inv_a2_2*_vec_normal[idim];
		_Aroe[(2+_Ndim+idim)*_nVar+_Ndim+1] += dpi2*varrho_2*(m_alp1)*inv_a1_2*_vec_normal[idim];
		for  (int i=0; i<_nVar; i++){
			_Aroe[ (1+idim)*_nVar+i]+=dpi1*varrho_2*(-m_alp1*(1-m_alp1)*inv_a2_2*b1*Delta_e1[i]+m_alp1*(1-m_alp1)*inv_a1_2*b2*Delta_e2[i])*_vec_normal[idim];
			_Aroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i]+=-dpi2*varrho_2*(-m_alp1*(1-m_alp1)*inv_a2_2*b1*Delta_e1[i]+m_alp1*(1-m_alp1)*inv_a1_2*b2*Delta_e2[i])*_vec_normal[idim];
		}
	}
	// update \Delta P
	for (int idim=0; idim<_Ndim; idim++){
		_Aroe[ (1+idim)*_nVar] += alpha*varrho_2*m_rho2*_vec_normal[idim];
		_Aroe[ (1+idim)*_nVar+_Ndim+1] += alpha* varrho_2*m_rho1*_vec_normal[idim];
		_Aroe[(2+_Ndim+idim)*_nVar] += (1-alpha)*varrho_2*m_rho2*_vec_normal[idim];
		_Aroe[(2+_Ndim+idim)*_nVar+_Ndim+1] += (1-alpha)* varrho_2*m_rho1*_vec_normal[idim];
		for  (int i=0; i<_nVar; i++){
			_Aroe[ (1+idim)*_nVar+i] += alpha*varrho_2*(-m_alp1*m_rho2*b1*Delta_e1[i] -(1-m_alp1)*m_rho1*b2*Delta_e2[i])*_vec_normal[idim];
			_Aroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i] += (1-alpha)*varrho_2*(-m_alp1*m_rho2*b1*Delta_e1[i] -(1-m_alp1)*m_rho1*b2*Delta_e2[i])*_vec_normal[idim];
		}
	}
	// last row (total energy)
	for (int i=0; i<_nVar; i++){
		_Aroe[(2*_Ndim+2)*_nVar +i] += A5[i];
	}

	int LDA = _nVar;
	int LDVL = _nVar;
	int LDVR = _nVar;
	int LWORK = 50*_nVar;
	char jobvl[]="N", jobvr[]="N";
	double WORK[LWORK], Aroe[_nVar*_nVar],egvaReal[_nVar],egvaImag[_nVar],
	egVectorL[_nVar*_nVar],egVectorR[_nVar*_nVar];
	int info=0;
	double maxvploc=0;
	std::vector< std::complex<double> > valeurs_propres_dist;
	int taille_vp ;

	for (int i=0; i<_nVar*_nVar; i++)
		Aroe[i] = _Aroe[i];

	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<endl<<"Matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++)
		{
			for(int j=0; j<_nVar;j++)
				cout << _Aroe[i*_nVar+j]<< " , ";
			cout<<endl;
		}
	}
	/******** Compute the eigenvalues and eigenvectors of Roe Matrix (using lapack)*********/
	Polynoms Poly;
	dgeev_(jobvl, jobvr,  &_nVar,
			Aroe,&LDA,egvaReal,egvaImag, egVectorL,
			&LDVL,egVectorR,
			&LDVR, WORK,&LWORK,
			&info);
	if(info <0){
		cout<<"FiveEqsTwoFluid::convectionMatrices: error dgeev_ : argument "<<-info<<" invalid"<<endl;
		*_runLogFile<<"FiveEqsTwoFluid::convectionMatrices: error dgeev_ : argument "<<-info<<" invalid"<<endl;
		throw CdmathException("FiveEqsTwoFluid::convectionMatrices: dgeev_ unsuccessful computation of the eigenvalues ");
	}
	else if(info <0)
	{
		cout<<"Warning FiveEqsTwoFluid::convectionMatrices: dgeev_ did not compute all the eigenvalues, trying Rusanov scheme "<<endl;
		cout<<"Converged eigenvalues are ";
		for(int i=info; i<_nVar; i++)
			cout<<"("<< egvaReal[i]<<","<< egvaImag[i]<<"), ";
		cout<<endl;

		_maxvploc=0;
		for(int i =info; i<_nVar; i++)
			if (fabs(egvaReal[i])>_maxvploc)
				_maxvploc=fabs(egvaReal[i]);
		if(_maxvploc>_maxvp)
			_maxvp=_maxvploc;

		valeurs_propres_dist=std::vector< std::complex<double> > (1,maxvploc);
		taille_vp =1;
	}
	else{
		if (_verbose && _nbTimeStep%_freqSave ==0)
		{
			for(int i=0; i<_nVar; i++)
				cout<<" Vp real part " << egvaReal[i]<<", Imaginary part " << egvaImag[i]<<endl;
		}

		std::vector< std::complex<double> > valeurs_propres(_nVar);

		bool complexEigenvalues = false;
		for(int j=0; j<_nVar; j++){
			if (max(_l[0],_r[0])<_precision && abs(egvaImag[j])<_precision )// Kieu test Monophase
				egvaImag[j] = 0.;
			else
				if (egvaImag[j] >_precision){// Kieu
					complexEigenvalues = true;
				}
			if (abs(_l[0]-_r[0])<_precision*_precision && fabs(egvaImag[j])<_precision)// Kieu interfacial pressure
				egvaImag[j] = 0.;
			valeurs_propres[j] = complex<double>(egvaReal[j],egvaImag[j]);
		}

		taille_vp =Poly.new_tri_selectif(valeurs_propres,valeurs_propres.size(),_precision);

		valeurs_propres_dist=vector< complex< double > >(taille_vp);
		for( int i=0 ; i<taille_vp ; i++)
			valeurs_propres_dist[i] = valeurs_propres[i];
		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<" Vp apres tri " << valeurs_propres_dist.size()<<endl;
			for(int ct =0; ct<taille_vp; ct++)
				cout<< "("<<valeurs_propres_dist[ct].real()<< ", " <<valeurs_propres_dist[ct].imag() <<")  ";
			cout<< endl;
		}

		for(int i =0; i<taille_vp; i++){
			if (fabs(valeurs_propres_dist[i].real())>maxvploc)
				maxvploc=fabs(valeurs_propres_dist[i].real());
			if(maxvploc>_maxvp)
				_maxvp=maxvploc;
		}

		int existVpCplx = 0,pos_conj;
		double vp_imag_iter;
		for (int ct=0; ct<taille_vp; ct++) {
			vp_imag_iter = valeurs_propres_dist[ct].imag();
			if ( fabs(vp_imag_iter) > 100*_precision ) {
				existVpCplx +=1;
				if ( _part_imag_max < fabs(vp_imag_iter))
					_part_imag_max = fabs(vp_imag_iter);
				//On cherhe le conjugue
				pos_conj = ct+1;
				while(pos_conj<taille_vp && fabs(valeurs_propres_dist[pos_conj].imag()+vp_imag_iter)>_precision)
					pos_conj++;
				if(pos_conj!=ct+1 && pos_conj<taille_vp )
				{
					tmp=valeurs_propres_dist[ct+1];
					valeurs_propres_dist[ct+1]=valeurs_propres_dist[pos_conj];
					valeurs_propres_dist[pos_conj] = tmp;
					ct++;
				}
			}
		}
		if (existVpCplx >0)
			_nbVpCplx +=1;
	}
	/******* Construction des matrices de decentrement *******/
	if(_spaceScheme == centered )
	{
		if(_entropicCorrection)
			throw CdmathException("FiveEqsTwoFluid::convectionMatrices: entropic scheme not available for centered scheme");
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i]=0;
		//  if(alpha<_precision || alpha>1-_precision)//rusanov
		// 	    for(int i=0; i<_nVar;i++)
		// 	      _absAroe[i*_nVar+i]=maxvploc;
	}
	if( _spaceScheme ==staggered){//To do: study entropic correction for staggered
		if(_entropicCorrection)//To do: study entropic correction for staggered
			throw CdmathException("FiveEqsTwoFluid::convectionMatrices: entropic scheme not yet available for staggered scheme");
		/******** Construction de la matrice de decentrement staggered *********/
		/***** Compute eta_n **************/
		eta_n = 0.;
		for (int idim=0; idim<_Ndim; idim++){
			eta_n += (m_alp_u1[idim]*(-1-m_P/m_rho1*inv_a1_2)+m_alp_u2[idim]*(-1-m_P/m_rho2*inv_a2_2))*_vec_normal[idim];
		}
		double eta_varrho_2n = eta_n*varrho_2;
		/**** compute jump flux (energy equation) A5 ****/

		double A5[_nVar];
		// initialize A5
		for (int i=0; i<_nVar; i++){
			A5[i]=0;
		}
		A5[0] = eta_varrho_2n*m_rho2; // mass gas
		A5[_Ndim+1] = eta_varrho_2n*m_rho1; // mass liquid
		// assign the value of A5 (last row of the Roe matrix)
		for (int idim=0; idim<_Ndim; idim++){
			for (int jdim=0; jdim<_Ndim;jdim++){
				A5[0] += 0.5*(0.5*(_l[idim+2]*_l[jdim+2]*_l[jdim+2]+_r[idim+2]*_r[jdim+2]*_r[jdim+2])*_vec_normal[idim]);// m_(uin*uj^2)
				A5[0] -= m_u1[idim]*m_u1[jdim]*m_u1[jdim]*_vec_normal[idim];// -m_uin * m_ujn^2
				A5[0] -= 0.5*(m_u1[idim]*_vec_normal[idim]*0.5*(_l[jdim+2]*_l[jdim+2]+_r[jdim+2]*_r[jdim+2]));//-0.5*m_uin*m_uj^2
				A5[_Ndim+1]+= 0.5*(0.5*(_l[_Ndim+idim+2]*_l[_Ndim+jdim+2]*_l[_Ndim+jdim+2]+_r[_Ndim+idim+2]*_r[_Ndim+jdim+2]*_r[_Ndim+jdim+2])*_vec_normal[idim]);// m_(uin*uj^2)
				A5[_Ndim+1] -= m_u2[idim]*m_u2[jdim]*m_u2[jdim]*_vec_normal[idim];// -m_uin * m_ujn^2
				A5[_Ndim+1] -= 0.5*(m_u2[idim]*_vec_normal[idim]*0.5*(_l[_Ndim+jdim+2]*_l[_Ndim+jdim+2]+_r[_Ndim+jdim+2]*_r[_Ndim+jdim+2]));//-0.5*m_uin*m_uj^2
			}
		}
		// final A5
		double coef_e1, coef_e2;
		coef_e1 = - eta_varrho_2n*m_alp1*m_rho2*b1;
		coef_e2 = - eta_varrho_2n*(1-m_alp1)*m_rho1*b2;
		for (int idim=0; idim<_Ndim; idim++){
			coef_e1 += (0.5*(_l[0]*rho1_l*_l[idim+2]+_r[0]*rho1_r*_r[idim+2]) - m_alp_u1[idim]*m_P*b1/m_rho1)*_vec_normal[idim];
			coef_e2 += (0.5*((1-_l[0])*rho2_l*_l[_Ndim+idim+2]+(1-_r[0])*rho2_r*_r[_Ndim+idim+2])-m_alp_u2[idim]*m_P*b2/m_rho2)*_vec_normal[idim];
		}
		for (int i =0; i< _nVar; i++){
			A5[i] += coef_e1*Delta_e1[i] + coef_e2*Delta_e2[i];
		}
		/** Début remplissage matrice décentrement staggered **/
		//lignes de masse
		for(int i=0; i<_nVar*_nVar;i++)
			_absAroe[i]=0;

		//Mass lines
		for(int idim=0; idim<_Ndim;idim++)
		{
			_absAroe[1+idim]=_vec_normal[idim];
			_absAroe[1+idim+_Ndim+1]=0;
			_absAroe[(_Ndim+1)*_nVar+1+idim]=0;
			_absAroe[(_Ndim+1)*_nVar+1+idim+_Ndim+1]=_vec_normal[idim];
		}
		//Contribution of convection (rho u\times u) in the momentum equations
		for(int idim=0; idim<_Ndim;idim++)
			for (int jdim=0; jdim<_Ndim;jdim++){
				// momentum gas (neglect delta alpha and delta P)
				_absAroe[ (1+idim)*_nVar] += (0.5*(_l[2+idim]*_l[2+jdim]+_r[2+idim]*_r[2+jdim])-2*m_u1[idim]*m_u1[jdim])*_vec_normal[jdim];
				_absAroe[(1+idim)*_nVar+jdim+1] += m_u1[idim]*_vec_normal[jdim];
				_absAroe[(1+idim)*_nVar+idim+1] += m_u1[jdim]*_vec_normal[jdim];
				// momentum liquid (neglect delta alpha and delta P)
				_absAroe[(2+_Ndim+idim)*_nVar+_Ndim+1] += (0.5*(_l[_Ndim+2+idim]*_l[_Ndim+2+jdim]+_r[_Ndim+2+idim]*_r[_Ndim+2+jdim])-2*m_u2[idim]*m_u2[jdim])*_vec_normal[jdim];
				_absAroe[(2+_Ndim+idim)*_nVar+_Ndim+1+jdim+1] += m_u2[idim]*_vec_normal[jdim];
				_absAroe[(2+_Ndim+idim)*_nVar+_Ndim+1+idim+1] += m_u2[jdim]*_vec_normal[jdim];
			}
		// contribution of interfacial pressure in momentum equation
		/*
		 * (alpha *rho2*varrho_2+dpi1*(1-alpha)*inv_a2_2*varrho_2)*_vec_normal[idim]
		 */
		for (int idim=0; idim<_Ndim; idim++){
			_absAroe[ (1+idim)*_nVar] += dpi1*varrho_2*(1-m_alp1)*inv_a2_2*_vec_normal[idim];
			_absAroe[ (1+idim)*_nVar+_Ndim+1] += -dpi1*varrho_2*(m_alp1)*inv_a1_2*_vec_normal[idim];
			_absAroe[(2+_Ndim+idim)*_nVar] += - dpi2*varrho_2*(1-m_alp1)*inv_a2_2*_vec_normal[idim];
			_absAroe[(2+_Ndim+idim)*_nVar+_Ndim+1] += dpi2*varrho_2*(m_alp1)*inv_a1_2*_vec_normal[idim];
			for  (int i=0; i<_nVar; i++){
				_absAroe[ (1+idim)*_nVar+i]+=dpi1*varrho_2*(-m_alp1*inv_a2_2*b1*Delta_e1[i]+(1-m_alp1)*inv_a1_2*b2*Delta_e2[i])*_vec_normal[idim];
				_absAroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i]+=-dpi2*varrho_2*(-m_alp1*inv_a2_2*b1*Delta_e1[i]+(1-m_alp1)*inv_a1_2*b2*Delta_e2[i])*_vec_normal[idim];
			}
		}
		// contribution of pressure gradient in momentum equation
		for (int idim=0; idim<_Ndim; idim++){
			_absAroe[ (1+idim)*_nVar] -= alpha*varrho_2*m_rho2*_vec_normal[idim];
			_absAroe[ (1+idim)*_nVar+_Ndim+1] -=alpha* varrho_2*m_rho1*_vec_normal[idim];
			_absAroe[(2+_Ndim+idim)*_nVar] -=  (1-alpha)*varrho_2*m_rho2*_vec_normal[idim];
			_absAroe[(2+_Ndim+idim)*_nVar+_Ndim+1] -= (1-alpha)* varrho_2*m_rho1*_vec_normal[idim];
			for  (int i=0; i<_nVar; i++){
				_absAroe[ (1+idim)*_nVar+i]-=alpha*varrho_2*(-m_alp1*m_rho2*b1*Delta_e1[i] -(1-m_alp1)*m_rho1*b2*Delta_e2[i])*_vec_normal[idim];
				_absAroe[(_Ndim+1)*_nVar+ (1+idim)*_nVar+i]-=(1-alpha)*varrho_2*(-m_alp1*m_rho2*b1*Delta_e1[i] -(1-m_alp1)*m_rho1*b2*Delta_e2[i])*_vec_normal[idim];
			}
		}
		// last row (total energy) To be changed soon
		for (int i=0; i<_nVar; i++){
			_Aroe[(2*_Ndim+2)*_nVar +i] = A5[i];
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
	if(_spaceScheme ==upwind || _spaceScheme==pressureCorrection || _spaceScheme ==lowMach)//calcul de la valeur absolue
	{		//on ordonne les deux premieres valeurs
		if(valeurs_propres_dist[1].real()<valeurs_propres_dist[0].real())
		{
			tmp=valeurs_propres_dist[0];
			valeurs_propres_dist[0]=valeurs_propres_dist[1];
			valeurs_propres_dist[1]=tmp;
		}
		vector< complex< double > > y (taille_vp,0);
		for( int i=0 ; i<taille_vp ; i++)
			y[i] = Poly.abs_generalise(valeurs_propres_dist[i]);

		if(_entropicCorrection)
		{
			double entShift0 = 0;
			double entShift1 = 0;
			entropicShift(_vec_normal,entShift0,entShift1);
			//cout<<"entShift0= "<<entShift0<<endl;
			for( int i=0 ; i<taille_vp ; i++)
			{
				//cout<<"y["<<i<<"]="<<y[i].real()<<endl;
				y[i] += max(entShift0,entShift1);
			}
		}

		if(_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<"valeurs propres"<<endl;
			for( int i=0 ; i<taille_vp ; i++)
				cout<<valeurs_propres_dist[i] <<", "<<endl;
			cout<<"valeurs à atteindre"<<endl;
			for( int i=0 ; i<taille_vp ; i++)
				cout<<y[i] <<", "<<endl;
		}
		Poly.abs_par_interp_directe(taille_vp,valeurs_propres_dist, _Aroe, _nVar,_precision, _absAroe,y);

		if( _spaceScheme ==pressureCorrection){
			for( int i=0 ; i<_Ndim ; i++)
				for( int j=0 ; j<_Ndim ; j++){
					_absAroe[(1+i)*_nVar+1+j]-=alpha*(valeurs_propres_dist[1].real()-valeurs_propres_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
					_absAroe[(2+_Ndim+i)*_nVar+2+_Ndim+j]-=(1-alpha)*(valeurs_propres_dist[1].real()-valeurs_propres_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
				}
		}
		else if( _spaceScheme ==lowMach){
			double M=max(fabs(u1_n),fabs(u2_n))/maxvploc;
			for( int i=0 ; i<_Ndim ; i++)
				for( int j=0 ; j<_Ndim ; j++){
					_absAroe[(1+i)*_nVar+1+j]-=(1-M)*alpha*(valeurs_propres_dist[1].real()-valeurs_propres_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
					_absAroe[(2+_Ndim+i)*_nVar+2+_Ndim+j]-=(1-M)*(1-alpha)*(valeurs_propres_dist[1].real()-valeurs_propres_dist[0].real())/2*_vec_normal[i]*_vec_normal[j];
				}
		}
	}

	//Calcul de la matrice signe pour VFFC, VFRoe et décentrement des termes source

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
		//_signAroe[_nVar*(_nVar-1)+_nVar-1] = signu;
	}
	else
		throw CdmathException("FiveEqsTwoFluid::convectionMatrices: well balanced option not treated");

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
		cout<<"Matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++){
			for(int j=0; j<_nVar;j++)
				cout<<_Aroe[i*_nVar+j]<<" , ";
			cout<<endl;
		}
		cout<<"Valeur absolue matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++){
			for(int j=0; j<_nVar;j++)
				cout<<_absAroe[i*_nVar+j]<<" , ";
			cout<<endl;
		}
		cout<<"Signe matrice de Roe"<<endl;
		for(int i=0; i<_nVar;i++){
			for(int j=0; j<_nVar;j++)
				cout<<_signAroe[i*_nVar+j]<<" , ";
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

void FiveEqsTwoFluid::jacobianDiff(const int &j, string nameOfGroup)
{
	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_JcbDiff[k] = 0;
	if (_limitField[nameOfGroup].bcType==Wall){
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		VecGetValues(_conservativeVars, _nVar, _idm, _Uj);

		double pression=_Vj[1];//pressure inside
		double T=_Vj[_nVar-1];//temperature outside
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		_JcbDiff[0] = 1;
		_JcbDiff[(1+_Ndim)*_nVar +1+_Ndim] = 1;
		_JcbDiff[_nVar]=_limitField[nameOfGroup].v_x[0];
		_JcbDiff[(2+_Ndim)*_nVar +1+_Ndim] =_limitField[nameOfGroup].v_x[1];
		double v2_v=_limitField[nameOfGroup].v_x[0]*_limitField[nameOfGroup].v_x[0];
		double v2_l=_limitField[nameOfGroup].v_x[1]*_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_JcbDiff[2*_nVar]=_limitField[nameOfGroup].v_y[0];
			_JcbDiff[(3+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_y[1];
			v2_v+=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			v2_l+=_limitField[nameOfGroup].v_y[1]*_limitField[nameOfGroup].v_y[1];
			if(_Ndim==3)
			{
				_JcbDiff[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				_JcbDiff[(4+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_z[1];
				v2_v+=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				v2_l+=_limitField[nameOfGroup].v_z[1]*_limitField[nameOfGroup].v_z[1];
			}
		}
		_JcbDiff[(_nVar-1)*_nVar]=         _fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho_v)+0.5*v2_v;
		_JcbDiff[(_nVar-1)*_nVar +1+_Ndim]=_fluides[1]->getInternalEnergy(_limitField[nameOfGroup].T,rho_l)+0.5*v2_l;
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){
		/*
		_JcbDiff[0] = 1;
		_JcbDiff[(1+_Ndim)*_nVar +1+_Ndim] = 1;
		_JcbDiff[_nVar]=_limitField[nameOfGroup].v_x[0];
		_JcbDiff[(2+_Ndim)*_nVar +1+_Ndim] =_limitField[nameOfGroup].v_x[1];
		double v2_v=_limitField[nameOfGroup].v_x[0]*_limitField[nameOfGroup].v_x[0];
		double v2_l=_limitField[nameOfGroup].v_x[1]*_limitField[nameOfGroup].v_x[1];
		if(_Ndim>1)
		{
			_JcbDiff[2*_nVar]=_limitField[nameOfGroup].v_y[0];
			_JcbDiff[(3+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_y[1];
			v2_v+=_limitField[nameOfGroup].v_y[0]*_limitField[nameOfGroup].v_y[0];
			v2_l+=_limitField[nameOfGroup].v_y[1]*_limitField[nameOfGroup].v_y[1];
			if(_Ndim==3)
			{
				_JcbDiff[3*_nVar]=_limitField[nameOfGroup].v_z[0];
				_JcbDiff[(4+_Ndim)*_nVar +1+_Ndim]= _limitField[nameOfGroup].v_z[1];
				v2_v+=_limitField[nameOfGroup].v_z[0]*_limitField[nameOfGroup].v_z[0];
				v2_l+=_limitField[nameOfGroup].v_z[1]*_limitField[nameOfGroup].v_z[1];
			}
		}
		_JcbDiff[(_nVar-1)*_nVar]=         _fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho_v)+0.5*v2_v;
		_JcbDiff[(_nVar-1)*_nVar +1+_Ndim]=_fluides[1]->getInternalEnergy(_limitField[nameOfGroup].T,rho_l)+0.5*v2_l;
		 */
	} else if (_limitField[nameOfGroup].bcType==Outlet){
		//extraction de l etat courant et primitives
		/*
		_idm[0] = j*_nVar;
		for(k=1; k<_nVar;k++)
		{_idm[k] = _idm[k-1] + 1;}
		VecGetValues(_conservativeVars, _nVar, _idm, _phi);
		VecGetValues(_primitiveVars, _nVar, _idm, _externalStates);
		 */
	}
	else if (_limitField[nameOfGroup].bcType!=Neumann && _limitField[nameOfGroup].bcType!=InletPressure){
		cout<<"Condition  limite non traitee pour le bord "<<nameOfGroup<< endl;
		throw CdmathException("FiveEqsTwoFluid::jacobianDiff: Condition  limite non traitee");
	}
}


void FiveEqsTwoFluid::setBoundaryState(string nameOfGroup, const int &j,double *normale){
	//To do controle signe des vitesses pour CL entree/sortie
	int k;
	double v1_2=0, v2_2=0, q1_n=0, q2_n=0, u1_n=0, u2_n=0;//quantités de mouvement et vitesses normales à la face limite;
	double v1[_Ndim], v2[_Ndim];

	_idm[0] = _nVar*j;
	for(k=1; k<_nVar; k++)
		_idm[k] = _idm[k-1] + 1;

	VecGetValues(_conservativeVars, _nVar, _idm, _externalStates);//On initialise l'état fantôme avec l'état interne

	for(k=0; k<_Ndim; k++){
		q1_n+=_externalStates[(k+1)]*normale[k];
		q2_n+=_externalStates[(k+1+1+_Ndim)]*normale[k];
		u1_n+=_Vj[(k+2)]*normale[k];
		u2_n+=_Vj[(k+2+_Ndim)]*normale[k];
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
		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		for(k=0; k<_Ndim; k++){
			_externalStates[(k+1)]-= 2*q1_n*normale[k];
			_externalStates[(k+1+1+_Ndim)]-= 2*q2_n*normale[k];
			_Vj[(k+2)]-= 2*u1_n*normale[k];
			_Vj[(k+2+_Ndim)]-= 2*u2_n*normale[k];
		}
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Vext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Vext, _nVar, _idm, _Vj, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Vext);

		//Pour la diffusion, paroi à vitesses et temperature imposees
		double pression=_Vj[1];//pressure inside
		double T=_Vj[_nVar-1];//temperature outside
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
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
		_externalStates[_nVar-1] = _externalStates[0]*(_fluides[0]->getInternalEnergy(_limitField[nameOfGroup].T,rho_v) + v1_2/2)
																																													+_externalStates[1+_Ndim]*(_fluides[1]->getInternalEnergy(_limitField[nameOfGroup].T,rho_l) + v2_2/2);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Neumann){
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);

		VecAssemblyBegin(_Vext);
		VecSetValues(_Vext, _nVar, _idm, _Vj, INSERT_VALUES);
		VecAssemblyEnd(_Vext);

		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==Inlet){
		_idm[0] = _nVar*j;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecGetValues(_primitiveVars, _nVar, _idm, _Vj);
		double alpha=_limitField[nameOfGroup].alpha;//void fraction outside
		double pression=_Vj[1];//pressure inside
		double T=_limitField[nameOfGroup].T;//temperature outside
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		//cout<<"Inlet alpha= "<<alpha<<" pression= "<<pression<<" temperature= "<<T<<" velocity gas "<<_limitField[nameOfGroup].v_x[0]<<" velocity liq "<<_limitField[nameOfGroup].v_x[1]<<endl;

		_externalStates[0]=alpha*rho_v;
		_externalStates[1 + _Ndim] = (1-alpha)*rho_l;
		v1[0] = _limitField[nameOfGroup].v_x[0];
		v2[0] = _limitField[nameOfGroup].v_x[1];
		if (_Ndim >1){
			v1[1] = _limitField[nameOfGroup].v_y[0];
			v2[1] = _limitField[nameOfGroup].v_y[1];
		}
		if (_Ndim >2){
			v1[2] = _limitField[nameOfGroup].v_z[0];
			v2[2] = _limitField[nameOfGroup].v_z[1];
		}
		for (int idim=0;idim<_Ndim;idim++){
			_externalStates[1 + idim] = v1[idim]* _externalStates[0]; // phase 1
			_externalStates[2 + _Ndim + idim] = v2[idim]* _externalStates[1+_Ndim]; // phase 2
			v1_2 += v1[idim]*v1[idim];
			v2_2 += v2[idim]*v2[idim];
			_Vj[2+idim] = v1[idim];
			_Vj[2+_Ndim+idim] = v2[idim];
		}
		_externalStates[_nVar-1] = _externalStates[0]      *(_fluides[0]->getInternalEnergy(T,rho_v) + v1_2/2)
    																																									+_externalStates[1+_Ndim]*(_fluides[1]->getInternalEnergy(T,rho_l) + v2_2/2);

		// _Vj external primitives
		_Vj[0] = alpha;
		_Vj[_nVar-1] = T;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;

		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Vext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Vext, _nVar, _idm, _Vj, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Vext);
		VecAssemblyEnd(_Uextdiff);
	}
	else if (_limitField[nameOfGroup].bcType==InletPressure){
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
		double alpha=_limitField[nameOfGroup].alpha;
		double pression=_limitField[nameOfGroup].p+hydroPress;
		double T=_limitField[nameOfGroup].T;
		double rho_v=_fluides[0]->getDensity(pression,T);
		double rho_l=_fluides[1]->getDensity(pression,T);
		_externalStates[0]=alpha*rho_v;
		_externalStates[1+_Ndim]=(1-alpha)*rho_l;

		for(int idim=0;idim<_Ndim;idim++){
			_externalStates[idim+1]=_externalStates[0]*_Vj[idim+2];
			_externalStates[idim+2+_Ndim]=_externalStates[1+_Ndim]*_Vj[idim+2+_Ndim];
			v1_2+=_Vj[2+idim]*_Vj[2+idim];
			v2_2+=_Vj[2+_Ndim+idim]*_Vj[2+_Ndim+idim];
		}
		_externalStates[_nVar-1]=    alpha *rho_v*(_fluides[0]->getInternalEnergy(T,rho_v)+v1_2/2)
				                																																				+(1-alpha)*rho_l*(_fluides[1]->getInternalEnergy(T,rho_l)+v2_2/2);

		// _Vj external primitives
		_Vj[0] = alpha;
		_Vj[1] = pression;
		_Vj[_nVar-1] = T;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Vext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Vext, _nVar, _idm, _Vj, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Vext);
		VecAssemblyEnd(_Uextdiff);
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
		double alpha=_Vj[0];
		//cout<<"Outlet alpha= "<<alpha<<" pression int= "<<pression_int<<" pression ext= "<<pression_ext<<" temperature= "<<T<<" velocity gas "<<_Uj[2]<<" velocity liq "<<_Uj[2+_Ndim]<<endl;
		for(int idim=0;idim<_Ndim;idim++){
			v1_2+=_Vj[2+idim]*_Vj[2+idim];
			v2_2+=_Vj[2+_Ndim+idim]*_Vj[2+_Ndim+idim];
		}
		_externalStates[_nVar-1]=alpha*rho_v_ext*(_fluides[0]->getInternalEnergy(T,rho_v_int)+v1_2/2)+(1-alpha)*rho_l_ext*(_fluides[1]->getInternalEnergy(T,rho_l_int)+v2_2/2);

		// _Vj external primitives
		_Vj[1] = pression_ext;

		_idm[0] = 0;
		for(k=1; k<_nVar; k++)
			_idm[k] = _idm[k-1] + 1;
		VecAssemblyBegin(_Uext);
		VecAssemblyBegin(_Vext);
		VecAssemblyBegin(_Uextdiff);
		VecSetValues(_Uext, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecSetValues(_Vext, _nVar, _idm, _Vj, INSERT_VALUES);
		VecSetValues(_Uextdiff, _nVar, _idm, _externalStates, INSERT_VALUES);
		VecAssemblyEnd(_Uext);
		VecAssemblyEnd(_Vext);
		VecAssemblyEnd(_Uextdiff);
	}
	else{
		cout<<"Boundary condition not set for boundary named "<<nameOfGroup<<endl;
		cout<<"Accepted boundary condition are Neumann, Wall, Inlet, and Outlet"<<endl;
		throw CdmathException("Unknown boundary condition");
	}
}

void FiveEqsTwoFluid::addDiffusionToSecondMember
(		const int &i,
		const int &j,
		bool isBord)
{
	double Tm=_Udiff[_nVar-1];
	double lambdal=_fluides[1]->getConductivity(Tm);
	double lambdav=_fluides[0]->getConductivity(Tm);
	double mu1 =_fluides[0]->getViscosity(Tm);
	double mu2 = _fluides[1]->getViscosity(Tm);

	if(mu1==0 && mu2 ==0 && lambdav==0 && lambdal==0 && _heatTransfertCoeff==0)
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
		lambdal=max(lambdal,_heatTransfertCoeff);//wall nucleate boing -> larger heat transfer

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
	double lambda = alpha*lambdav+(1-alpha)*lambdal;
	//on n'a pas de contribution sur la masse
	_phi[0]=0;
	_phi[_Ndim+1]=0;
	//contribution visqueuse sur la quantite de mouvement
	for(int k=1; k<_Ndim+1; k++)
	{
		_phi[k]         = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*mu1*alpha    *(_Vj[2+k] - _Vi[2+k]);//attention car primitif=alpha p u1 u2 T
		_phi[k+_Ndim+1] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*mu2*(1-alpha)*(_Vj[2+k+_Ndim] - _Vi[2+k+_Ndim]);
	}
	_phi[_nVar-1] = _inv_dxi*2/(1/_inv_dxi+1/_inv_dxj)*lambda  *(_Vj[_nVar-1] - _Vi[_nVar-1]);
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

void FiveEqsTwoFluid::jacobian(const int &j, string nameOfGroup, double * normale)
{
	int k;
	for(k=0; k<_nVar*_nVar;k++)
		_Jcb[k] = 0;//No implicitation at this stage

	// loop of boundary types
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
		throw CdmathException("FiveEqs::jacobianDiff: This boundary condition is not treated");
	}
}


void FiveEqsTwoFluid::primToCons(const double *P, const int &i, double *W, const int &j){
	//P= alpha, p, u1, u2, T
	//W=m1,q1,m2,q2,rhoE =alpha1*rho1*(e1+u1^2/2)+alpha2*rho2*(e2+u2^2/2)
	double alpha=P[i*_nVar];
	double pression=P[i*_nVar+1];
	double temperature=P[i*_nVar+_nVar-1];
	double rho_v=_fluides[0]->getDensity(pression,temperature);
	double rho_l=_fluides[1]->getDensity(pression,temperature);
	double u1_sq=0, u2_sq=0;

	W[j*_nVar] = alpha*rho_v;
	W[j*_nVar+1+_Ndim] = (1-alpha)*rho_l;
	// alpha*rho*u
	for(int k=0; k<_Ndim; k++)
	{
		W[j*_nVar+(k+1)] = W[j*_nVar]*P[i*_nVar+(k+2)];//alpha1*rho1*u1
		W[j*_nVar+(k+1)+1+_Ndim] = W[j*_nVar+1+_Ndim]*P[i*_nVar+(k+2)+_Ndim];//alpha2*rho2*u2
	}
	// total energy
	W[j*_nVar+_nVar-1] = W[j*(_nVar)]* _fluides[0]->getInternalEnergy(temperature,rho_v)+W[j*(_nVar)+1+_Ndim]* _fluides[1]->getInternalEnergy(temperature,rho_l);
	for(int k=0; k<_Ndim; k++){
		u1_sq+=P[i*_nVar+(k+2)]*P[i*_nVar+(k+2)];
		u2_sq+=P[i*_nVar+(k+2)+_Ndim]*P[i*_nVar+(k+2)+_Ndim];
	}
	W[j*_nVar+_nVar-1] += (W[j*_nVar]*u1_sq+W[j*_nVar+1+_Ndim]*u2_sq)*0.5;
}

void FiveEqsTwoFluid::consToPrim(const double *Wcons, double* Wprim,double porosity)//To do: treat porosity
{
	//Wprim= alpha, p, u1, u2, T
	//Wcons=m1,q1,m2,q2,rhoE
	double m_v=Wcons[0];
	double m_l=Wcons[1+_Ndim];
	double q1_sq = 0,q2_sq = 0;
	_minm1=min(m_v,_minm1);
	_minm2=min(m_l,_minm2);

	if(m_v<-_precision || m_l<-_precision){
		_nbMaillesNeg+=1;
	}
	for(int k=0;k<_Ndim;k++){
		q1_sq += Wcons[k+1]*Wcons[k+1];
		q2_sq += Wcons[k+1+1+_Ndim]*Wcons[k+1+1+_Ndim];
	}
	if(Wcons[0]>0)//_precision*_precision*_precision)
		q1_sq /= Wcons[0];	//alpha1 rho1 u1²
	else
		q1_sq = 0;
	if(Wcons[1+_Ndim]>0)//_precision*_precision*_precision)
		q2_sq /= Wcons[1+_Ndim];	//alpha2 rho2 u1²
	else
		q2_sq = 0;
	double rho_m_e_m=Wcons[_nVar-1] -0.5*(q1_sq+q2_sq);
	//calcul de la temperature et de la pression pour une loi stiffened gas
	double temperature= (rho_m_e_m-m_v*static_cast<StiffenedGas*>(_fluides[0])->getInternalEnergy(0)-m_l*static_cast<StiffenedGas*>(_fluides[1])->getInternalEnergy(0))/(m_v*_fluides[0]->constante("cv")+m_l*_fluides[1]->constante("cv"));
	double e_v=static_cast<StiffenedGas*>(_fluides[0])->getInternalEnergy(temperature);
	double e_l=static_cast<StiffenedGas*>(_fluides[1])->getInternalEnergy(temperature);
	double gamma_v=_fluides[0]->constante("gamma");
	double gamma_l=_fluides[1]->constante("gamma");
	double Pinf_v=- gamma_v*_fluides[0]->constante("p0");
	double Pinf_l=- gamma_l*_fluides[1]->constante("p0");
	double a=1;
	double b=-(Pinf_v+m_v*(gamma_v-1)*e_v+Pinf_l+m_l*(gamma_l-1)*e_l);
	double c=Pinf_v*Pinf_l+Pinf_v*m_l*(gamma_l-1)*e_l+ Pinf_l*m_v*(gamma_v-1)*e_v;
	double delta= b*b-4*a*c;
	if(delta<0){
		cout<<"delta= "<<delta<<" <0"<<endl;
		*_runLogFile<<"delta= "<<delta<<" <0"<<endl;
		throw CdmathException("FiveEqsTwoFluid::consToPrim: Failed to compute pressure");
	}
	double pression=(-b+sqrt(delta))/(2*a);
	if (pression < 1){
		cout << "pressure = "<< pression << " < 1 Pa " << endl;
		cout << "Conservative state = ";
		for(int k=0; k<_nVar; k++){
			cout<<Wcons[k]<<", ";
		}
		cout<<endl;
		*_runLogFile << "FiveEqsTwoFluid::consToPrim: Failed to compute pressure = "<< pression << " < 1 Pa " << endl;
		throw CdmathException("FiveEqsTwoFluid::consToPrim: Failed to compute pressure");
	}

	double rho_v=_fluides[0]->getDensity(pression,temperature);
	double alpha=m_v/rho_v;
	Wprim[0]= alpha;
	Wprim[1] =  pression;
	for(int k=0;k<_Ndim;k++){//vitesses
		if(Wcons[0]>0)//_precision*_precision*_precision)
			Wprim[k+2] = Wcons[k+1]/Wcons[0];
		else
			Wprim[k+2] = Wcons[k+2+_Ndim]/Wcons[1+_Ndim];
		if(Wcons[1+_Ndim]>0)//_precision*_precision*_precision)
			Wprim[k+2+_Ndim] = Wcons[k+2+_Ndim]/Wcons[1+_Ndim];
		else
			Wprim[k+2+_Ndim] = Wcons[k+1]/Wcons[0];
	}
	Wprim[_nVar-1] = temperature;
}

void FiveEqsTwoFluid::entropicShift(double* n, double& vpcorr0, double& vpcorr1)
{

	// parameters of function dgeev_ (compute the eigenvalues)
	int LDA, LDVL,LWORK, SDIM,LDVR;
	LDA = _nVar;
	LDVL = _nVar;
	LDVR = _nVar;
	LWORK = 20*_nVar;
	char jobvl[]="N", jobvr[]="N";
	double WORK[LWORK], JacoMat[_nVar*_nVar],egvaReal[_nVar],egvaImag[_nVar],egVectorL[_nVar*_nVar],egVectorR[_nVar*_nVar];
	int info_l = 0, info_r = 0;

	/******** Left: Compute the eigenvalues and eigenvectors of Jacobian Matrix (using lapack)********/
	convectionJacobianMatrix(_l, n);
	for (int i=0; i<_nVar*_nVar; i++){
		JacoMat[i] = _JacoMat[i];
	}
	dgeev_(jobvl, jobvl,  &_nVar,
			JacoMat,&LDA,egvaReal,egvaImag, egVectorL,
			&LDVL,egVectorR,
			&LDVR, WORK,&LWORK,
			&info_l);

	//	/******** Right: Compute the eigenvalues and eigenvectors of Jacobian Matrix (using lapack)********/
	convectionJacobianMatrix(_r, n);
	for (int i=0; i<_nVar*_nVar; i++){
		JacoMat[i] = _JacoMat[i];
	}
	dgeev_(jobvl, jobvl,  &_nVar,
			JacoMat,&LDA,egvaReal,egvaImag, egVectorL,
			&LDVL,egVectorR,
			&LDVR, WORK,&LWORK,
			&info_r);

	if (info_l < 0 || info_r < 0)
	{
		cout<<"Warning FiveEqsTwoFluid::entropicShift: dgeev_ did not compute all the eigenvalues, trying heuristic entropy correction "<<endl;
		double u1l_n=0, u2l_n=0, u1r_n=0, u2r_n=0;
		for (int idim=0; idim<_Ndim; idim++){
			u1l_n= _l[2+idim]      *n[idim];
			u2l_n= _l[2+idim+_Ndim]*n[idim];
			u1r_n= _r[2+idim]      *n[idim];
			u2r_n= _r[2+idim+_Ndim]*n[idim];
		}

		vpcorr0 =max(fabs(u1l_n-u1r_n),fabs(u2l_n-u2r_n));
		vpcorr1 = vpcorr0;
	}
	else
	{
		std::vector< std::complex<double> > eigValuesLeft(_nVar);
		std::vector< std::complex<double> > eigValuesRight(_nVar);
		for(int j=0; j<_nVar; j++){
			eigValuesLeft[j] = complex<double>(egvaReal[j],egvaImag[j]);
			eigValuesRight[j] = complex<double>(egvaReal[j],egvaImag[j]);
		}
		Polynoms Poly;
		int sizeLeft =  Poly.new_tri_selectif(eigValuesLeft, eigValuesLeft.size(), _precision);
		int sizeRight =  Poly.new_tri_selectif(eigValuesRight, eigValuesRight.size(), _precision);
		if (_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<" Eigenvalue of JacoMat Left: " << endl;
			for(int i=0; i<sizeLeft; i++)
				cout<<eigValuesLeft[i] << ", "<<endl;
		}
		if (_verbose && _nbTimeStep%_freqSave ==0)
		{
			cout<<" Eigenvalue of JacoMat Right: " << endl;
			for(int i=0; i<sizeRight; i++)
				cout<<eigValuesRight[i] << ", "<<endl;
		}
		vpcorr0 = 0;
		for (int i=1; i<min(sizeLeft,sizeRight)-1; i++)
			vpcorr0 = max(vpcorr0, abs(eigValuesRight[i]-eigValuesLeft[i]));// Kieu
		vpcorr1 = vpcorr0;
	}
}

void FiveEqsTwoFluid::entropicShift(double* n)
{

	// parameters of function dgeev_ (compute the eigenvalues)
	int LDA, LDVL,LWORK, SDIM,LDVR;
	LDA = _nVar;
	LDVL = _nVar;
	LDVR = _nVar;
	LWORK = 20*_nVar;
	char jobvl[]="N", jobvr[]="N";
	double WORK[LWORK], JacoMat[_nVar*_nVar],egvaReal[_nVar],egvaImag[_nVar],egVectorL[_nVar*_nVar],egVectorR[_nVar*_nVar];
	int info = 0;
	/******** Left: Compute the eigenvalues and eigenvectors of Jacobian Matrix (using lapack)********/
	convectionJacobianMatrix(_l, n);
	for (int i=0; i<_nVar*_nVar; i++){
		JacoMat[i] = _JacoMat[i];
	}
	dgeev_(jobvl, jobvl,  &_nVar,
			JacoMat,&LDA,egvaReal,egvaImag, egVectorL,
			&LDVL,egVectorR,
			&LDVR, WORK,&LWORK,
			&info);
	if (info != 0){
		*_runLogFile<<"FiveEqsTwoFluid::JacoMat: dgeev_ unsuccessful computation of the eigenvalues of Jacobian Matrix (left)"<<endl;
		throw CdmathException(
				"FiveEqsTwoFluid::JacoMat: dgeev_ unsuccessful computation of the eigenvalues of Jacobian Matrix (left)");
	}

	std::vector< std::complex<double> > eigValuesLeft(_nVar,0.);
	for(int j=0; j<_nVar; j++){
		eigValuesLeft[j] = complex<double>(egvaReal[j],egvaImag[j]);
	}
	//	/******** Right: Compute the eigenvalues and eigenvectors of Jacobian Matrix (using lapack)********/
	convectionJacobianMatrix(_r, n);
	for (int i=0; i<_nVar*_nVar; i++){
		JacoMat[i] = _JacoMat[i];
	}
	dgeev_(jobvl, jobvl,  &_nVar,
			JacoMat,&LDA,egvaReal,egvaImag, egVectorL,
			&LDVL,egVectorR,
			&LDVR, WORK,&LWORK,
			&info);
	if (info != 0)
		{
		*_runLogFile<<"FiveEqsTwoFluid::entropicShift: dgeev_ unsuccessful computation of the eigenvalues of Jacobian Matrix (right)"<<endl;
		throw CdmathException(
				"FiveEqsTwoFluid::entropicShift: dgeev_ unsuccessful computation of the eigenvalues of Jacobian Matrix (right)");
		}

	std::vector< std::complex<double> > eigValuesRight(_nVar,0.);
	for(int j=0; j<_nVar; j++){
		eigValuesRight[j] = complex<double>(egvaReal[j],egvaImag[j]);
	}
	Polynoms Poly;
	int sizeLeft =  Poly.new_tri_selectif(eigValuesLeft, eigValuesLeft.size(), _precision);
	int sizeRight =  Poly.new_tri_selectif(eigValuesRight, eigValuesRight.size(), _precision);
	if (_verbose && _nbTimeStep%_freqSave ==0)
	{
		cout<<" Eigenvalue of JacoMat Left: " << endl;
		for(int i=0; i<sizeLeft; i++)
			cout<<eigValuesLeft[i] << ", "<<endl;
		cout<<" Eigenvalue of JacoMat Right: " << endl;
		for(int i=0; i<sizeRight; i++)
			cout<<eigValuesRight[i] << ", "<<endl;
	}
	for (int i=0; i<min(sizeLeft,sizeRight)-1; i++)
		_entropicShift[i]= abs(eigValuesRight[i]-eigValuesLeft[i]);
}

Vector FiveEqsTwoFluid::staggeredVFFCFlux()
{
	if(_spaceScheme!=staggered || _nonLinearFormulation!=VFFC)
		throw CdmathException("IsothermalTwoFluid::staggeredVFFCFlux: staggeredVFFCFlux method should be called only for VFFC formulation and staggered upwinding");
	else//_spaceScheme==staggered
	{
		Vector Fij(_nVar);
		Fij(_nVar-1)=0;
		double alpha_roe = _Uroe[0];//Toumi formula
		// interfacial pressure term (hyperbolic correction)
		double dpi = _Uroe[_nVar];

		double u1ijn=0, u2ijn=0, phialphaq1n=0, phialphaq2n=0,vitesse1n=0,vitesse2n=0;
		for(int idim=0;idim<_Ndim;idim++){//URoe = alpha, p, u1, u2, T, dpi
			u1ijn+=_vec_normal[idim]*_Uroe[2+idim];
			u2ijn+=_vec_normal[idim]*_Uroe[2+_Ndim+idim];
		}
		if(u1ijn>=0)
		{
			for(int idim=0;idim<_Ndim;idim++){
				phialphaq1n+=_vec_normal[idim]*_Ui[1+idim];//phi alpha rho u n
				vitesse1n  +=_vec_normal[idim]*_Vi[2+idim];
			}
			Fij(0)=phialphaq1n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phialphaq1n*_Vi[2+idim]+(alpha_roe*_Vj[1]*_porosityj+dpi*_Vi[0]*_porosityi)*_vec_normal[idim];

			double pressioni=_Vi[1];
			double Temperaturei= _Vi[_nVar-1];
			double rho1=_fluides[0]->getDensity(pressioni,Temperaturei);
			double e1_int=_fluides[0]->getInternalEnergy(Temperaturei,rho1);
			Fij(_nVar-1)+=_Ui[0]*(e1_int+0.5*vitesse1n*vitesse1n+_Vj[1]/rho1)*vitesse1n;
		}
		else
		{
			for(int idim=0;idim<_Ndim;idim++){
				phialphaq2n+=_vec_normal[idim]*_Uj[1+idim];//phi alpha rho u n
				vitesse1n  +=_vec_normal[idim]*_Vj[2+idim];
			}
			Fij(0)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(1+idim)=phialphaq2n*_Vj[2+idim]+(alpha_roe*_Vi[1]*_porosityi+dpi*_Vj[0]*_porosityj)*_vec_normal[idim];

			double pressionj=_Vj[1];
			double Temperaturej= _Vj[_nVar-1];
			double rho1=_fluides[0]->getDensity(pressionj,Temperaturej);
			double e1_int=_fluides[0]->getInternalEnergy(Temperaturej,rho1);
			Fij(_nVar-1)+=_Uj[0]*(e1_int+0.5*vitesse1n*vitesse1n+_Vi[1]/rho1)*vitesse1n;
		}

		if(u2ijn>=0)
		{
			for(int idim=0;idim<_Ndim;idim++){
				phialphaq2n+=_vec_normal[idim]*_Ui[2+_Ndim+idim];//phi alpha rho u n
				vitesse2n  +=_vec_normal[idim]*_Vi[2+idim+_Ndim];
			}
			Fij(1+_Ndim)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+_Ndim+idim)=phialphaq2n*_Vi[2+_Ndim+idim]+((1-alpha_roe)*_Vj[1]*_porosityj-dpi*_Vi[0]*_porosityi)*_vec_normal[idim];

			double pressioni=_Vi[1];
			double Temperaturei= _Vi[_nVar-1];
			double rho2=_fluides[1]->getDensity(pressioni,Temperaturei);
			double e2_int=_fluides[1]->getInternalEnergy(Temperaturei,rho2);
			Fij(_nVar-1)+=_Ui[1+_Ndim]*(e2_int+0.5*vitesse2n*vitesse2n+_Vj[1]/rho2)*vitesse2n;
		}
		else
		{
			for(int idim=0;idim<_Ndim;idim++){
				phialphaq2n+=_vec_normal[idim]*_Uj[2+_Ndim+idim];//phi alpha rho u n
				vitesse2n  +=_vec_normal[idim]*_Vj[2+idim+_Ndim];
			}
			Fij(1+_Ndim)=phialphaq2n;
			for(int idim=0;idim<_Ndim;idim++)
				Fij(2+_Ndim+idim)=phialphaq2n*_Vj[2+_Ndim+idim]+((1-alpha_roe)*_Vi[1]*_porosityi-dpi*_Vj[0]*_porosityj)*_vec_normal[idim];

			double pressionj=_Vj[1];
			double Temperaturej= _Vj[_nVar-1];
			double rho2=_fluides[1]->getDensity(pressionj,Temperaturej);
			double e2_int=_fluides[1]->getInternalEnergy(Temperaturej,rho2);
			Fij(_nVar-1)+=_Uj[1+_Ndim]*(e2_int+0.5*vitesse2n*vitesse2n+_Vi[1]/rho2)*vitesse2n;
		}
		return Fij;
	}
}

void FiveEqsTwoFluid::applyVFRoeLowMachCorrections(bool isBord, string groupname)
{
	if(_nonLinearFormulation!=VFRoe)
		throw CdmathException("FiveEqsTwoFluid::applyVFRoeLowMachCorrections: applyVFRoeLowMachCorrections method should be called only for VFRoe formulation");
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
				throw CdmathException("FiveEqsTwoFluid::applyVFRoeLowMachCorrections pressure correction order can be only 1 or 2 for five equation two-fluid model");

			double norm_uij=0, uij_n=0, ui_n=0, uj_n=0;//mean velocities
			double rho1 =  _fluides[0]->getDensity(_Uroe[1],_Uroe[_nVar-1]);
			double rho2 =  _fluides[1]->getDensity(_Uroe[1],_Uroe[_nVar-1]);
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
				_Vij[_nVar-1]=_Vi[_nVar-1];
			}
			else{
				_Vij[0]=_Vj[0];
				_Vij[1]=_Vi[1];
				for(int i=0;i<_Ndim;i++)
				{
					_Vij[2+i]		=_Vj[2+i];
					_Vij[2+i+_Ndim] =_Vj[2+i+_Ndim];
				}
				_Vij[_nVar-1]=_Vj[_nVar-1];
			}
			primToCons(_Vij,0,_Uij,0);
		}
	}
}

void FiveEqsTwoFluid::computeScaling(double maxvp)
{
	_blockDiag[0]=1;//alphaScaling;
	_invBlockDiag[0]=1;//_blockDiag[0];
	_blockDiag[1+_Ndim]=1;//-alphaScaling;
	_invBlockDiag[1+_Ndim]=1.0;//_blockDiag[1+_Ndim];
	for(int q=1; q<_Ndim+1; q++)
	{
		_blockDiag[q]=1/maxvp;
		_invBlockDiag[q]=1/_blockDiag[q];
		_blockDiag[q+1+_Ndim]=1/maxvp;
		_invBlockDiag[q+1+_Ndim]=1/_blockDiag[q+1+_Ndim];
	}
	_blockDiag[_nVar - 1]=1/(maxvp*maxvp);//1
	_invBlockDiag[_nVar - 1]=  1./_blockDiag[_nVar - 1] ;// 1.;//
}

void FiveEqsTwoFluid::testConservation()
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
			else if( i == _nVar-1)
				cout << "Energie totale " <<" (J.m^-3): ";
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
		I =  i;
		for(int j=0; j<_Nmailles; j++)
		{
			VecGetValues(_old, 1, &I, &x);//on recupere la valeur du champ
			SUM += x*_mesh.getCell(j).getMeasure();
			VecGetValues(_newtonVariation, 1, &I, &x);//on recupere la variation du champ
			DELTA += x*_mesh.getCell(j).getMeasure();
			I += _nVar;
		}
		if(fabs(SUM)>_precision)
			cout << SUM << ", variation relative: " << fabs(DELTA /SUM)  << endl;
		else
			cout << " a une somme nulle,  variation absolue: " << fabs(DELTA) << endl;
	}
}

void FiveEqsTwoFluid::save(){
	string prim(_path+"/FiveEqsTwoFluidPrim_");
	string cons(_path+"/FiveEqsTwoFluidCons_");
	prim+=_fileName;
	cons+=_fileName;

	PetscInt Ii, lj;
	for (PetscInt i = 0; i < _Nmailles; i++){
		/* j = 0 : void fraction
			   j = 1 : pressure
			   j = 2, 3, 4: velocity phase 1
			   j = 5, 6, 7: velocity phase 2
			   j = 8 : temperature */
		for (int j = 0; j < _nVar; j++){
			Ii = i*_nVar +j;
			VecGetValues(_primitiveVars,1,&Ii,&_VV(i,j));
		}
	}
	if(_saveConservativeField){
		for (long i = 0; i < _Nmailles; i++){
			for (int j = 0; j < _nVar; j++){
				Ii = i*_nVar +j;
				VecGetValues(_conservativeVars,1,&Ii,&_UU(i,j));
			}
		}
		_UU.setTime(_time,_nbTimeStep+1);
	}
	_VV.setTime(_time,_nbTimeStep+1);

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
		if(_saveConservativeField){
			_UU.setInfoOnComponent(0,"Partial_density1");// (kg/m^3)
			_UU.setInfoOnComponent(1,"Momentum1_x");// phase1 (kg/m^2/s)
			if (_Ndim>1)
				_UU.setInfoOnComponent(2,"Momentum1_y");// phase1 (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(3,"Momentum1_z");// phase1  (kg/m^2/s)
			_UU.setInfoOnComponent(1+_Ndim,"Partial_density2");// phase2 (kg/m^3)
			_UU.setInfoOnComponent(2+_Ndim,"Momentum2_x");// phase2 (kg/m^2/s)

			if (_Ndim>1)
				_UU.setInfoOnComponent(3+_Ndim,"Momentum2_y");// phase2 (kg/m^2/s)
			if (_Ndim>2)
				_UU.setInfoOnComponent(4+_Ndim,"Momentum2_z");// phase2 (kg/m^2/s)
			_UU.setInfoOnComponent(_nVar-1,"Total_energy");

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
