#include "Fluide.h"
#include <iostream>

//Perfect gas EOS Loi d'etat gaz parfait
StiffenedGas::StiffenedGas( double gamma, double cv, double T_ref, double e_ref): Fluide()
{
	if(gamma -1<=0)
		throw CdmathException("StiffenedGas::StiffenedGas: gamma<1");
	_gamma=gamma;
	_Cv=cv;
	_Cp=_gamma*_Cv;
	_p0=0;
	_Tref=T_ref;
	_e_ref=e_ref;
	_q=0;
	cout<<"Perfect gas EOS P=(gamma - 1) * rho e with parameter"<< " gamma= " << _gamma<<endl;
	cout<<"Linearised internal energy law e(T)=  e_ref+ cv_ref (T-Tref), around temperature Tref= "<< _Tref<<" K, internal energy e_ref= "<<_e_ref<<" J/Kg, with specific heat cv_ref= "<< _Cv<<" J/Kg/K"<<endl;
}
//Stiffened gas fitted using sound speed
StiffenedGas::StiffenedGas(double rho_ref, double p_ref, double T_ref, double e_ref, double c_ref, double cv_ref): Fluide()
{
	//Old formula
	//_gamma=(1+sqrt(1+4*c_ref*c_ref/(cv_ref*T_ref)))/2;
	//New formula
	_e_ref=e_ref;
	_gamma=1+c_ref*c_ref/(_e_ref+p_ref/rho_ref);
	if(_gamma -1<=0)
		throw CdmathException("StiffenedGas::setEOS: gamma<1");
	_Tref=T_ref;
	_Cv=cv_ref;
	_Cp=_gamma*_Cv;
	_p0=((_gamma-1)*rho_ref*_e_ref-p_ref)/_gamma;
	_q=0;
	cout<<"Stiffened gas EOS P=(gamma - 1) * rho (e(T)-q) - gamma*p0 with parameters"<< " p0= " << _p0 << " gamma= " << _gamma<<endl;
	cout<<"Calibrated around pressure= "<<p_ref<< " Pa, density= "<< rho_ref<< " internal energy= " << e_ref << " J/Kg, sound speed= " << c_ref<< endl;
	cout<<"Linearised internal energy law e(T)=  e_ref+ cv_ref (T-Tref), around temperature Tref= "<< _Tref<<" K, internal energy e_ref= "<<_e_ref<<" J/Kg, specific heat cv_ref= "<< _Cv<<" J/Kg/K"<<endl;
}
// Loi d'etat stiffened gas S. Dellacherie
StiffenedGasDellacherie::StiffenedGasDellacherie( double gamma, double p0, double q, double cv_ref)
{
	if(gamma -1<=0)
		throw CdmathException("StiffenedGas::StiffenedGas: gamma<1");
	_gamma=gamma;
	_Cv=cv_ref;
	_Cp=_gamma*_Cv;
	_p0=p0;
	_q=q;
	_Tref=0;
	_h_ref=q;

	cout<<"S. Dellacherie Stiffened gas EOS P=(gamma - 1)/gamma * rho (h(T)-q) - p0 with parameters"<< " p0=" << _p0 << " gamma= " << _gamma<< " q= " << _q<< endl;
	cout<<"Linearised internal energy law h(T)=  h_ref+ cp_ref (T-Tref) around temperature "<< _Tref<<" K, enthalpy "<<_h_ref<<"  J/Kg,, specific heat cp_ref= "<< _Cp<<" J/Kg/K"<<endl;

}
double StiffenedGas::getInternalEnergy(double T, double rho)
{  
	return _Cv*(T-_Tref)+_e_ref;
}
double StiffenedGasDellacherie::getInternalEnergy(double T, double rho)
{
	double h= getEnthalpy(T);//h=e+p/rho=e+(gamma-1)(e-q)-gamma p0/rho=gamma(e- p0/rho)-(gamma-1)q
	return (h+(_gamma-1)*_q)/_gamma+_p0/rho;
}
double StiffenedGas::getEnthalpy(double T, double rho)
{
	double e=getInternalEnergy( T, rho);
	return _gamma*(e-_p0/rho)-(_gamma-1)*_q;
}
double StiffenedGasDellacherie::getEnthalpy(double T, double rho)
{
	return _Cp*(T-_Tref)+_h_ref;
}
double StiffenedGas::getTemperatureFromPressure(const double  p, const double rho)
{
	//p=(gamma-1)rho(e-q)-gamma p0
	double e=_q +(p+_gamma*_p0)/((_gamma-1)*rho);
	return (e-_e_ref)/_Cv + _Tref;
}
double StiffenedGasDellacherie::getTemperatureFromPressure(const double  p, const double rho)
{
	//P=(gamma - 1)/gamma * rho (h(T)-q) - p0
	double h=_q +_gamma*(p+_p0)/((_gamma-1)*rho);
	return (h-_h_ref)/_Cp + _Tref;
}

double StiffenedGas::getTemperatureFromEnthalpy(const double  h, const double rho)
{
	//h=e+p/rho et p=(gamma-1)rho(e-q)-gamma p0
	double e=(h+(_gamma-1)*_q)/_gamma + _p0/rho;
	return (e-_e_ref)/_Cv + _Tref;
}
double StiffenedGasDellacherie::getTemperatureFromEnthalpy(const double  h, const double rho)
{
	return (h-_h_ref)/_Cp + _Tref;
}

double StiffenedGas::getDensity(double p, double T)
{
	return (p+_gamma*_p0)/((_gamma-1)*(getInternalEnergy(T)-_q));
}
double StiffenedGasDellacherie::getDensity(double p, double T)
{
	return _gamma*(p+_p0)/((_gamma-1)*(getEnthalpy(T)-_q));
}

double StiffenedGas::getJumpDensPress(const double e_l, const double e_r){
	double inv_a_2;
    inv_a_2 = (e_l+e_r)/((_gamma-1)*2*e_l*e_r);
	return inv_a_2;
}
double StiffenedGas::getJumpDensInternalEnergy(const double p_l,const double p_r,const double e_l,const double e_r){
	double b, p_inf;
	p_inf = - _gamma*_p0;
	b = (p_inf-0.5*(p_l+p_r))/((_gamma-1)*e_l*e_r);
	return b;
}
double StiffenedGas::getJumpInternalEnergyTemperature(){
	return _Cv;
}
// function to compute partial rho/ partial p (e constant), partial rho/partial e (p constant)
// use to compute the Jacobian matrix
double StiffenedGas::getDiffDensPress(const double e){
	double inv_a_2;
    inv_a_2 = 1/((_gamma-1)*e);
	return inv_a_2;
}
double StiffenedGas::getDiffDensInternalEnergy(const double p,const double e){
	double b, p_inf;
	p_inf = - _gamma*_p0;
	b = (p_inf-p)/((_gamma-1)*e*e);
	return b;
}
double StiffenedGas::getDiffInternalEnergyTemperature(){
	return _Cv;
}
// function to compute partial rho / partial h (p constant), partial rho/ partial p (h constant)
// use to compute p_x stationary
double StiffenedGas::getDiffDensEnthalpyPressconstant(const double p, const double h){
	double  p_inf = - _gamma*_p0;
	return (p_inf - _gamma*p)/((_gamma-1)*h*h);
}
double StiffenedGas::getDiffDensPressEnthalpyconstant(const double h){
	return _gamma/((_gamma-1)*h);
}
