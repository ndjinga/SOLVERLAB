#ifndef FLUIDE_H
#define FLUIDE_H

#include <string>
#include <assert.h>
#include "EosException.hxx"
#include <iostream>

/*! \class Fluide Fluide.hxx "Fluide.hxx"
 *  \brief Fluid thermodynamics properties (fluid may be compressible or incompressible)
 *  \details Provides fluid density \f$\rho(P,T)\f$,  viscosity \f$\mu(P,T)\f$ and conductivity \f$\lambda(P,T)\f$ laws 
 */

using namespace std;

class Fluide{
 protected:
  double _mu;/* Constant viscosity coeff */
  double _lambda;/* Constant conductity coeff */
  double _dragCoeff;/* Constant drag coefficient */
  double _Tref, _Pref;/* Reference value of the temperature and pressure */
  double _Cv;/* Constant specific heat at constant volume */
  double _Cp;/* Constant specific heat at constant pressure */
 public:
  Fluide( double Pref=0, double Tref=0)
  {
   _Tref=Tref, _Pref=Pref;
   _mu=0; _lambda=0; _dragCoeff=0; 
   _Cv=0; _Cp=0; 
  };
  virtual ~Fluide(){};

  //order 0 derivatives
  virtual double getDensity(double p, double T)=0;
  virtual double getInternalEnergy(double T, double rho=0)=0;
  virtual double getTemperatureFromPressure(double  p, double rho)=0;
  virtual double getTemperatureFromEnthalpy(const double  h, const double rho)=0;
  virtual double getEnthalpy(double T, double rho)=0;
  virtual double getPressure(double  rhoe,const double  rho) = 0;
  virtual double getPressureFromEnthalpy(double  h,const double  rho)  = 0;
  
  double getViscosity(double T=0) {return _mu;};
  double getConductivity(double T=0) {return _lambda;};
  double getDragCoeffs(double T=0) { return _dragCoeff;};
  void setViscosity(double mu) { _mu=mu;};
  void setDragCoeff(double dragCoeff) {_dragCoeff=dragCoeff;};
  void setConductivity(double lambda) { _lambda= lambda;};

  //order 1 derivatives
  virtual double getInverseSoundSpeed(double p, double T)=0;
  virtual double getDrhoDT_P(double P,double T)=0;// Drho/DT at constant pressure (cf->alpha coefficient de dilatation isobare)
  virtual double getDpDT_rho(double P,double T)=0;// Dp/DT at constant density (cf->beta coefficient de compressibilité isochore) ->relates to the internal pressure
  virtual double getDhDT_P(double P,double T)=0;// Dh/DT at constant pressure->specific heat
  virtual double getDeDT_rho(double P,double T)=0;// De/DT at constant density->specific heat
  virtual double getDrhoDe_P(double P,double T)=0;
  virtual double getDrhoDP_e(double P,double T)=0;
  virtual double getDpDe_rho(double P,double T)=0;
  virtual double getDrhoDP_h(double P,double T)=0;
  virtual double getDrhoDh_p(double P,double T)=0;
  double getDpDT_h  (double P,double T){ double rho=getDensity(P,T); return -_Cp/(1/rho+T/(rho*rho)*getDrhoDT_P(P,T));}// Dh/DT at constant pressure
  double getDrhoDT_e(double P,double T){ double rho=getDensity(P,T); return  _Cv*rho*rho/(T*getDpDT_rho(P,T)-P);}// Dh/DT at constant pressure
  virtual double getDeDp_h(double P,double T)=0;
  virtual double getDeDh_p(double P,double T)=0;
  
  //return constants mu, lambda, dragCoeff
  double constante(string name)
  {
  	if(name=="dragCoeff")
  		return _dragCoeff;
  	else if (name == "mu"||name == "viscosity")
  		return _mu;
  	else if (name == "lambda"||name == "conductivity")
  		return _lambda;
  	else if (name == "cv"||name == "Cv")
  		return _Cv;
  	else if (name == "cp"||name == "Cp")
  		return _Cp;
  	else
		throw EosException("Unknown constant: "+name);
  }
};

/*! \class CompressibleFluid Fluide.h "Fluide.h"
 *  \brief Class implementing compressible fluid laws with finite speed of sound (could be stiffened gas or IAPWS)
 *  \details Provides EOS laws that are function of density \f$h(T,rho) \f$, \f$T(P,rho) \f$, \f$T(h,rho) \f$,  \f$e(T,rho) \f$, as well as functions to compute the sound speed
 */
class CompressibleFluid:public Fluide{
 protected:
  double _gamma;
  
 public:
  CompressibleFluid():Fluide()
  { 
    _gamma=0;/* ratio of cp and cv */
  }
  
  //order 0 derivatives
  virtual double getInternalEnergy(double T, double rho=0)=0;
  virtual double getTemperatureFromPressure(double  p, double rho)=0;
  virtual double getTemperatureFromEnthalpy(const double  h, const double rho)=0;
  virtual double getEnthalpy(double T, double rho)=0;
  virtual double getPressure(double  rhoe,const double  rho) = 0;
  virtual double getPressureFromEnthalpy(double  h,const double  rho)  = 0;
  
  /*For the newton scheme in the IsothermalTwoFluid model */
  virtual double getPressureDerivativeRhoE()  = 0;
  virtual double getDensityFromEnthalpy(double p, double h) = 0;

  //order 1 derivatives
  virtual double getDrhoDe_P(double P,double T)=0;
  virtual double getDrhoDP_e(double P,double T)=0;
  virtual double getDpDe_rho(double P,double T)=0;
  virtual double getDrhoDT_P(double P,double T)=0;// Drho/DT at constant pressure (cf->coefficient de dilatation isobare)
  virtual double getDpDT_rho(double P,double T)=0;// Dp/DT at constant density (cf->coefficient de compressibilité isochore) ->relates to the internal pressure
  virtual double getDhDT_P(double P,double T)=0;// Dh/DT at constant pressure->specific heat
  virtual double getDeDT_rho(double P,double T)=0;// De/DT at constant density->specific heat
  virtual double getDrhoDP_h(double P,double T)=0;
  virtual double getDrhoDh_p(double P,double T)=0;
  virtual double getDeDp_h(double P,double T)=0;
  virtual double getDeDh_p(double P,double T)=0;

  double getInverseSoundSpeed(double P, double T)
  {
  	assert(P>0);
  	assert(T>0);
  	return 1./vitesseSonPressure( P, T);
  }
  virtual double vitesseSonEnthalpie(double h) = 0;
  virtual double vitesseSonTemperature(const double T, const double rho)
  {
  	assert(rho>0);
  	assert(T>0);
  	double h= getEnthalpy(T,rho);
  	return vitesseSonEnthalpie(h);
  }
  virtual double vitesseSonPressure(const double P, const double T)
  {
  	assert(P>0);
  	assert(T>0);
  	double rho=getDensity(P, T);
  	return vitesseSonTemperature(T,rho);
  }
  
  //return constants gamma, cp, cv or Fluide class constants
  double constante(string name)
  {
  	if(name== "gamma")
  		return _gamma;
  	else
		return Fluide::constante(name);
  }
};

/*! \class IncompressibleFluid Fluide.h "Fluide.h"
 *  \brief Class implementing incompressible fluid laws with infinite speed of sound
 *  \details The density is constant
 */
class IncompressibleFluid:public Fluide{
 protected:
  double _rho;
  //Reference values
  double _href,_eref;
  
 public:
  IncompressibleFluid(double rho, double pref=0, double Tref=0, double eref=0):Fluide(pref,Tref)
  { 
	_rho=rho;
	std::cout<<"Incompressible fluid with density rho = "<< _rho <<std::endl;
	_eref=eref;
	_href=eref+pref/rho;
  }
  
  double getDensity(double p, double T){return _rho;}
  double getInternalEnergy(double T, double rho=0){ return _eref; }
  double getTemperatureFromPressure(double  p, double rho){ return _Tref;}
  double getTemperatureFromEnthalpy(const double  h, const double rho){ return _Tref;}
  double getEnthalpy(double T, double rho){ return _href;}
  double getPressure(double  rhoe,const double  rho){ return _Pref;}
  double getPressureFromEnthalpy(double  h,const double  rho){ return _Pref;}
  
  //order 1 derivatives
  double getInverseSoundSpeed(double P, double T){return 0;}
  double getDhDT_P(double P,double T){return _Cp;}// Dh/DT at constant pressure
  double getDeDT_rho(double P,double T){return _Cv;}// De/DT at constant density
  double getDrhoDe_P(double P,double T){return 0.;};
  double getDrhoDP_e(double P,double T){return 0.;};
  double getDpDe_rho(double P,double T){return 0.;};//not sure what to return here
  double getDrhoDT_P(double P,double T){return 0.;}// Drho/DT at constant pressure (cf->coefficient de dilatation isobare)
  double getDpDT_rho(double P,double T){return 0.;}//not sure what to return here
  double getDrhoDP_h(double P,double T){return 0.;};
  double getDrhoDh_p(double P,double T){return 0.;};
  double getDeDp_h(double P,double T){return 0.;}//not sure what to return here
  double getDeDh_p(double P,double T){return 0.;}//not sure what to return here
};

#endif
