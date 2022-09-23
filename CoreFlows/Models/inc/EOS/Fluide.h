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
 public:
  Fluide( double Pref=0, double Tref=0)
  {
   _Tref=Tref, _Pref=Pref;
   _mu=0; _lambda=0; _dragCoeff=0; 
  };
  virtual ~Fluide(){};

  double getViscosity(double T=0) {return _mu;};
  double getConductivity(double T=0) {return _lambda;};
  double getDragCoeffs(double T=0) { return _dragCoeff;};
  void setViscosity(double mu) { _mu=mu;};
  void setDragCoeff(double dragCoeff) {_dragCoeff=dragCoeff;};
  void setConductivity(double lambda) { _lambda= lambda;};

  virtual double getDensity(double p, double T)=0;
  virtual double getInverseSoundSpeed(double p, double T)=0;

  //return constants mu, lambda, dragCoeff
  double constante(string name)
  {
  	if(name=="dragCoeff")
  		return _dragCoeff;
  	else if (name == "mu"||name == "viscosity")
  		return _mu;
  	else if (name == "lambda"||name == "conductivity")
  		return _lambda;
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
  double _Cv;/* Constant specific heat at constant volume */
  double _Cp;/* Constant specific heat at constant pressure */
  
 public:
  CompressibleFluid()
  { 
  _Cv=0; _Cp=0; 
  _Tref=0; _Pref=0; 
  _gamma=0;
  }
  
  virtual double getInternalEnergy(double T, double rho=0)=0;
  virtual double getTemperatureFromPressure(double  p, double rho)=0;
  virtual double getTemperatureFromEnthalpy(const double  h, const double rho)=0;
  virtual double getEnthalpy(double T, double rho)=0;
  virtual double getPressure(double  rhoe,const double  rho) = 0;
  virtual double getPressureFromEnthalpy(double  h,const double  rho)  = 0;
  
  /*For the newton scheme in the IsothermalTwoFluid model */
  virtual double getPressureDerivativeRhoE()  = 0;
  virtual double getDensityFromEnthalpy(double p, double h) = 0;
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
  
  double getInverseSoundSpeed(double P, double T)
  {
  	assert(P>0);
  	assert(T>0);
  	return 1./vitesseSonPressure( P, T);
  }
  
  //return constants gamma, cp, cv or Fluide class constants
  double constante(string name)
  {
  	if (name == "cv"||name == "Cv")
  		return _Cv;
  	else if (name == "cp"||name == "Cp")
  		return _Cp;
  	else if(name== "gamma")
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
  
 public:
  IncompressibleFluid(double rho)
  { 
	_rho=rho;
	std::cout<<"Incompressible fluid with density rho = "<< _rho <<std::endl;
  }
  
  double getDensity(double p, double T)
  {
  	return _rho;
  }
  double getInverseSoundSpeed(double P, double T)
  {
  	return 0;
  }
};

#endif
