#ifndef FLUIDE_H
#define FLUIDE_H

#include <string>
#include "math.h"
#include "CdmathException.hxx"

/*! \class Fluide Fluide.hxx "Fluide.hxx"
 *  \brief Fluid thermodynamics properties
 *  \details Provides pressure, density, temperature, internal energy, enthalpy, viscosity and conductivity 
 */

using namespace std;

class Fluide{
 protected:
  double _mu, _lambda,_Cv, _Cp,_Tref,_gamma,  _p0, _q,_dragCoeff;
 public:
  Fluide(){_mu=0; _lambda=0;_Cv=0;_Cp=0;_Tref=0;_dragCoeff=0;_gamma=0;  _p0=0; _q=0;}
  virtual ~Fluide(){};

  //Stiffened gas equation of state
  double getPressure(double  rhoe,const double  rho) {
  	return (_gamma - 1) * (rhoe - rho*_q) - _gamma*_p0;
  };
  double getPressureFromEnthalpy(double  h,const double  rho) {
  	return (_gamma - 1)/_gamma * rho * (h - _q) - _p0;
  };
  /*For the newton scheme in the IsothermalTwoFluid model */
  double getPressureDerivativeRhoE()  { return _gamma - 1; }
  double getDensityFromEnthalpy(double p, double h)
  {
  	return _gamma*(p+_p0)/((_gamma-1)*(h-_q));
  }
  double vitesseSonEnthalpie(double h) {  return sqrt((_gamma-1)*h);  };
  double vitesseSonTemperature(const double T, const double rho)
  {
  	double h= getEnthalpy(T,rho);
  	return vitesseSonEnthalpie(h);
  }

  double getViscosity(double T) {return _mu;};
  double getConductivity(double T) {return _lambda;};
  double getDragCoeffs(double T) { return _dragCoeff;};
  void setViscosity(double mu) { _mu=mu;};
  void setConductivity(double lambda) { _lambda= lambda;};
  void setDragCoeffs(double dragCoeff) {_dragCoeff=dragCoeff;};
  //return constants gamma, cp, cv, p0, q
  double constante(string name)
  {
  	if(name== "gamma")
  		return _gamma;
  	else if (name == "cv"||name == "Cv")
  		return _Cv;
  	else if (name == "cp"||name == "Cp")
  		return _Cp;
  	else if(name=="p0")
  		return _p0;
  	else if(name=="q")
  		return _q;
  	else
  		throw CdmathException("Unknown constant: "+name);
  }
  virtual double getDensity(double p, double T)=0;
  virtual double getTemperatureFromPressure(const double  p, const double rho)=0;
  virtual double getTemperatureFromEnthalpy(const double  h, const double rho)=0;
  virtual double getInternalEnergy(double T, double rho)=0;
  virtual double getEnthalpy(double T, double rho)=0;

};

#endif
