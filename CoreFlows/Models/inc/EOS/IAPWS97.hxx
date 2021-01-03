#ifndef IAPWS97_H
#define IAPWS97_H

#include <string>
#include "CdmathException.hxx"
#include "Fluide.h"

/*! \class IAPWS97 IAPWS97.hxx "IAPWS97.hxx"
 *  \brief IAPWS97 release of water and steam thermodynamics properties (freesteam library)
 *  \details Provides pressure, density, temperature, internal energy, enthalpy, viscosity and conductivity 
 */

using namespace std;

//The IAPWS-IF97 law implemented by the freesteam team
class FluideIAPWS97:public Fluide{
 protected:

 public:
  FluideIAPWS97(){}
  virtual ~FluideIAPWS97(){};

  //Stiffened gas equation of state
  double getPressure(double  rhoe,const double  rho) {
  	return (_gamma - 1) * (rhoe - rho*_q) - _gamma*_p0;
  };
  double getPressureFromEnthalpy(double  h,const double  rho) {
  	return (_gamma - 1)/_gamma * rho * (h - _q) - _p0;
  };
  //For the newton scheme in the IsothermalTwoFluid model
  double getPressureDerivativeRhoE()  { return _gamma - 1; }
  double getDensityFromEnthalpy(double p, double h)  {  	return rhomass_phmass(p,h);  }
  double vitesseSonEnthalpie(double h, double p) {  return speed_sound(T_phmass(p,h), p);  };
  double vitesseSonTemperature(const double T, const double rho)  {  	return speed_sound(T,rho);  }

  double getViscosity(double T, rho) {return visc(T, rho);};
  double getConductivity(double T, double p) {return tcond(T,p,getDensity(T,p);};
  //return constants gamma, cp, cv, p0, q
  double constante(string name, double T, double p)
  {
  	if(name== "gamma")
  		return cpmass(T,p)/cvmass(T,p);
  	else if (name == "cv"||name == "Cv")
  		return cvmass(T,p);
  	else if (name == "cp"||name == "Cp")
  		return cpmass(T,p);
  	else
  		throw CdmathException("Unknown constant: "+name);
  }
  double getDensity(double p, double T)  {  	return rhomass(T,p);  }
  double getTemperatureFromPressure(const double  p, const double rho) {  	return T_phmass(p,h);  }
  double getTemperatureFromEnthalpy(const double  h, const double rho) {  	return T_phmass(p,h);  }
  double getInternalEnergy(double T, double rho) {  	return umass(T,rho);  }
  double getEnthalpy(double T, double rho)       {      return hmass(T,rho);  }

};

#endif
