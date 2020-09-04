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

// A standard stiffened gas class

/*! \class StiffenedGas Fluide.hxx "Fluide.hxx"
 *  \brief Class implementing a standard stiffened gas law between pressure, density and internal energy
 *  \details  
 */
class StiffenedGas:public Fluide{
 private:
  double  _e_ref;//Stiffened gas law : P=(gamma - 1) * rho e(T) - _gamma*_p0
 public:
  StiffenedGas():Fluide(){_e_ref=0;};
  //Perfect gas EOS
  StiffenedGas( double gamma, double cv, double T_ref, double e_ref);
  //Stiffened gas law fitting reference pressure, density and sound speed
  StiffenedGas(double rho_ref, double p_ref, double T_ref, double e_ref, double soundSpeed_ref, double heatCap_ref);

  double getInternalEnergy(double T, double rho=0);
  double getEnthalpy(double T, double rho);
  double getTemperatureFromPressure(double  p, double rho);
  double getTemperatureFromEnthalpy(const double  h, const double rho);
  double getDensity(double p, double T);

  // Functions used to compute the Roe matrix for the five equation model (Kieu)
  /* get differential of the density rho = rho(P,e)
   * wrt the pressure (const e) wrt the internal energy (const P) */
  double getJumpDensPress(const double e_l, const double e_r);
  double getJumpDensInternalEnergy(const double p_l,const double p_r,const double e_l,const double e_r);
  double getJumpInternalEnergyTemperature();
  double getDiffDensPress(const double e);
  double getDiffDensInternalEnergy(const double p,const double e);
  double getDiffInternalEnergyTemperature();
  /* get differential of the density rho = rho(P,h)
     * wrt the pressure (const h) wrt the enthalpy (const P) */
  double getDiffDensEnthalpyPressconstant(const double p, const double h);
  double getDiffDensPressEnthalpyconstant(const double h);
};

// S. Dellacherie stiffened gas class

/*! \class StiffenedGasDellacherie Fluide.hxx "Fluide.hxx"
 *  \brief Class implementing a particular stiffened gas law including saturation properties
 *  \details
 */
class StiffenedGasDellacherie:public Fluide{
 private:
  double _h_ref;//Stiffened gas law according to S. Dellacherie : P=(gamma - 1) * rho (e(T)-q) - _gamma*_p0
 public:
  StiffenedGasDellacherie():Fluide(){_h_ref=0;};
  /* Loi des gaz raidis avec coefficients imposés suivant S. Dellacherie*/
  StiffenedGasDellacherie( double gamma, double p0, double q, double cv);

  double getInternalEnergy(double T, double rho);
  double getEnthalpy(double T, double rho=0);
  double getTemperatureFromPressure(double  p, double rho);
  double getTemperatureFromEnthalpy(const double  h, const double rho=0);
  double getDensity(double p, double T);

 };
#endif
