#ifndef STIFFENEDGAS_H
#define STIFFENEDGAS_H

#include <string>
#include "math.h"
#include "CdmathException.hxx"
#include "Fluide.h"

/*! \class StiffenedGas StiffenedGas.hxx "StiffenedGas.hxx"
 *  \brief Stiffened Gas law approximating water and steam : P = (gamma - 1) * (rho * e - rho * q) - gamma * p0
 *  \details Provides pressure, density, temperature, internal energy, enthalpy, viscosity and conductivity 
 */

using namespace std;

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
