#ifndef STIFFENEDGAS_H
#define STIFFENEDGAS_H

#include "math.h"
#include "Fluide.h"

using namespace std;

// A stiffened gas class

/*! \class StiffenedGas Fluide.hxx "Fluide.hxx"
 *  \brief Class implementing a generalised stiffened gas law between pressure, density and internal energy
 *  \details Provides viscosity, conductivity laws as well as EOS  \f$P=(\gamma - 1) * rho (e(T)-q) - \gamma*p0 \f$
 */
class StiffenedGas:public CompressibleFluid{
 private:
  double  _e_ref;//Stiffened gas law : P=(gamma - 1) * rho (e(T)-q) - _gamma*_p0
  double _p0;//coefficient of the stiffened gas law
  double _q ;//coefficient of the stiffened gas law
 public:
  StiffenedGas():CompressibleFluid(){_e_ref=0;_p0=0;_q =0;};
  //Perfect gas EOS
  StiffenedGas( double gamma, double cv, double T_ref, double e_ref);
  //Stiffened gas law fitting reference pressure, density and sound speed
  StiffenedGas(double rho_ref, double p_ref, double T_ref, double e_ref, double soundSpeed_ref, double heatCap_ref);

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
  double vitesseSonEnthalpie(double h) { assert(h>0);  return sqrt((_gamma-1)*h);  };
  double getDhDT_P(double P,double T){return _Cp;}// Dh/DT at constant pressure
  double getDeDT_rho(double P,double T){return _Cv;}// De/DT at constant density
  double getDrhoDe_P(double P,double T){ double e=getInternalEnergy(T); return    -(P+_gamma*_p0)/((_gamma-1.)*(e - _q)*(e - _q));};
  double getDrhoDP_e(double P,double T){ return 1./((_gamma-1.)*getInternalEnergy(T));};
  double getDpDe_rho(double P,double T){ return (_gamma-1.)*getDensity(P,T);};
  double getDrhoDT_P(double P,double T){ double e=getInternalEnergy(T); return -_Cv*(P+_gamma*_p0)/((_gamma-1.)*(e - _q)*(e - _q));};
  double getDpDT_rho(double P,double T){ return (_gamma-1.)*getDensity(P,T)*_Cv;};
  double getDrhoDP_h(double P,double T)  { double rho=getDensity( P, T); double h=getEnthalpy(T,rho); return _gamma/((_gamma-1.)*(h - _q));};
  double getDrhoDh_p(double P,double T){ double rho=getDensity( P, T); double h=getEnthalpy(T,rho); return -(P+_p0)*_gamma/((_gamma-1.)*(h - _q)*(h - _q));};
  double getDeDp_h(double P,double T){ double rho=getDensity( P, T); double h=getEnthalpy(T,rho); return ((h+(_gamma-1)*_q)*_gamma*(_p0+P))-_gamma*(h*(_gamma*_p0+P)+(_gamma-1)*P*_q)/(_gamma*(_p0+P)*_gamma*(_p0+P));}
  double getDeDh_p(double P,double T){return (_gamma*_p0+P)/(_gamma*(_p0+P));}

  double getInternalEnergy(double T, double rho=0);//e is a linear function of T
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

  //return constants p0, q or CompressibleFluid constants
  double constante(string name)
  {
  	if (name == "p0"||name == "P0")
  		return _p0;
  	else if (name == "q"||name == "Q")
  		return _q;
  	else
		return CompressibleFluid::constante(name);
  }
};

// S. Dellacherie stiffened gas class

/*! \class StiffenedGasDellacherie Fluide.hxx "Fluide.hxx"
 *  \brief Class implementing a particular stiffened gas law including saturation properties
 *  \details \f$P=(\gamma - 1)/\gamma * rho (h(T)-q) - p0\f$, formally equivalent to \f$P=(\gamma - 1) * rho (e(T)-q) - \gamma*p0 \f$
 */
class StiffenedGasDellacherie:public CompressibleFluid{
 private:
  double _h_ref;//Stiffened gas law according to S. Dellacherie : P=(gamma - 1)/gamma * rho (h(T)-q) - p0
  double _p0;//coefficient of the stiffened gas law
  double _q ;//coefficient of the stiffened gas law
 public:
  StiffenedGasDellacherie():CompressibleFluid(){_h_ref=0;_p0=0;_q =0;};
  /* Loi des gaz raidis avec coefficients imposés suivant S. Dellacherie*/
  StiffenedGasDellacherie( double gamma, double p0, double q, double cv);

  double getInternalEnergy(double T, double rho);//Problem here since e cannot depend on the sole temperature
  double getEnthalpy(double T, double rho=0);//e is a linear function of T
  double getTemperatureFromPressure(double  p, double rho);
  double getTemperatureFromEnthalpy(const double  h, const double rho=0);
  double getDensity(double p, double T);

  //Stiffened gas equation of state
  double getPressure(double  rhoe,const double  rho) {
  	return (_gamma - 1) * (rhoe - rho*_q) - _gamma*_p0;
  }
  double getPressureFromEnthalpy(double  h,const double  rho) {
  	return (_gamma - 1)/_gamma * rho * (h - _q) - _p0;
  }
  /*For the newton scheme in the IsothermalTwoFluid model */
  double getPressureDerivativeRhoE()  { return _gamma - 1; }
  double getDensityFromEnthalpy(double p, double h)
  {
	assert(h-_q>0);
  	return _gamma*(p+_p0)/((_gamma-1)*(h-_q));
  }
  double vitesseSonEnthalpie(double h) { assert(h>0); return sqrt((_gamma-1)*h);  }
  double getDhDT_P(double P,double T){return _Cp;}// Dh/DT at constant pressure
  double getDeDT_rho(double P,double T){return _Cv;}// De/DT at constant density
  double getDrhoDe_P(double P,double T){ double h=getEnthalpy(T); return -_Cp/_Cv*(P+_p0)*_gamma/((_gamma-1.)*(h - _q)*(h - _q));};
  double getDrhoDP_e(double P,double T){ double rho=getDensity(P,T); return 1./((_gamma-1.)*getInternalEnergy(T, rho));};
  double getDpDe_rho(double P,double T){ return (_gamma-1.)/_gamma*getDensity(P,T)*_Cp/_Cv;};
  double getDrhoDT_P(double P,double T){ double h=getEnthalpy(T); return -_Cp*(P+_p0)*_gamma/((_gamma-1.)*(h - _q)*(h - _q));};
  double getDpDT_rho(double P,double T){ return (_gamma-1.)/_gamma*getDensity(P,T)*_Cp;};
  double getDrhoDP_h(double P,double T){ double h=getEnthalpy(T); return _gamma/((_gamma-1.)*(h - _q));};
  double getDrhoDh_p(double P,double T){ double h=getEnthalpy(T); return -(P+_p0)*_gamma/((_gamma-1.)*(h - _q)*(h - _q));};
  double getDeDp_h(double P,double T){ double h=getEnthalpy(T);return ((h+(_gamma-1)*_q)*_gamma*(_p0+P))-_gamma*(h*(_gamma*_p0+P)+(_gamma-1)*P*_q)/(_gamma*(_p0+P)*_gamma*(_p0+P));}
  double getDeDh_p(double P,double T){return (_gamma*_p0+P)/(_gamma*(_p0+P));}
};

#endif
