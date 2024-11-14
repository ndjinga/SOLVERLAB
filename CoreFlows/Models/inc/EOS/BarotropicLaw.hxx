#ifndef BAROTROPICLAW_H
#define BAROTROPICLAW_H

#include "math.h"
#include "Fluide.h"

using namespace std;

// A stiffened gas class

/*! \class BarotropicLaw Fluide.hxx "Fluide.hxx"
 *  \brief Class implementing a generalised stiffened gas law between pressure, density and internal energy
 *  \details Provides viscosity, conductivity laws as well as EOS  \f$P=(\gamma - 1) * rho (e(T)-q) - \gamma*p0 \f$
 */
class BarotropicLaw:public Fluide{

 public:
  BarotropicLaw( double a, double gamma);

  double getPressure(const double  rho);
  double vitesseSon(double rho) ;
  double constante(string name);
  double getDensity(double p, double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getInternalEnergy(double T, double rho=0){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getTemperatureFromPressure(double  p, double rho){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getTemperatureFromEnthalpy(const double  h, const double rho){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getEnthalpy(double T, double rho){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getPressure(double  rhoe,const double  rho){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getPressureFromEnthalpy(double  h,const double  rho){ cout << "This function is not used in the barotropic law case" << endl; return 0; };

  double getInverseSoundSpeed(double p, double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDrhoDT_P(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDpDT_rho(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDhDT_P(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDeDT_rho(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDrhoDe_P(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDrhoDP_e(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDpDe_rho(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDrhoDP_h(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDrhoDh_p(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDeDp_h(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };
  double getDeDh_p(double P,double T){ cout << "This function is not used in the barotropic law case" << endl; return 0; };

    private:
  double  _gamma;//Stiffened gas law : P=_a * rho**(_gamma)
  double _a;//coefficient of the barotropic law
};



#endif
