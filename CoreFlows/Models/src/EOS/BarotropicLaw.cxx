#include "BarotropicLaw.hxx"

//Perfect gas EOS with given gamma
BarotropicLaw::BarotropicLaw( double a, double gamma):Fluide()
{
	if(gamma -1<=0)
		throw EosException("BarotropicLaw::BarotropicLaw: gamma<=1");
	_gamma=gamma;
	_a = a;
	cout<<"Barotropic Law P= a * rho**gamma  with parameter" << " gamma= " << _gamma<<endl;
}

double BarotropicLaw::getPressure(const double  rho) {
	return _a * pow(rho, _gamma); 
}

double BarotropicLaw::vitesseSon(double rho) {
	assert(rho>0);  
	if (_gamma > 1.0)
		return sqrt(_a * _gamma * pow(rho,_gamma-1) ); 
	else if (_gamma == 1.0)
		return sqrt(_a ); 
}
  
double BarotropicLaw::constante(string name){
  	if (name == "a"||name == "A")
  		return _a;
  	else if (name == "gamma"||name == "GAMMA")
  		return _gamma;
  	else
		  return Fluide::constante(name);
  }


