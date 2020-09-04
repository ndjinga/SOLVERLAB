#include "Fluide.h"
#include <cstdlib>

#include <iostream>

using namespace std;

int main(int argc, char** argv) {
	double _Tsat = 656; //saturation temperature used in Dellacherie EOS
	StiffenedGasDellacherie fluid1 = StiffenedGasDellacherie(1.43, 0,
			2.030255e6, 1040.14); //stiffened gas law for Gas from S. Dellacherie
	StiffenedGasDellacherie fluid2 = StiffenedGasDellacherie(2.35, 1e9,
			-1.167056e6, 1816.2); //stiffened gas law for water from S. Dellacherie

	double P = 155e6;
	double T = 500;
	double h = 0;

	double rho1 = fluid1.getDensity(P, T);
	double Tvalid1 = fluid1.getTemperatureFromPressure(P, rho1);
	double h1 = fluid1.getEnthalpy(T, rho1);

	cout << endl;
	cout << "density fluide 1 = " << rho1 << endl;
	cout << "Tvalid1 fluide 1 = " << Tvalid1 << endl;
	cout << "h1 fluide 1 = " << h1 << endl;

	return  EXIT_SUCCESS;
}
