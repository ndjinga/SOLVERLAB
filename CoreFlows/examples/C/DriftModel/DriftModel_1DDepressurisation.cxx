#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building regular mesh " << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;
	int spaceDim = 1;

	//Initial data
	double initialConc=0;
	double initialVelocityX =0;
	double initialTemperature=600;
	double initialPressure=155e5;

	//Boundary data
	double wallVelocityX=0;
	double wallTemperature=563;
	double outletPressure=50e5;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();

	// Prepare the initial condition
	vector<double> VV_Constant(nVar);
	VV_Constant[1] = initialConc;
	VV_Constant[1] = initialPressure;
	VV_Constant[2] = initialVelocityX;
	VV_Constant[3] = initialTemperature;

	//Initial field creation
	cout << "Building initial data " << endl;
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"wall","Outlet");

	//set the boundary conditions
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet",outletPressure,vector<double>(1,xsup));

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setEntropicCorrection(true);

	// name file save
	string fileName = "1DDepressurisation";

	//Numerical parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 1;
	double maxTime = 1;
	double precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	bool ok;

	// evolution
	myProblem.initialize();

	ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
