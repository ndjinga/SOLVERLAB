#include "FiveEqsTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh data
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;

	int spaceDim=1;

	double inletVoidFraction=0;
	vector<double>inletVelocityX(2,2);
	double inletTemperature=563;

	double outletPressure=155e5;

	// physical constants
	double heatPower=1e8;	
	int nbPhase=2;

	FiveEqsTwoFluid  myProblem(around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = inletVoidFraction;
	VV_Constant[1] = outletPressure;
	VV_Constant[2] = inletVelocityX[0];
	VV_Constant[3] = inletVelocityX[1];
	VV_Constant[2+spaceDim*nbPhase] = inletTemperature;

	cout << "Building initial data" << endl;
	// generate initial condition
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"inlet","outlet");

	//set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet",inletVoidFraction,inletTemperature,inletVelocityX);
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);

	// physical parameters
	myProblem.setHeatSource(heatPower);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);
	myProblem.setEntropicCorrection(true);

	// name file save
	string fileName = "1DBoilingChannel";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 5;
	double precision = 1e-8;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
