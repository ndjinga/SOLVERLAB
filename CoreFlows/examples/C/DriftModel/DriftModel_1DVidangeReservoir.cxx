#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv) {
	//setting mesh and groups
	cout << "Building a regular grid " << endl;
	double xinf = 0.0;
	double xsup = 4.2;
	int nx = 2; //50;
	Mesh M(xinf, xsup, nx);
	double eps = 1.E-8;
	M.setGroupAtPlan(xinf, 0, eps, "Outlet");
	M.setGroupAtPlan(xsup, 0, eps, "Inlet");
	int spaceDim = M.getSpaceDimension();

	// setting boundary conditions
	double inletConc = 1;
	double inletTemperature = 300;
	double outletPressure = 1e5;

	double initialConcTop = 1;
	double initialConcBottom = 0.0001;
	double initialVelocityX = 1;
	double initialPressure = 1e5;
	double initialTemperature = 300;

	// setting physical parameters
	vector<double> gravite(spaceDim, 0.);
	gravite[0] = -10;

	DriftModel myProblem(around1bar300K, spaceDim, false);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	Vector VV_top(nVar), VV_bottom(nVar);

// top and bottom vectors
	VV_top[0] = initialConcTop;
	VV_top[1] = initialPressure;
	VV_top[2] = initialVelocityX;
	VV_top[3] = initialTemperature;

	VV_bottom[0] = initialConcBottom;
	VV_bottom[1] = initialPressure;
	VV_bottom[2] = initialVelocityX;
	VV_bottom[3] = initialTemperature;

	//Initial field creation
	cout << "Setting initial data " << endl;
	myProblem.setInitialFieldStepFunction(M, VV_bottom, VV_top, .8, 0);

	//set the boundary conditions
	myProblem.setInletPressureBoundaryCondition("Inlet", outletPressure,inletTemperature, inletConc, vector<double>(1, xinf));
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure,vector<double>(1, xsup));

	// physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);
	myProblem.setNonLinearFormulation(VFFC);

	// name the result file
	string fileName = "Driftmodel_1DVidangeReservoir";

	// setting numerical parameters
	unsigned MaxNbOfTimeStep = 1;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 1;
	double precision = 1e-5;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.usePrimitiveVarsInNewton(true);
	myProblem.saveAllFields(true);
	myProblem.setVerbose(true);
	myProblem.displayConditionNumber();
	myProblem.setSaveFileFormat(CSV);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation " << fileName << " is successful !" << endl;
	else
		cout << "Simulation " << fileName << "  failed ! " << endl;

	cout << "------------ End of calculation -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
