#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	int spaceDim = M.getSpaceDimension();

	// set the initial field
	double initialPressure=155e5;
	double initialVelocityX=0;
	double initialTemperature=573;

	//set the boundary data for each boundary
	double outletPressure=80e5;
	double wallVelocityX=0;
	double wallTemperature=573;

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();

	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	Vector VV_constant(nVar);
	VV_constant(0) = initialPressure ;
	VV_constant(1) = initialVelocityX;
	VV_constant(2) = initialTemperature	;

	cout << "Building initial data" << endl;
	Field VV("Primitive", CELLS, M, nVar);

	myProblem.setInitialFieldConstant(M,VV_constant);

	//set the boundary conditions
	myProblem.setWallBoundaryCondition("Wall", wallTemperature, wallVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure,vector<double>(1,xsup));

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setEntropicCorrection(true);

	// name file save
	string fileName = "1DDepressurisation";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 5;
	double precision = 1e-8;

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

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
