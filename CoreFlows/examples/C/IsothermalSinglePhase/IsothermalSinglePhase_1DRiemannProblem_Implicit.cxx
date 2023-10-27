#include "IsothermalSinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building Cartesian mesh " << endl;
	double xinf=0.0;
	double xsup=1.0;
	int nx=10;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"LeftBoundary");
	M.setGroupAtPlan(xinf,0,eps,"RightBoundary");
	int spaceDim = M.getSpaceDimension();

	//initial data
	double initialVelocity_Left=1;
	double initialPressure_Left=155e5;
	double initialVelocity_Right=1;
	double initialPressure_Right=155.0001e5;

	IsothermalSinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	Vector VV_Left(nVar),VV_Right(nVar);
	// left and right constant vectors
	VV_Left[0] = initialPressure_Left;
	VV_Left[1] = initialVelocity_Left;

	VV_Right[0] = initialPressure_Right;
	VV_Right[1] = initialVelocity_Right;

	//Initial field creation
	double discontinuity = (xinf+xsup)/2.;

	cout << "Building initial data " << endl;
	Field VV("Primitive", CELLS, M, nVar);

	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

	//set the boundary conditions
	myProblem.setNeumannBoundaryCondition("LeftBoundary");
	myProblem.setNeumannBoundaryCondition("RightBoundary");

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// name file save
	string fileName = "1DRiemannProblem_implicit";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 5;
	double precision = 1e-6; 
	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveConservativeField(true);
	myProblem.setSaveFileFormat(CSV);
	myProblem.setNewtonSolver(precision,20);

	// set display option to monitor the calculation 
	myProblem.setVerbose( true);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();

	//Ecriture des résultats en CSV
	Field pV = myProblem.getPrimitiveField();
	pV.writeCSV(fileName);


	return EXIT_SUCCESS;
}
