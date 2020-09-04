#include "IsothermalTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	cout << "Building Cartesian mesh " << endl;
	int spaceDim=1;
	double xinf=0.0;
	double xsup=1.0;
	int nx=50;

	vector<double> wallVelocityX(2,0);

	// physical constants
	vector<double> gravite(spaceDim,0.) ;
	gravite[0]=-10;

	IsothermalTwoFluid  myProblem(around1bar300K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar,0.);
	// constant vector
	VV_Constant[0] = 0.5;
	VV_Constant[1] = 1e5;
	VV_Constant[2] = 0;
	VV_Constant[3] = 0;

	//Initial field creation
	cout << "Building initial data " << endl;
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"wall","wall");
	myProblem.setWallBoundaryCondition("wall",wallVelocityX);


	// physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);
	myProblem.setEntropicCorrection(true);

	// name file save
	string fileName = "1DSedimentation";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 1;
	double maxTime = 5;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.displayConditionNumber();
	myProblem.setSaveFileFormat(CSV);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of simulation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
