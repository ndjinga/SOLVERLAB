#include "FiveEqsTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building Cartesian mesh " << endl;
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=50;
	int ny=50;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny);
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Wall");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	M.setGroupAtPlan(yinf,1,eps,"Wall");
	M.setGroupAtPlan(ysup,1,eps,"Wall");
	int spaceDim = M.getSpaceDimension();

	// set the limit field for each boundary
	vector<double> wallVelocityX(2,0);
	vector<double> wallVelocityY(2,0);
	double wallTemperature=300;
	
	// physical constants
	vector<double> gravite(spaceDim,0.) ;
	gravite[1]=-7;
	gravite[0]=7;

	FiveEqsTwoFluid  myProblem(around1bar300K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();
	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 0.5;
	VV_Constant(1) = 1e5;
	VV_Constant(2) = 0;
	VV_Constant(3) = 0;
	VV_Constant(4) = 0;
	VV_Constant(5) = 0;
	VV_Constant(6) = wallTemperature;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setWallBoundaryCondition("Wall",wallTemperature,wallVelocityX,wallVelocityY);

	// set physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// name file save
	string fileName = "2DInclinedSedimentation";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 0.1;
	double maxTime = 5;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveVelocity();
	myProblem.displayConditionNumber();
	myProblem.setSaveFileFormat(CSV);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
