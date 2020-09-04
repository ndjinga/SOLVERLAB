#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh " << endl;
	double xinf=0.0;
	double xsup=1;
	double yinf=0.0;
	double ysup=1;
	int nx=50;
	int ny=50;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny);
	int spaceDim = M.getSpaceDimension();

	// physical constants
	vector<double> viscosite(1)  ;
	viscosite[0]= 0.025;

	//	 set the limit field for each boundary
	double fixedWallVelocityX=0;
	double fixedWallVelocityY=0;
	double fixedWallTemperature=273;

	double movingWallVelocityX=1;
	double movingWallVelocityY=0;
	double movingWallTemperature=273;


	SinglePhase  myProblem(Gas,around1bar300K,spaceDim);
	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = 1e5;
	VV_Constant[1] = 0;
	VV_Constant[2] = 0;
	VV_Constant[3] = 273;

	// name output file
	string fileName = "2DLidDrivenCavityStructuredCentered1bar";

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"fixedWall","fixedWall",yinf,ysup,ny,"fixedWall","movingWall");

	//set the boundary conditions
	myProblem.setWallBoundaryCondition("fixedWall", fixedWallTemperature, fixedWallVelocityX, fixedWallVelocityY);
	myProblem.setWallBoundaryCondition("movingWall", movingWallTemperature, movingWallVelocityX, movingWallVelocityY);

	// physical parameters
	myProblem.setViscosity(viscosite);

	// set the numerical method
	myProblem.setNumericalScheme(staggered, Implicit);

	// set the Petsc resolution
	myProblem.setLinearSolver(GMRES,LU,true);

	//Numerical parameters
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 1;
	double maxTime = 100000;
	double precision = 1e-9;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision*1e8,20);
	myProblem.saveVelocity();

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
