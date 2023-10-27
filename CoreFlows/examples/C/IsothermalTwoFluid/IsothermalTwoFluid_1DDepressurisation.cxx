#include "IsothermalTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=1.0;
	int nx=50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	int spaceDim = M.getSpaceDimension();

	// physical constants
	double dHsatl_over_dp=0.05;
	double Psat=85e5;
	double latentHeat=1e6;

	// set the limit field for each boundary
	LimitField limitOutlet, limitWall;
	map<string, LimitField> boundaryFields;
	limitOutlet.bcType=Outlet;
	limitOutlet.p = 80e5;
	boundaryFields["Outlet"] = limitOutlet;

	limitWall.bcType=Wall;
	limitWall.v_x = vector<double>(2,0);
	boundaryFields["Wall"]= limitWall;
	IsothermalTwoFluid  myProblem(around155bars600K,spaceDim);
	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	Vector VV_Constant(nVar);
	VV_Constant(0) = 0.;
	VV_Constant(1) = 155e5;
	VV_Constant(2) = 0;
	VV_Constant(3) = 0;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	//set physical parameters
//	myProblem.setSatPressure( Psat, dHsatl_over_dp);
//	myProblem.setLatentHeat(latentHeat);

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
	double precision = 1e-6;

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

	cout << "------------ End of simulation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
