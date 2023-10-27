#include "FiveEqsTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building regular mesh " << endl;
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=20;
	int ny=20;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny);
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Wall");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	M.setGroupAtPlan(yinf,1,eps,"Inlet");
	M.setGroupAtPlan(ysup,1,eps,"Outlet");
	int spaceDim = M.getSpaceDimension();

	// physical constants
	vector<double> gravite(spaceDim,0.) ;
	gravite[1]=-7;
	gravite[0]=7;

	// set the limit field for each boundary
	LimitField limitWall;
	map<string, LimitField> boundaryFields;
	limitWall.bcType=Wall;
	limitWall.T = 563;
	limitWall.v_x = vector<double>(2,0);
	limitWall.v_y = vector<double>(2,0);
	boundaryFields["Wall"]= limitWall;

	LimitField limitInlet;
	limitInlet.bcType=Inlet;
	limitInlet.T = 563;
	limitInlet.alpha = 0;
	limitInlet.v_x = vector<double>(2,0);
	limitInlet.v_y = vector<double>(2,1);
	boundaryFields["Inlet"]= limitInlet;

	LimitField limitOutlet;
	limitOutlet.bcType=Outlet;
	limitOutlet.p = 155e5;
	boundaryFields["Outlet"]= limitOutlet;

	// physical constants
	double heatPower=1e8;

	FiveEqsTwoFluid  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();
	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 0;
	VV_Constant(1) = 155e5;
	VV_Constant(2) = 0;
	VV_Constant(3) = 1;
	VV_Constant(4) = 0;
	VV_Constant(5) = 1;
	VV_Constant(6) = 563;

	//Initial field creation
	cout << "Building initial data " << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	// set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	// name file save
	string fileName = "2DInclinedBoilingChannel";

	//numerical parameters
	myProblem.setNumericalScheme(upwind, Explicit);
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 5;
	double cfl = 0.5;
	double maxTime = 5;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
