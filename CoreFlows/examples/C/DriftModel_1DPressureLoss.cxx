#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	cout << "Building a regular grid " << endl;
	int spaceDim=1;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;

	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");
	M.setGroupAtPlan(xinf,0,eps,"Inlet");

	double inletConc=0;
	double inletVelocityX =1;
	double inletTemperature=563;

	double outletPressure=155e5;

	// physical parameters
	Field pressureLossField("pressureLoss", FACES, M, 1);
	pressureLossField(nx/4)=50;
	pressureLossField(nx/2)=100;
	pressureLossField(3*nx/4)=150;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = inletConc;
	VV_Constant[1] = outletPressure;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletTemperature;

	cout << "Building initial data " << endl;

	// generate initial condition
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"inlet","outlet");

	//set the boundary conditions
	myProblem.setInletBoundaryCondition("inlet",inletTemperature,inletConc,inletVelocityX);
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,vector<double>(1,xsup));

	// physical parameters
	myProblem.setPressureLossField(pressureLossField);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);

	// name file save
	string fileName = "1DPressureLossUpwindWB";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 5;
	double precision = 1e-5;

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
