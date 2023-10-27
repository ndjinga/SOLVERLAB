#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

	// set the limit field for each boundary
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallTemperature=563;

	double inletConcentration=0;
	double inletVelocityX=0;
	double inletVelocityY=1;
	double inletTemperature=563;

	double outletPressure=155e5;

	// physical constants
	vector<double> gravite(spaceDim,0.) ;
	gravite[1]=-7;
	gravite[0]=7;
	double heatPower=1e8;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	//Prepare for the mesh
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=20;
	int ny=20;

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = 0;
	VV_Constant[1] = 155e5;
	VV_Constant[2] = 0;
	VV_Constant[3] = 1;
	VV_Constant[4] = 563;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"wall","wall",yinf,ysup,ny,"inlet","outlet");

	//set the boundary conditions
	vector<double>pressure_reference_point(2);
	pressure_reference_point[0]=xsup;
	pressure_reference_point[1]=ysup;
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,pressure_reference_point);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConcentration, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

	// set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);

	// name of result file
	string fileName = "DriftModel_2DInclinedBoilingChannel";

	// computation parameters
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 5;
	double precision = 1e-4;

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
