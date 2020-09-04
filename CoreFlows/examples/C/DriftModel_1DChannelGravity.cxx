#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//setting mesh and groups
	cout << "Building a regular grid " << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;//50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");
	M.setGroupAtPlan(xinf,0,eps,"Inlet");
	int spaceDim = M.getSpaceDimension();

	// setting boundary conditions 
	double inletConc=0;
	double inletVelocityX=1;
	double inletEnthalpy=1.3e6;
	double outletPressure=155e5;

	// setting physical parameters 
	vector<double> gravite(spaceDim,0.) ;
	gravite[0]=-10;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 0.;
	VV_Constant(1) = 155e5;
	for (int idim=0; idim<spaceDim;idim++)
		VV_Constant(2+idim) = 1;
	VV_Constant(nVar-1) = 578;

	//Initial field creation
	cout << "Setting initial data " << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setInletEnthalpyBoundaryCondition("Inlet",inletEnthalpy,inletConc,inletVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure,vector<double>(1,xsup));

	// physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);
	myProblem.setWellBalancedCorrection(true);
	myProblem.setNonLinearFormulation(VFRoe);

	// name the result file
	string fileName = "Driftmodel_1DChannelGravity";

	// setting numerical parameters
	unsigned MaxNbOfTimeStep =3 ;
	int freqSave = 1;
	double cfl = 100;
	double maxTime = 1;
	double precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.usePrimitiveVarsInNewton(true);
	myProblem.saveAllFields(true);
	myProblem.displayConditionNumber();
	myProblem.setSaveFileFormat(CSV);

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
