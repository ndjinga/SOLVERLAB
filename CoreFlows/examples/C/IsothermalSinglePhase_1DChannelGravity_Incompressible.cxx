#include "IsothermalSinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//setting mesh and groups
	cout << "Building a regular grid " << endl;
	double xinf=0.0;
	double xsup=1.;
	int nx=2;//50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Top");
	M.setGroupAtPlan(xinf,0,eps,"Bottom");
	int spaceDim = M.getSpaceDimension();

	// setting physical parameters 
	vector<double> gravite(spaceDim,0.) ;
	gravite[0]=-10;

	IsothermalSinglePhase  myProblem(Liquid,around1bar300K,spaceDim, false);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 1e5;
	for (int idim=0; idim<spaceDim;idim++)
		VV_Constant(1+idim) = 0;

	//Initial field creation
	cout << "Setting initial data " << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setNeumannBoundaryCondition("Top");
	myProblem.setNeumannBoundaryCondition("Bottom");

	// physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(staggered, Implicit);
	myProblem.setLinearSolver(GMRES, LU);
	
	// name the result file
	string fileName = "1DChannelGravity_Incompressible";

	// setting numerical parameters
	unsigned MaxNbOfTimeStep =1 ;
	int freqSave = 1;
	double cfl = 1;
	double maxTime = 1;
	double precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
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
