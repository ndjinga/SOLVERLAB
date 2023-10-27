#include "FiveEqsTwoFluid.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");//Neumann
	M.setGroupAtPlan(xinf,0,eps,"Wall");//
	int spaceDim = M.getSpaceDimension();

	// set the limit field for each boundary
	LimitField limitOutlet, limitWall;
	map<string, LimitField> boundaryFields;
	limitOutlet.bcType=Outlet;
	limitOutlet.p = 100e5;
	boundaryFields["Outlet"] = limitOutlet;

	limitWall.bcType=Wall;
	limitWall.T = 600;
	limitWall.v_x = vector<double>(2,0);
	boundaryFields["Wall"]= limitWall;

	// physical constants
	double latentHeat=1e6;
	double satTemp=618;
	double dHsatl_over_dp=0.05;
	double Psat=85e5;

	FiveEqsTwoFluid  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();

	//Initial field creation
	Vector VV_Constant(nVar);
	VV_Constant(0) = 0.;
	VV_Constant(1) = 155e5;
	for (int idim=0; idim<spaceDim;idim++){
		VV_Constant(2+idim) = 0;
		VV_Constant(2+idim +spaceDim) =0;
	}
	VV_Constant(2+spaceDim*nbPhase) = 600;

	cout << "Number of Phases = " << nbPhase << endl;
	cout << "Building initial data " << endl;

	// generate initial condition
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);
	/* set physical parameters*/
//	myProblem.setLatentHeat(latentHeat);
//	myProblem.setSatPressure( Psat, dHsatl_over_dp);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setEntropicCorrection(true);

	// name file save
	string fileName = "1DDepressurisation";

	// set numerical parameters
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 0.5;
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

	cout << "------------ End of calculation -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
