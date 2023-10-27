#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	/* preprocessing: mesh and group creation */
	cout << "Loading unstructured mesh for test SinglePhase_2DLidDrivenCavity_unstructured()" << endl;
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	Mesh M("resources/BoxWithMeshWithTriangularCells.med");
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"wall");
	M.setGroupAtPlan(xinf,0,eps,"wall");
	M.setGroupAtPlan(yinf,1,eps,"wall");
	M.setGroupAtPlan(ysup,1,eps,"MovingWall");
	int spaceDim = M.getSpaceDimension();

	// physical constants
	vector<double> viscosite(1)  ;
	viscosite[0]= 0.025;

	/* set the limit field for each boundary*/
	LimitField limitWall;
	map<string, LimitField> boundaryFields;
	limitWall.bcType=Wall;
	limitWall.T = 273;
	limitWall.p = 1e5;
	limitWall.v_x = vector<double>(1,0);
	limitWall.v_y = vector<double>(1,0);
	limitWall.v_z = vector<double>(1,0);
	boundaryFields["wall"]= limitWall;

	LimitField limitMovingWall;
	limitMovingWall.bcType=Wall;
	limitMovingWall.T = 273;
	limitMovingWall.p = 1e5;
	limitMovingWall.v_x = vector<double>(1,1);
	limitMovingWall.v_y = vector<double>(1,0);
	limitMovingWall.v_z = vector<double>(1,0);
	boundaryFields["MovingWall"]= limitMovingWall;


	SinglePhase  myProblem(Liquid,around1bar300K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 1e5;
	VV_Constant(1) = 0;
	VV_Constant(2) = 0;
	VV_Constant(3) = 273;

	//Initial field creation
	cout << "Setting initial data " << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	// physical parameters
	myProblem.setViscosity(viscosite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// set the Petsc resolution
	myProblem.setLinearSolver(GMRES,ILU,true);

	// name file save
	string fileName = "2DLidDrivenCavity_unstructured";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 5;
	double maxTime = 5;
	double precision = 1e-8;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
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
