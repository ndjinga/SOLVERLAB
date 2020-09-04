#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	/*Preprocessing: mesh and group creation*/
	double xinf=0;
	double xsup=1;
	double yinf=0;
	double ysup=1;
	cout << "Loading unstuctured mesh for test SinglePhase_2DHeatDrivenCavity_unstructured" << endl;
	Mesh M("resources/BoxWithMeshWithTriangularCells.med");
	double eps=1.E-6;
	M.setGroupAtPlan(xinf,0,eps,"coldWall");
	M.setGroupAtPlan(xsup,0,eps,"hotWall");
	M.setGroupAtPlan(yinf,1,eps,"coldWall");
	M.setGroupAtPlan(ysup,1,eps,"hotWall");
	int spaceDim = M.getSpaceDimension();

	// physical constants
	vector<double> viscosite(1), conductivite(1);
	viscosite[0]= 8.85e-5;
	conductivite[0]=1000;//transfert de chaleur du à l'ébullition en paroi.
	vector<double> gravite(spaceDim,0.) ;
	gravite[1]=-10;
	gravite[0]=0;

	// set the limit field for each boundary
	LimitField limitColdWall, limitHotWall;
	map<string, LimitField> boundaryFields;
	limitColdWall.bcType=Wall;
	limitColdWall.T = 590;//Temperature de la parois froide
	limitColdWall.v_x = vector<double>(1,0);
	limitColdWall.v_y = vector<double>(1,0);
	boundaryFields["coldWall"]= limitColdWall;

	limitHotWall.bcType=Wall;
	limitHotWall.T = 560;//Temperature des parois chauffantes
	limitHotWall.v_x = vector<double>(1,0);
	limitHotWall.v_y = vector<double>(1,0);
	boundaryFields["hotWall"]= limitHotWall;

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	//Initial field creation
	cout << "Construction de la condition initiale" << endl;
	/* First case constant initial data */
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 155e5;
	VV_Constant(1) = 0;
	VV_Constant(2) = 0;
	VV_Constant(3) = 573;

	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	// physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// set the Petsc resolution
	myProblem.setLinearSolver(GMRES,ILU,true);

	// name result file
	string fileName = "2DHeatDrivenCavity_unstructured";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 1;
	double maxTime = 50;
	double precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,50);
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
