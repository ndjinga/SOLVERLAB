#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	double xinf=0.0;
	double xsup=4.2;
	int nx=50;
	cout << "Building a regular mesh of "<< nx<< " cells " << endl;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");//Neumann
	M.setGroupAtPlan(xinf,0,eps,"Inlet");//
	int spaceDim = M.getSpaceDimension();

	// set the limit field for each boundary
	LimitField limitInlet, limitOutlet;
	map<string, LimitField> boundaryFields;

	limitInlet.T =573.;
	limitInlet.bcType=Inlet;
	limitInlet.v_x = vector<double>(1,5);
	boundaryFields["Inlet"] = limitInlet;

	limitOutlet.bcType=Outlet;
	limitOutlet.p = 155e5;
	boundaryFields["Outlet"] = limitOutlet;

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();
	Field VV("Primitive", CELLS, M, nVar);//Field of primitive unknowns

	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	VV_Constant(0) = 155e5;//pression initiale
	VV_Constant(1) = 5;//vitesse initiale
	VV_Constant(2) = 573;//temperature initiale

	cout << "Number of Phases = " << nbPhase << endl;
	cout << "Construction de la condition initiale ... " << endl;
	//set the initial field
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	//Physical parameters
	double heatPower=1e8;
	myProblem.setHeatSource(heatPower);
	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);

	// name file save
	string fileName = "1DHeatedChannel";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 5;
	double precision = 1e-7;

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

