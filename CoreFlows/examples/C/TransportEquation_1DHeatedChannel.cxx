#include "TransportEquation.hxx"

using namespace std;


int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	double xinf=0.0;
	double xsup=4.2;
	int nx=10;
	cout << "Building a 1D mesh with "<<nx<<" cells" << endl;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Neumann");
	M.setGroupAtPlan(xinf,0,eps,"Inlet");
	int spaceDim = M.getSpaceDimension();

	// Boundary conditions
	map<string, LimitFieldTransport> boundaryFields;

	LimitFieldTransport limitNeumann;
	limitNeumann.bcType=NeumannTransport;
	boundaryFields["Neumann"] = limitNeumann;

	LimitFieldTransport limitInlet;
	limitInlet.bcType=InletTransport;
	limitInlet.h =1.3e6;//Inlet water enthalpy
	boundaryFields["Inlet"] = limitInlet;

	//Set the fluid transport velocity
	vector<double> transportVelocity(1,5);//fluid velocity vector

	TransportEquation  myProblem(LiquidPhase,around155bars600KTransport,transportVelocity);
	Field VV("Enthalpy", CELLS, M, 1);

	//Set rod temperature and heat exchamge coefficient
	double rodTemp=623;//Rod clad temperature
	double heatTransfertCoeff=1000;//fluid/solid exchange coefficient 
	myProblem.setRodTemperature(rodTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);

	//Initial field creation
	Vector VV_Constant(1);//initial enthalpy
	VV_Constant(0) = 1.3e6;

	cout << "Building the initial data " << endl;

	// generate initial condition
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	// set the numerical method
	myProblem.setTimeScheme( Explicit);

	// name result file
	string fileName = "1DFluidEnthalpy";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
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

	// set display option to monitor the calculation
	bool computation=true;
	bool system=true;
	myProblem.setVerbose( computation, system);
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
