#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Reading mesh with two branches and two forks" << endl;
	Mesh M("resources/BifurcatingFlow2BranchesEqualSections.med");
	cout << "Reading power and coss sectional area fields " << endl;
	Field Sections("resources/BifurcatingFlow2BranchesEqualSections", CELLS,"Section area");
	Field heatPowerField("resources/BifurcatingFlow2BranchesEqualSections", CELLS,"Heat power");

	heatPowerField.writeVTK("heatPowerField");
	Sections.writeVTK("crossSectionPowerField");

	M.getFace(0).setGroupName("Inlet");//z==0
	M.getFace(31).setGroupName("Outlet");//z==4.2
	cout<<"F0.isBorder() "<<M.getFace(0).isBorder()<<endl;
	int meshDim = 1;//M.getSpaceDimension();

	// set the limit values for each boundary
	double inletTemperature =573.;
	double inletVelocityX = 5;
	double outletPressure = 155e5;

	SinglePhase  myProblem(Liquid,around155bars600K,meshDim);
	int nVar = myProblem.getNumberOfVariables();

	//Set heat source
	myProblem.setHeatPowerField(heatPowerField);
	//Set gravity force
	vector<double> gravite(1,-10);
	myProblem.setGravity(gravite);

	//Set section field
	myProblem.setSectionField(Sections);
	// Prepare the initial condition
	Vector VV_Constant(nVar);
	VV_Constant(0) = 155e5;
	VV_Constant(1) = 5;
	VV_Constant(2) = 573;

	cout << "Building initial data " << endl;

	// generate initial condition
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setInletBoundaryCondition("Inlet",inletTemperature,inletVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);

	// name file save
	string fileName = "2BranchesHeatedChannels";

	// parameters calculation
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

