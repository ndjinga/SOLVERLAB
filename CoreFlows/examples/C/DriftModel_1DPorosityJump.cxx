#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=4.2;
	int nx=100;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	int spaceDim = M.getSpaceDimension();

	//Initial data
	double initialConc=0;
	double initialVelocityX =1;
	double initialTemperature=600;
	double initialPressure=155e5;

	// physical parameters
	Field porosityField("Porosity", CELLS, M, 1);
	for(int i=0;i<M.getNumberOfCells();i++){
		double x=M.getCell(i).x();
		if (x> (xsup-xinf)/3 && x< 2*(xsup-xinf)/3)
			porosityField[i]=0.5;
		else
			porosityField[i]=1;
	}
	porosityField.writeVTK("PorosityField",true);		


	DriftModel myProblem(around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[1] = initialConc;
	VV_Constant[1] = initialPressure;
	VV_Constant[2] = initialVelocityX;
	VV_Constant[3] = initialTemperature;

	cout << "Building initial data " << endl;

	// generate initial condition
	myProblem.setInitialFieldConstant( spaceDim, VV_Constant, xinf, xsup, nx,"Inlet","Outlet");

	//set the boundary conditions
	myProblem.setInletBoundaryCondition("Inlet",initialTemperature,initialConc,initialVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet",initialPressure,vector<double>(1,xsup));

	// physical parameters
	myProblem.setPorosityField(porosityField);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);
    myProblem.setNonLinearFormulation(VFFC) ;
    
	// name file save
	string fileName = "1DPorosityJumpUpwindWB";


	/* set numerical parameters */
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

