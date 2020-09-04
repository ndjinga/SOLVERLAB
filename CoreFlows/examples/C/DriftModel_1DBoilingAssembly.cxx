#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//setting mesh and groups
	cout << "Building a regular grid " << endl;
	double xinf=0.0;
	double xsup=4.2;
	double xinfcore=1.4;
	double xsupcore=2.8;
	
	int nx=3;//50;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Outlet");
	M.setGroupAtPlan(xinf,0,eps,"Inlet");
	int spaceDim = M.getSpaceDimension();

	// setting boundary conditions 
	double inletConc=0;
	double inletVelocityX=1;
	double inletTemperature=565;
	double outletPressure=155e5;

	// setting physical parameters 
	Field heatPowerField=Field("heatPowerField",CELLS, M, 1);
	int nbCells=M.getNumberOfCells();

	for(int i=0;i<nbCells;i++){
		double x=M.getCell(i).x();

		if (x> xinfcore && x< xsupcore)
			heatPowerField[i]=1e8;
		else
			heatPowerField[i]=0;
	}
	heatPowerField.writeVTK("heatPowerField",true);		

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();
	Field VV("Primitive", CELLS, M, nVar);

	// Prepare for the initial condition
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 0.;
	VV_Constant(1) = 155e5;
	for (int idim=0; idim<spaceDim;idim++)
		VV_Constant(2+idim) = 1;
	VV_Constant(nVar-1) = 565;

	//Initial field creation
	cout << "Building initial field " << endl;
	myProblem.setInitialFieldConstant( M, VV_Constant);

	//set the boundary conditions
	myProblem.setInletBoundaryCondition("Inlet",inletTemperature,inletConc,inletVelocityX);
	myProblem.setOutletBoundaryCondition("Outlet", outletPressure,vector<double>(1,xsup));

	// physical parameters
	myProblem.setHeatPowerField(heatPowerField);


	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);

	// name the result file
	string fileName = "DriftModel1DBoilingAssembly";

	// setting numerical parameters
	unsigned MaxNbOfTimeStep =3 ;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 1;
	double precision = 1e-7;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.saveAllFields(true);
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
