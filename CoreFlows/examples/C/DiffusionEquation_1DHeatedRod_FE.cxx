#include "DiffusionEquation.hxx"

using namespace std;

#define PI 3.14159265

void power_field_diffusionTest(Field & Phi){
	double L=4.2;
	double lambda=0.2;
	double phi=1e5;
	double x;
	Mesh M = Phi.getMesh();
	int nbNodes = M.getNumberOfNodes();
	for (int j = 0; j < nbNodes; j++) {
		x=M.getNode(j).x();
		Phi(j) = phi*cos(PI*(x-L/2)/(L+lambda));
	}
}

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	double xinf=0.0;
	double xsup=4.2;
	int nx=10;
	cout << "Building of a 1D mesh with "<<nx<<" cells" << endl;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Neumann");
	M.setGroupAtPlan(xinf,0,eps,"Neumann");
	int spaceDim = M.getSpaceDimension();


	//Solid parameters
	double cp_ur=300;//Uranium specific heat
	double rho_ur=10000;//Uranium density
	double lambda_ur=5;
 
    bool FEcalculation=true;
	DiffusionEquation  myProblem(spaceDim,FEcalculation,rho_ur,cp_ur,lambda_ur);
	Field VV("Solid temperature", NODES, M, 1);

	//Set fluid temperature (temperature du fluide)
	double fluidTemp=573;//fluid mean temperature
	double heatTransfertCoeff=1000;//fluid/solid exchange coefficient
	myProblem.setFluidTemperature(fluidTemp);
	myProblem.setHeatTransfertCoeff(heatTransfertCoeff);
	//Set heat source
	Field Phi("Heat power field", NODES, M, 1);
	power_field_diffusionTest(Phi);
	myProblem.setHeatPowerField(Phi);
	Phi.writeVTK("1DheatPowerField");

	//Initial field creation
	Vector VV_Constant(1);
	VV_Constant(0) = 623;//Rod clad temperature

	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setNeumannBoundaryCondition("Neumann");

	// set the numerical method
	myProblem.setTimeScheme( Explicit);

	// name result file
	string fileName = "1DRodTemperature_FE";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 1000000;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

	// set display option to monitor the calculation
	myProblem.setVerbose( true);
	//set file saving format
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
