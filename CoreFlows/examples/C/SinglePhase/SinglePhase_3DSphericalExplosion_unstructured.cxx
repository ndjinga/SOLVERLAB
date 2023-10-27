#include "SinglePhase.hxx"

using namespace std;


int main(int argc, char** argv)
{
	// preprocessing: mesh and group creation
	cout << "Loading unstructured mesh for test SinglePhase_3DSphericalExplosion_unstructured()" << endl;
	string inputfile="resources/meshCube.med";

	double xinf=0;
	double xsup=1;
	double yinf=0;
	double ysup=1;
	double zinf=0;
	double zsup=1;
	Mesh M(inputfile);
	double eps=1.E-6;
	M.setGroupAtPlan(xinf,0,eps,"GAUCHE");
	M.setGroupAtPlan(xsup,0,eps,"DROITE");
	M.setGroupAtPlan(yinf,1,eps,"ARRIERE");
	M.setGroupAtPlan(ysup,1,eps,"AVANT");
	M.setGroupAtPlan(zinf,2,eps,"BAS");
	M.setGroupAtPlan(zsup,2,eps,"HAUT");

	/* Initial field data */
	int spaceDim = 3;
	int nVar=2+spaceDim;
	double radius=0.5;
	Vector Center(3);//default value is (0,0,0)
	Vector Vout(nVar), Vin(nVar);
	Vin[0]=1.1;
	Vin[1]=0;
	Vin[2]=0;
	Vin[3]=0;
	Vin[4]=300;
	Vout[0]=1;
	Vout[1]=0;
	Vout[2]=0;
	Vout[3]=0;
	Vout[4]=300;


	SinglePhase  myProblem(Gas,around1bar300K,spaceDim);

	/*Setting mesh and Initial */
	cout << "Setting initial data " << endl;
	myProblem.setInitialFieldSphericalStepFunction( M, Vout, Vin, radius, Center);

	//set the boundary conditions
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallVelocityZ=0;
	double wallTemperature=563;
	myProblem.setWallBoundaryCondition("GAUCHE", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("DROITE", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("HAUT", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("BAS" , wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("AVANT", wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);
	myProblem.setWallBoundaryCondition("ARRIERE" , wallTemperature, wallVelocityX, wallVelocityY, wallVelocityZ);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);

	// name file save
	string fileName = "3DSphericalExplosion_unstructured";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 0.3;
	double maxTime = 5;
	double precision = 1e-6;

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
