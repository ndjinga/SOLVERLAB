#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 3;
	/*Preprocessing: mesh data*/
	double xinf=0;
	double xsup=1;
	double yinf=0;
	double ysup=1;
	double zinf=0;
	double zsup=1;
	int nx=10;
	int ny=10;
	int nz=10;

	/* set the limit field for each boundary*/
	double coldWallVelocityX=0;
	double coldWallVelocityY=0;
	double coldWallVelocityZ=0;
	double coldWallTemperature=563;

	double hotWallVelocityX=0;
	double hotWallVelocityY=0;
	double hotWallVelocityZ=0;
	double hotWallTemperature=613;

	/* physical constants*/
	vector<double> gravite(spaceDim,0.) ;
	gravite[2]=-10;
	gravite[1]=0;
	gravite[0]=0;
	vector<double> viscosite(1), conductivite(1);
	viscosite[0]= 8.85e-5;
	conductivite[0]=1000;//nucleate boiling heat transfert coefficient

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	//Initial field creation
	cout << "Construction de la condition initiale" << endl;
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = 155e5;
	VV_Constant[1] = 0;
	VV_Constant[2] = 0;
	VV_Constant[3] = 0;
	VV_Constant[4] = 573;
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,xinf,xsup,nx,"hotWall","hotWall",
															yinf,ysup,ny,"hotWall","hotWall",
															zinf,zsup,nz,"hotWall","coldWall");

	//set the boundary conditions
	myProblem.setWallBoundaryCondition("coldWall", coldWallTemperature, coldWallVelocityX, coldWallVelocityY, coldWallVelocityZ);
	myProblem.setWallBoundaryCondition("hotWall", hotWallTemperature, hotWallVelocityX, hotWallVelocityY, hotWallVelocityZ);


	// physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// set the Petsc resolution
	myProblem.setLinearSolver(GMRES,ILU,false);

	// name result file
	string fileName = "3DHeatDrivenCavity";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 10;
	double maxTime = 50;
	double precision = 1e-6;

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
