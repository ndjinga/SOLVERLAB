#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

	//Prepare for the mesh
	double xinf=0.0;
	double xsup=.6;
	double yinf=0.0;
	double ysup=2.0;
	int nx=3;
	int ny=100;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny);

	//Set the barriers
	double xcloison1=xinf+(xsup-xinf)/3;
	double xcloison2=xinf+2*(xsup-xinf)/3;
	Field barrierField("Barrier Field", FACES, M, 1);
	double eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"wall");
	M.setGroupAtPlan(xinf,0,eps,"wall");
	M.setGroupAtPlan(ysup,1,eps,"outlet");
	M.setGroupAtPlan(yinf,1,eps,"inlet");
	double dy=(ysup-yinf)/ny;
	int ncloison=3*ny/4;
	int i=0;
	while( i<= ncloison+1)
	{
		M.setGroupAtFaceByCoords(xcloison1,yinf+((ysup-yinf)/4)+(i+0.5)*dy,0,eps,"wall");
		M.setGroupAtFaceByCoords(xcloison2,yinf+((ysup-yinf)/4)+(i+0.5)*dy,0,eps,"wall");
		i++;
	}

	int nbFaces=M.getNumberOfFaces();
	for( i=0;i<nbFaces;i++)
	{
		double x=M.getFace(i).x();
		double y=M.getFace(i).y();
		if (((y> yinf+(ysup-yinf)/4) && (abs(x-xcloison1)< eps or abs(x-xcloison2)< eps) ) || abs(x-xinf)< eps || abs(x-xsup)< eps)
			barrierField[i]=1;
		else
			barrierField[i]=0;
	}

	barrierField.writeVTK("barrierField",true);

	// set the limit field for each boundary
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallTemperature=563;

	double inletConcentration=0;
	double inletVelocityX=0;
	double inletVelocityY=1;
	double inletTemperature=563;

	double outletPressure=155e5;

	// physical constants
	vector<double> gravite(spaceDim,0.) ;
	gravite[1]=-7;
	gravite[0]=7;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	// Prepare for the initial condition
	vector<double> VV_Constant(nVar);
	// constant vector
	VV_Constant[0] = 0;
	VV_Constant[1] = 155e5;
	VV_Constant[2] = 0;
	VV_Constant[3] = 1;
	VV_Constant[4] = 563;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	vector<double>pressure_reference_point(2);
	pressure_reference_point[0]=xsup;
	pressure_reference_point[1]=ysup;
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,pressure_reference_point);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConcentration, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);

	// set physical parameters
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(true);
	myProblem.setNonLinearFormulation(VFFC);

	// name of result file
	string fileName = "2DInclinedChannelGravityBarriers";

	// computation parameters
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 500;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.usePrimitiveVarsInNewton(true);

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
