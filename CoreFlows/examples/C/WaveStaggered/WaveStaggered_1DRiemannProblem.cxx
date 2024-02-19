#include "WaveStaggered.hxx"
#include "math.h"

using namespace std;

// set the boundary conditions
	double boundPressure(double x)
	{
	 return 155e7;
	}

	double boundVelocity(double x)
	{
		return 1;
	}

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building Cartesian mesh " << endl;
	double xinf=0.0;
	double xsup=1.0;
	int nx=2;
	Mesh M(xinf,xsup,nx);
	int spaceDim = M.getSpaceDimension();

	double kappa =2;
	double rho = 5;
	WaveStaggered myProblem = WaveStaggered(spaceDim, kappa, rho );

    // Prepare for the initial condition
	std::vector<double> initialVelocity;
	initialVelocity.push_back(1);
	std::vector<double> initialPressure;
	initialPressure.push_back(155e7);

	//Initial field creation
	cout << "Building initial data " << endl; 
	myProblem.setInitialFieldConstant(M, initialVelocity, FACES);
	myProblem.setInitialFieldConstant(M, initialPressure, CELLS);

	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		bool isBoundary = Fj.isBorder();
		if (isBoundary == true){
			wallPressureMap[j] = boundPressure(Fj.x()) ;
			wallVelocityMap[j] = boundVelocity(Fj.x()) ;
		}
	}

	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit);
    
    // name of result file
	string fileName = "WaveStaggered_1DRiemannProblem";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 0.2;
	double maxTime = 5;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(VTK);
	
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
