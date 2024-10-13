#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building Cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=1.0;
	double discontinuity = (xinf+xsup)/2.;
	int nx=4;
	Mesh M(xinf,xsup,nx);
	int spaceDim = M.getSpaceDimension();

	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(Gas, around1bar300K, spaceDim );

	double initialVelocity_Left=1;
	double initialPressure_Left=155e5;

	double initialVelocity_Right=1;
	double initialPressure_Right=150e5;


    // Prepare for the initial condition
	Vector Pressure_Left(1);
	Vector Pressure_Right(1);
	Vector Velocity_Left(1);
	Vector Velocity_Right(1);
	
	// left and right constant vectors		
	Pressure_Left[0] = initialPressure_Left;
	Pressure_Right[0] = initialPressure_Right;
	Velocity_Left[0] = initialVelocity_Left;
	Velocity_Right[0] = initialVelocity_Right;

    //Initial field creation
	cout << "Building initial data " <<endl; 
	int direction =0; // TODO : what is it ?
	myProblem.setInitialFieldStepFunction(M,Pressure_Left,Pressure_Right,discontinuity, direction, CELLS);
	myProblem.setInitialFieldStepFunction(M,Velocity_Left,Velocity_Right,discontinuity, direction, FACES);
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		bool isBoundary = Fj.isBorder();
		if (Fj.getNumberOfCells()==1){
			if (Fj.x() < discontinuity) {
				wallPressureMap[j] = initialPressure_Left ;
				wallVelocityMap[j] = initialVelocity_Left ;
			}
			else{
				wallPressureMap[j] = initialPressure_Right;
				wallVelocityMap[j] = initialVelocity_Right ;
			}
		}
	}

	
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit);
    
    // name of result file
	string fileName = "EulerBarotropicStaggered_1DRiemannProblem";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 4;
	int freqSave = 1;
	double cfl = 0.2;
	double maxTime = 30;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(CSV);
	myProblem.setVerbose(false); //TODO _A n'est pas initalisée or le code demande tout de meme à l'afficher en mode verbose
	
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
