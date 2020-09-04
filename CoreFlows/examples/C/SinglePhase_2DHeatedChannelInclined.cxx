#include "SinglePhase.hxx"
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

    // Prepare for the mesh
	double xinf = 0 ;
	double xsup=3.0;
	double yinf=0.0;
	double ysup=5.0;
	int nx=10;
	int ny=10; 

    // set the limit field for each boundary
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallTemperature=573;
	double inletVelocityX=0;
	double inletVelocityY=0.5;
	double inletTemperature=563;
	double outletPressure=155e5;

    // physical constants
	vector<double> gravite (spaceDim);
    
	gravite[1]=-7;
	gravite[0]=7;

	double 	heatPower=1e8;

	SinglePhase myProblem(Liquid,around155bars600K,spaceDim);
	int nVar =myProblem.getNumberOfVariables();

    // Prepare for the initial condition
	vector<double> VV_Constant (nVar);

	// constant vector
	VV_Constant[0] = outletPressure ;
	VV_Constant[1] = inletVelocityX;
	VV_Constant[2] = inletVelocityY;
	VV_Constant[3] = inletTemperature ;

    //Initial field creation
	cout<<"Building initial data"<<endl;
	myProblem.setInitialFieldConstant(spaceDim,VV_Constant,
                                          xinf,xsup,nx,"wall","wall",
					  yinf,ysup,ny,"inlet","outlet", 
					  0.0,0.0,  0,  "", "");

    // the boundary conditions
	vector<double>pressure_reference_point(2);
	pressure_reference_point[0]=xsup;
	pressure_reference_point[1]=ysup;
	myProblem.setOutletBoundaryCondition("outlet", outletPressure,pressure_reference_point);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletVelocityX, inletVelocityY);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY);
    
	// set physical parameters
	myProblem.setHeatSource(heatPower);
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(staggered, Implicit);
	myProblem.setNonLinearFormulation(VFFC);
    
	// name file save
	string fileName = "2DInclinedHeatedChannel";

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
	myProblem.usePrimitiveVarsInNewton(true);
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

