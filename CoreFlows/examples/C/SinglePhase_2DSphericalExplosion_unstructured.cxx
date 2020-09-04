#include "SinglePhase.hxx"

using namespace std;

//Function that generates and save a spherical initial data
/*void initialField(){
	int spaceDim = 2;
	int nVar=2+spaceDim;
	double x,y;
	cout << "Loading unstructured mesh " << endl;
	Mesh M("../examples/resources/BoxWithMeshWithTriangularCells.med");
	Field VV("Primitive variables for spherical explosion", CELLS, M, nVar);
	vector<double>Vout(nVar), Vin(nVar);
	Vin[0]=1.1;
	Vin[1]=0;
	Vin[2]=0;
	Vin[3]=300;
	Vout[0]=1;
	Vout[1]=0;
	Vout[2]=0;
	Vout[3]=300;

	for(int i=0;i<M.getNumberOfCells();i++){
		x=M.getCell(i).x();
		y=M.getCell(i).y();
		if((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<0.25*0.25)
			for(int j=0;j<nVar;j++)
				VV(i,j)=Vin[j];
		else
			for(int j=0;j<nVar;j++)
				VV(i,j)=Vout[j];
	}
	//VV.writeMED("../examples/ressources/BoxWithMeshWithTriangularCells",false);
}*/
int main(int argc, char** argv)
{
	// preprocessing: mesh and group creation
	cout << "Loading unstructured mesh and initial data for test SinglePhase_2DSphericalExplosion_unstructured()" << endl;
	string inputfile="resources/BoxWithMeshWithTriangularCells";
	string fieldName="Initial variables for spherical explosion";
	int spaceDim=2;

	SinglePhase  myProblem(Gas,around1bar300K,spaceDim);

	//Initial field creation
	cout << "Loading unstructured mesh and initial data for test SinglePhase_2DSphericalExplosion_unstructured()" << endl;
	myProblem.setInitialField(inputfile,fieldName,0);

	//set the boundary conditions
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallTemperature=563;
	myProblem.setWallBoundaryCondition("GAUCHE", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("DROITE", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("HAUT", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("BAS", wallTemperature, wallVelocityX, wallVelocityY);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);

	// name file save
	string fileName = "2DSphericalExplosion_unstructured";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3 ;
	int freqSave = 5;
	double cfl = 0.5;
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
