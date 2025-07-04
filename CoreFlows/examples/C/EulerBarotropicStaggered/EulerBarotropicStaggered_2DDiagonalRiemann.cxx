#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <cstdlib>
#include <chrono>

using namespace std;

double initialDensity( double x, double y , double nx){
	if (y < (x+0.1 - 1.0/(4 * nx) ) )	return 3.0 ;//12;
	else 									return 1.0;
}

std::vector<double> initialVelocity(double x, double y, double nx){
	std::vector<double> vec(2);
	double ul = 5;//1.5 ; 
	double ur = -1.0;//-3 ; 
	if (y < (x+0.1  -1.0/(4 * nx) )){
		vec[0] = -ul/sqrt(2.0);
		vec[1] = ul/sqrt(2.0);
	}
	else {
		vec[0] = -ur/sqrt(2.0);
		vec[1] = ur/sqrt(2.0);
	}
	
	return vec;
}

double dotprod(std::vector<double> vector, std::vector<double> normal){
	assert(vector.size() == normal.size());
	double dotprod =0;
	for (int n =0; n< vector.size(); n++){
		dotprod += vector[n] * normal[n];
	}
	return dotprod;
}

int main(int argc, char** argv)
{	
	auto start = std::chrono::high_resolution_clock::now();

	int spaceDim = 2;
	Mesh M;
	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );
	int nx;
	double inf = 0.0;
	double sup = 1.0;
	if (argc>1  ){
		// ./resources/AnnulusSpiderWeb5x16.med or ./resources/AnnulusTriangles60.med
		cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
		cout << "Loading of a mesh named "<<argv[1] << endl;
		string filename = argv[1];
		nx=55;
		M=Mesh(filename);
	}
	else{
		PetscInitialize(&argc, &argv, 0,0);
		// Prepare for the mesh
		cout << "Building mesh" << endl;
		cout << "Construction of a cartesian mesh" << endl;
	
		int ny;
		nx=10	;
		ny=10;
		M=Mesh(inf,sup,nx,inf,sup,ny);
	}
	//Preprocessing: mesh and group creation

	// Prepare for the initial condition
	// set the boundary conditions
	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallDensityMap;
	std::map<int ,double> wallVelocityMap ;
	Field Density0("Density", CELLS, M, 1);
	Field Momentum0("velocity", FACES, M, 1);
	std::vector<double> wallVelocityVector(spaceDim);

	//myProblem.setPeriodicFaces(M, 'd', nx, inf , sup);
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		std::vector<int> idCells = Fj.getCellsId();
		std::vector<double> vec_normal_sigma(2) ; 
		Cell Ctemp1 = M.getCell(idCells[0]);
		for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
			if (j == Ctemp1.getFacesId()[l]){
				for (int idim = 0; idim < spaceDim; ++idim)
					vec_normal_sigma[idim] =  Ctemp1.getNormalVector(l,idim)  ;
			}
		}
		myProblem.setOrientation(j,vec_normal_sigma);

		Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y(), nx);
		if(Fj.getNumberOfCells()==2 ){ 
			myProblem.setInteriorIndex(j);
			Cell Ctemp2 = M.getCell(idCells[1]);
			Density0[idCells[1]] = initialDensity(Ctemp2.x(),Ctemp2.y(), nx);
			Momentum0[j] = dotprod(initialVelocity(Fj.x(), Fj.y(),nx ) ,vec_normal_sigma )*(Density0[idCells[0]]+Density0[idCells[1]])/2.0;
		}
		else if (Fj.getNumberOfCells()==1  ){ 
			Momentum0[j] = dotprod(initialVelocity(Fj.x(), Fj.y(),nx),vec_normal_sigma )*Density0[idCells[0]];
			myProblem.setNeumannBoundIndex(j);	
		}
	}
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);

	// set the numerical method
	myProblem.setTimeScheme(Explicit);
	
	// name of result file
	string fileName = "EulerBarotropicStaggered_2DDiagonalRiemann";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 1000000;
	int freqSave = 1;
	double cfl = 0.99;
	double maxTime = 0.07;
	double precision = 1e-10;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(VTK);
	myProblem.saveVelocity(true);
	myProblem.savePressure(true);
	myProblem.setVerbose(false);
	
	// evolution
	myProblem.initialize();
	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> duration = end - start;
	std::cout << "Execution time: " << duration.count() << " ms" << std::endl;

	myProblem.terminate();
		
	return EXIT_SUCCESS;
}
