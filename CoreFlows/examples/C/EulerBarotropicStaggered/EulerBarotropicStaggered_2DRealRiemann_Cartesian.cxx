#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <cstdlib>

using namespace std;

double initialPressure( double x, double y, double discontinuity_x, double discontinuity_y){
	if (x < discontinuity_x){
		if (y < discontinuity_y){
			return 1;
		}
		else{
			return 4;
		}
	}
	else {
		if (y < discontinuity_y){
			return 2;
		}
		else{
			return 0.5;
		}
	}
}

std::vector<double> initialVelocity(double x, double y, double discontinuity_x, double discontinuity_y){
	std::vector<double> vec(2, 0.0); 
	if (x < discontinuity_x){
		if (y < discontinuity_y){
			vec[0] = 3;
			vec[1] = 0;
		}
		else{
			vec[0] = 0;
			vec[1] = -1;
		}
	}
	else {
		if (y < discontinuity_y){
			vec[0] = -2;
			vec[1] = 0;
		}
		else{
			vec[0] = 3;
			vec[1] = 0;
		}
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
	
	//Preprocessing: mesh and group creation
	int spaceDim = 2;
	PetscInitialize(&argc, &argv, 0,0);
	// Prepare for the mesh
	cout << "Building mesh" << endl;
	cout << "Construction of a cartesian mesh" << endl;
	double inf = 0.0;
	double sup = 1.0;
	double discontinuity_x, discontinuity_y;
	int nx, ny;
	nx=25;
	ny=25;
	discontinuity_x = (inf + sup)/2.0 +  0.75/nx;
	discontinuity_y = (inf + sup)/2.0 +  0.75/ny;

	Mesh M=Mesh(inf,sup,nx,inf,sup,ny);
	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );

	// Prepare for the initial condition
	// set the boundary conditions
	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	Field Pressure0("pressure", CELLS, M, 1);
	Field Velocity0("velocity", FACES, M, 1);

	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		std::vector<int> idCells = Fj.getCellsId();
		std::vector<double> vec_normal_sigma(2) ; 
		Cell Ctemp1 = M.getCell(idCells[0]);
		for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
			if (j == Ctemp1.getFacesId()[l]){
				for (int idim = 0; idim < spaceDim; ++idim)
					vec_normal_sigma[idim] = fabs( Ctemp1.getNormalVector(l,idim) ) ;
			}
		}
		myProblem.setOrientation(j,vec_normal_sigma);
		if(Fj.getNumberOfCells()==2 ){ 
			myProblem.setInteriorIndex(j);
			Cell Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(), Ctemp1.y(),discontinuity_x , discontinuity_y);
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(), Ctemp2.y(),discontinuity_x , discontinuity_y);
			Velocity0[j] = dotprod(initialVelocity(Fj.x(), Fj.y(), discontinuity_x, discontinuity_y),vec_normal_sigma );
		}
		else if (Fj.getNumberOfCells()==1  ){ 
			myProblem.setWallBoundIndex(j);	
			wallVelocityMap[j] = 0;
		}
	}
	
	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

	// set the numerical method
	myProblem.setTimeScheme(Explicit);
	
	// name of result file
	string fileName = "EulerBarotropicStaggered_2DRiemann_StructuredSquares";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 200000000;
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
	myProblem.terminate();
		
	return EXIT_SUCCESS;
}
