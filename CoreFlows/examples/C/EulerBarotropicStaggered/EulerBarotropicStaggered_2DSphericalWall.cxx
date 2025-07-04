#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <cstdlib>
#include <chrono>

using namespace std;

double initialDensity( double x, double y){
	return 2;
}

std::vector<double> initialVelocity(double x, double y){
	std::vector<double> vec(2);
	double theta = atan2(x,y); 
	if (theta < 0) theta += 2 * M_PI;
	double vr = 1;
	vec[0] =    cos(theta) * 0 + sin(theta) * vr  ;
	vec[1] =  - sin(theta) * 0 + cos(theta) * vr  ;
	
	/* vec[0] =   	 -sin(theta) * 0 + cos(theta) * vr   ;
	vec[1] =  	  cos(theta) * 0 + sin(theta) * vr  ; */
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
	char wall = *(argv[1]);
	int spaceDim = 2;
	Mesh M;
	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );
	int nx;
	if (argc<1 || (*(argv[1]) != 'l' && *(argv[1]) != 'r')  ){
		cout << "ERROR : you have to give the side of the boundary condition, either 'l' for left side or 'r' for the right side";
	}
	if (argc>1  ){
		// ./resources/AnnulusSpiderWeb5x16.med or ./resources/AnnulusSpiderWeb60x20.med
		cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
		cout << "Loading of a mesh named "<<argv[1] << endl;
		string filename = "./resources/AnnulusSpiderWeb60x20.med"; //argv[1];
		nx=40;
		M=Mesh(filename);
	}

	// Prepare for the initial condition
	// set the boundary conditions
	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallDensityMap;
	std::map<int ,double> wallVelocityMap ;
	Field Density0("Density", CELLS, M, 1);
	Field Momentum0("velocity", FACES, M, 1);
	std::vector<double> wallVelocityVector(spaceDim);

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
		//TODO pb orientation ?
		if (  Fj.x() >1e-10 && fabs( atan(Fj.y()/Fj.x()) ) <1e-10 ){ 
			vec_normal_sigma[0] *= -1;
			vec_normal_sigma[1] *= -1;
		} 
		myProblem.setOrientation(j,vec_normal_sigma);

		Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y());
		if(Fj.getNumberOfCells()==2 ){ 
			myProblem.setInteriorIndex(j);
			Cell Ctemp2 = M.getCell(idCells[1]);
			Density0[idCells[1]] = initialDensity(Ctemp2.x(),Ctemp2.y());
			Momentum0[j] = dotprod(initialVelocity(Fj.x(), Fj.y() ) ,vec_normal_sigma )*(Density0[idCells[0]]+Density0[idCells[1]])/2.0;
		}
		else if (Fj.getNumberOfCells()==1  ){ 
			Momentum0[j] = dotprod(initialVelocity(Fj.x(), Fj.y()),vec_normal_sigma )*Density0[idCells[0]];
			double r =sqrt(pow(Fj.x(),2) + pow(Fj.y(),2));
			if (wall =='l' && r<(0.8+ 6.0)/2.0 ){
				myProblem.setWallBoundIndex(j);
				wallVelocityMap[j] = 0;
				for (int idm = 0 ;idm <spaceDim; idm ++)
					wallVelocityVector[idm] = 0;
				myProblem.setboundaryVelocityVector(j, wallVelocityVector);
			}
			else if (wall =='r' && r>(0.8+ 6.0)/2.0 ){
				myProblem.setWallBoundIndex(j);
				wallVelocityMap[j] = 0;
				for (int idm = 0 ;idm <spaceDim; idm ++)
					wallVelocityVector[idm] = 0;
				myProblem.setboundaryVelocityVector(j, wallVelocityVector);
			}
			else myProblem.setNeumannBoundIndex(j);	
			
		}
	}
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryVelocity(wallVelocityMap);

	// set the numerical method
	myProblem.setTimeScheme(Explicit);
	
	// name of result file
	string fileName = "EulerBarotropicStaggered_2DSphericalWall";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 1000000;
	int freqSave = 40;
	double cfl = 0.99;
	double maxTime = 0.1;
	double precision = 1e-10;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(VTK);
	myProblem.saveVelocity(true);
	myProblem.saveSphericalVelocity(true);
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
