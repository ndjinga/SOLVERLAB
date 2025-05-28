#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <chrono>

#define M_ref 1e-3
#define rho_b 2
#define gamma 2
#define kappa 1
using namespace std;

std::vector<double> ExactVelocity(double r, double theta, double r1, double r0){
	std::vector<double> vec(2);
	vec[0] = sqrt(gamma * kappa* pow(rho_b, gamma-1) ) * M_ref *r1*r1/(r1*r1 -r0*r0)* (1 - r0*r0/(r*r) * cos(2*theta)); 
	vec[1] = sqrt(gamma * kappa* pow(rho_b, gamma-1) ) * M_ref *r1*r1/(r1*r1 -r0*r0)* (  - r0*r0/(r*r) * sin(2*theta));  
	return vec;
}

double initialDensity( double x, double y){
	return rho_b;
}

// sqrt(p'(rho_0)) M_\infty
std::vector<double> initialVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = 0; 
	vec[1] = 0;
	return vec;
}

std::vector<double> initialBoundVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = sqrt(gamma * kappa* pow(rho_b, gamma-1) ) *M_ref ; 
	vec[1] = 0;
	return vec;
}

double dotprod(std::vector<double> vector1, std::vector<double> vector2){
	assert(vector1.size() == vector2.size());
	double dotprod =0;
	for (int n =0; n< vector1.size(); n++)
		dotprod += vector1[n] * vector2[n];
	return dotprod;
}

int main(int argc, char** argv)
{
    auto start = std::chrono::high_resolution_clock::now();
   
	//Preprocessing: mesh and group creation
	int spaceDim = 2;
	
	// Prepare for the mesh
	cout << "Building mesh" << endl;
	double r0 = 0.8;
	double r1 = 6;

	Mesh M;
	// ./resources/AnnulusSpiderWeb5x16.med or ./resources/AnnulusTriangles60.med
	cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
	cout << "Loading of a mesh named "<<argv[1] << endl;
	string filename = argv[1];
	M=Mesh(filename);

	double a = 1.0;
	//gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );

	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallDensityMap;
	std::map<int ,double> wallMomentumMap ;
	std::vector<double> wallVelocityVector(spaceDim);
	Field Density0("Density", CELLS, M, 1);
	Field Momentum0("velocity", FACES, M, 1);
	Field ExactVelocityAtFaces("ExactVelocityAtFaces", FACES, M, 1);
	
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		std::vector<int> idCells = Fj.getCellsId();
		std::vector<double> vec_normal_sigma(2) ; 
		Cell Ctemp1 = M.getCell(idCells[0]);
		for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
			if (j == Ctemp1.getFacesId()[l]){
				for (int idim = 0; idim < spaceDim; ++idim)
					vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
			}
		}
		if ( Fj.x() >1e-10 && fabs( atan(Fj.y()/Fj.x()) ) <1e-10 ){
			vec_normal_sigma[0] *= -1;
			vec_normal_sigma[1] *= -1;

		}  
		myProblem.setOrientation(j,vec_normal_sigma);
		if(Fj.getNumberOfCells()==2){
			Cell Ctemp2 = M.getCell(idCells[1]);
			myProblem.setInteriorIndex(j);
			Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y());
			Density0[idCells[1]] = initialDensity(Ctemp2.x(),Ctemp2.y());
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * ( initialDensity(Ctemp1.x(),Ctemp1.y()) + initialDensity(Ctemp2.x(),Ctemp2.y())  )/2;
		}
		else if (Fj.getNumberOfCells()==1){
			Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y());
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * initialDensity(Fj.x(),Fj.y());
			if (( sqrt( Fj.x()*Fj.x()+ Fj.y()*Fj.y() )  ) <= (r0 +r1)/2.0 ){// if face is on interior (wallbound condition) r_int = 1.2 ou 0.8 selon le maillage
				myProblem.setWallBoundIndex(j);
				wallMomentumMap[j] =  0;
				for (int idm = 0; idm <spaceDim; idm ++)
					wallVelocityVector[idm] = 0;
			}
			else {		
				myProblem.setSteggerBoundIndex(j);								
				wallMomentumMap[j] = dotprod( initialBoundVelocity( Fj.x(),Fj.y()), vec_normal_sigma );
				wallDensityMap[j] = initialDensity(Fj.x(),Fj.y());
				for (int idm = 0; idm <spaceDim; idm ++)
					wallVelocityVector[idm] = initialBoundVelocity(Fj.x(), Fj.y())[idm];
			} 
			ExactVelocityAtFaces[j] = wallMomentumMap[j];
			myProblem.setboundaryVelocityVector(j, wallVelocityVector);
		}
	}

		
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallMomentumMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit	);
    
    // name of result file
	string fileName = "EulerBarotropicStaggered_2DCylinderDeflection";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 100000000	;
	int freqSave = 300	;
	double cfl = 0.99;
	double maxTime = 50;
	double precision = 1e-8;

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
	Field ExactVelocityInterpolate("ExactVelocityInterpolate", CELLS, M, 3);
	//myProblem.InterpolateFromFacesToCells(ExactVelocityAtFaces, ExactVelocityInterpolate);
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
