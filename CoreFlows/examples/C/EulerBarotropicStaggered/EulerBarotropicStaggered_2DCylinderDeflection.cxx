#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <chrono>

#define M_ref 1e-4
#define rho_b 2
#define gamma 2
#define kappa 1
#define c_b sqrt(gamma * kappa* pow(rho_b, gamma-1) )
#define u0 M_ref * c_b
using namespace std;

std::vector<double> ExactVelocity(double r, double theta, double r1, double r0){
	std::vector<double> vec(2);
	vec[0] = u0 * pow(r1,2)/(pow(r1,2) - pow(r0,2))* (1 - pow(r0,2)/pow(r,2) * cos(2*theta)); 
	vec[1] = u0 * pow(r1,2)/(pow(r1,2) - pow(r0,2))* (  - pow(r0,2)/pow(r,2) * sin(2*theta));  
	return vec;
}

double initialDensity( double x, double y){ return rho_b; }


std::vector<double> initialVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = 0; 
	vec[1] = 0;
	return vec;
}

std::vector<double> initialBoundVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = u0 ; 
	vec[1] = 0 ;
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
	std::map<int ,double> wallVelocityMap ;
	std::vector<double> wallVelocityVector(spaceDim);
	Field Density0("Density", CELLS, M, 1);
	Field Momentum0("velocity", FACES, M, 1);
	std::vector<double> ExactVelocityAtFaces(M.getNumberOfFaces());
	
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
		//TODO at theta=0; changing the sign of the basis function seems to give a better metric
		if (  Fj.x() >1e-10 && fabs( atan(Fj.y()/Fj.x()) ) <1e-10 )  vec_normal_sigma[1] *= -1;

		myProblem.setOrientation(j,vec_normal_sigma);
		double r =  sqrt(Fj.x()*Fj.x() + Fj.y()*Fj.y());
		double theta = atan2(Fj.y(),Fj.x()); 
		if (theta < 0) theta += 2 * M_PI;
		ExactVelocityAtFaces[j] = dotprod( ExactVelocity(r, theta, r1, r0), vec_normal_sigma); 

		Density0[idCells[0]] = initialDensity(Ctemp1.x(),Ctemp1.y());
		if(Fj.getNumberOfCells()==2){
			myProblem.setInteriorIndex(j);
			Density0[idCells[1]] = initialDensity(M.getCell(idCells[1]).x(),M.getCell(idCells[1]).y());
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * ( Density0[idCells[0]] + Density0[idCells[1]]  )/2;
		}
		else if (Fj.getNumberOfCells()==1){
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * initialDensity(Fj.x(),Fj.y());
			if (  sqrt( pow(Fj.x(),2) + pow(Fj.y(),2) )  <= (r0 +r1)/2.0 ){// if face is on interior (wallbound condition) r_int = 1.2 ou 0.8 selon le maillage
				myProblem.setWallBoundIndex(j);
				for (int idm = 0; idm <spaceDim; idm ++)	wallVelocityVector[idm] = 0;
			}
			else {		
				myProblem.setSteggerBoundIndex(j);								
				wallDensityMap[j] = initialDensity(Fj.x(),Fj.y());
				for (int idm = 0; idm <spaceDim; idm ++)    wallVelocityVector[idm] = initialBoundVelocity(Fj.x(), Fj.y())[idm];
			} 
			wallVelocityMap[j] = dotprod( wallVelocityVector, vec_normal_sigma );
			myProblem.setboundaryVelocityVector(j, wallVelocityVector);
		}
	}

	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit);
	myProblem.setLinearSolver(GMRES, LU, 50); //If Implicit
	double cfl = 1;
    
    // name of result file
	string fileName = "EulerBarotropicStaggered_2DCylinderDeflection";
    // parameters calculation
	unsigned MaxNbOfTimeStep = 10000000	;
	double precision = 1e-11;
	int freqSave = 1000;
	double maxTime = 50;
	
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(VTK);
	myProblem.saveVelocity(true);
	myProblem.savePressure(true);
	myProblem.setVerbose(false);
	
	// evolution
	myProblem.initialize();
	myProblem.InterpolateFromFacesToCells(ExactVelocityAtFaces);
	bool ok = myProblem.run();
	myProblem.computeOrder2Density(rho_b, M_ref);

	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;

	cout << "For "<< M.getNumberOfCells() <<" cells, Velocity error  AT FACES = "<<myProblem.ErrorVelocity(ExactVelocityAtFaces)[0]<< " AT CELLS X = "<<myProblem.ErrorVelocity(ExactVelocityAtFaces)[1]<<" AT CELLS Y = "<<myProblem.ErrorVelocity(ExactVelocityAtFaces)[2]<<endl;
	auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    cout << "Execution time: " << duration.count() << " ms" << " or "<< duration.count()/60000.0<<" minutes "<< endl;

	myProblem.terminate();
	

	return EXIT_SUCCESS;
}
