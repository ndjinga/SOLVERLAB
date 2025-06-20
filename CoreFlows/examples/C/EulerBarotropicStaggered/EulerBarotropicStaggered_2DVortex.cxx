#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <cstdlib>
#include <chrono>

# define M_ref 1e-2
# define rho0  2
# define u0  M_ref * sqrt(2*rho0)
# define alpha 2
# define r0 	0.45
# define rbarmax      sqrt(- pow(alpha,2) + sqrt(1 + pow(alpha,4)) )
# define xc 0.5
# define yc 0.5
# define lambda_max     2 * alpha *rbarmax/(1-pow(rbarmax,2)) * exp(- pow(alpha, 2)/(1 - pow(rbarmax,2) ) ) 

using namespace std;


//******************** STATIONNARY VORTEX *************//
double rhoVortex( double x, double y ){
	double rbarmax2 = - pow(alpha,2) + sqrt(1 + pow(alpha,4));
	double r = sqrt(pow(x- 0.5,2) +pow(y- 0.5,2) );
	double rbar = r/r0; 
	double expr = exp( - pow(alpha,2)/(1-pow(rbar,2)) );
	double umax =  2 * alpha * rbarmax / (1-rbarmax2) * exp( - alpha*alpha/(1-rbarmax2) ) ;
	return  rho0 * (1 - pow(M_ref * expr/umax,2) ) ;

}
std::vector<double> MomentumVortex(double x, double y ){
	std::vector<double> vec(2,0.0);
	double rbarmax2 = - pow(alpha,2) + sqrt(1 + pow(alpha,4));
	double r = sqrt(pow(x- 0.5,2) +pow(y- 0.5,2) );
	double rbar = r/r0; 
	double expr = exp( - alpha*alpha/(1-rbar*rbar) );
	double umax =  2 * alpha * rbarmax / (1-rbarmax2) * exp( - alpha*alpha/(1-rbarmax2) ) ;
	double rho = rho0 * (1 - pow(M_ref * expr/umax,2) ) ;
	vec[0] =  rho * (y-0.5)  * u0 * 2 * alpha/r0 * expr / ( (1-rbar*rbar) * umax) ;
	vec[1] = -rho * (x-0.5)  * u0 * 2 * alpha/r0 * expr / ( (1-rbar*rbar) * umax); 
	return vec;
}

//******************** LOW-MACH WAVE *****************//
double rhoLowMach( double xbar){  return rho0* (1 + M_ref * exp(1.0 - 1.0/(1-pow(xbar,2))) ) ; }

std::vector<double> MomentumLowMach(double xbar ){
	std::vector<double> vec(2,0.0); 
	double rhoLow = rho0* (1 + M_ref * exp(1.0 - 1.0/(1-pow(xbar,2))) );
	double du = 2 *( sqrt( 2*rhoLow ) - sqrt(2*rho0) );
	double drho = rho0 *M_ref * exp(1.0 - 1.0/(1.0-pow(xbar,2))) ;
	vec[0] = rho0 * du + drho * u0 + drho *du;
	return vec;
}

//******************** INITIAL CONDITION*****************//
double InitialDensity( double x, double y ){
	//if ( sqrt(pow(x- 0.5,2) +pow(y- 0.5,2) )<r0 ) return rhoVortex( x, y); 
	if ( fabs(x/0.05)<1 )  						  return rhoLowMach(x/0.05); 
	else                       					  return rho0;
}
std::vector<double> InitialMomentum(double x, double y ){
	std::vector<double> vec(2,0.0); 
	//if ( sqrt(pow(x- 0.5,2) +pow(y- 0.5,2) )<r0 ) 	return MomentumVortex( x, y); 
	if ( fabs(x/0.05)<1 )   						return MomentumLowMach( x/0.05); 
	else 											return vec;
}

double dotprod(std::vector<double> vector, std::vector<double> normal){
	assert(vector.size() == normal.size());
	double dotprod =0;
	for (int n =0; n< vector.size(); n++){
		dotprod += vector[n] * normal[n];
	}
	return dotprod;
}




int main(int argc, char** argv){	
	auto start = std::chrono::high_resolution_clock::now();
	
	//Preprocessing: mesh and group creation
	int spaceDim = 2;
	PetscInitialize(&argc, &argv, 0,0);
	// Prepare for the mesh
	cout << "Building mesh" << endl;
	cout << "Construction of a cartesian mesh" << endl;

	double infx = -0.1;
	double supx = 1.1;
	double infy = 0.0;
	double supy = 1.0;
	int nx =100;
	int ny =2;
	Mesh M=Mesh(infx, supx, nx, infy, supy, ny);
	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );

	// Prepare for the initial condition
	// set the boundary conditions
	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallDensityMap;
	std::map<int ,double> wallVelocityMap ;
	Field Density0("Density", CELLS, M, 1);
	Field Momentum0("velocity", FACES, M, 1);
	std::vector<double> wallVelocityVector(spaceDim);


	myProblem.setPeriodicFaces(M, 'x', nx, infx, supx); 

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
		Density0[idCells[0]] = InitialDensity(Ctemp1.x(),Ctemp1.y());
		if(Fj.getNumberOfCells()==2 ){ 
			Cell Ctemp2 = M.getCell(idCells[1]);
			myProblem.setInteriorIndex(j);
			Density0[idCells[1]] = InitialDensity(Ctemp2.x(),Ctemp2.y());
			Momentum0[j] = dotprod(InitialMomentum(Fj.x(),Fj.y()),vec_normal_sigma);
		}
		else if (Fj.getNumberOfCells()==1  ){ 
			Momentum0[j] = dotprod(InitialMomentum(Fj.x(),Fj.y()),vec_normal_sigma);
			// if periodic check that the boundary face is the computed (avoid passing twice ) 
			if  (myProblem.IsFaceBoundaryNotComputedInPeriodic(j) == false && myProblem.IsFaceBoundaryComputedInPeriodic(j) == false)
				myProblem.setSteggerBoundIndex(j);	
			// Boundary normal velocity, Density and full velocity vector
			for (int idm = 0; idm <spaceDim; idm ++)
				wallVelocityVector[idm] = InitialMomentum(Fj.x(),Fj.y())[idm]/Density0[idCells[0]];
			myProblem.setboundaryVelocityVector(j, wallVelocityVector);
			wallVelocityMap[j] = dotprod(wallVelocityVector,vec_normal_sigma );
			wallDensityMap[j]  = Density0[idCells[0]];
		}
	}
	myProblem.setInitialField(Density0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallDensityMap);
	myProblem.setboundaryVelocity(wallVelocityMap);
	
	// set the numerical method
	myProblem.setTimeScheme(Explicit);
	myProblem.setLinearSolver(GMRES, LU, 70);
	
	// name of result file
	string fileName = "EulerBarotropicStaggered_2DVortex";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 1000000;
	int freqSave = 20;
	double cfl = 1;
	double maxTime = 3.5;
	double precision = 1e-9;

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
