#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>
#include <cstdlib>

using namespace std;

double initialPressure( double z, double discontinuity){
	if (z < discontinuity)
		return 2;
	else
		return 2;
}

std::vector<double> initialVelocity(double z, double discontinuity, char Direction){
	std::vector<double> vec(2);
	double ul = 1;
	double ur = 1;
	if (Direction == 'x'){
		if (z < discontinuity){
			vec[0] = ul;
			vec[1] = 0;
		}
		else {
			vec[0] = ur;
			vec[1] = 0;
		}
	}
	else if (Direction == 'y'){
		if (z < discontinuity){
			vec[0] = 0;
			vec[1] = ul;
		}
		else {
			vec[0] = 0;
			vec[1] = ur;
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
	if (argc<3 ||  ( (*(argv[1]) != 'x' && *(argv[1]) != 'y') && (*(argv[2]) != 'l' && *(argv[2]) != 'r') ) ){
		cout << "ERROR : you have to give a direction for the pseudo 1d Riemann problem, either 'x' or 'y' AND the side of the boundary condition, either 'l' for left side or 'r' for the right side";
	}
	else{
		char Direction = *(argv[1]);
		char wall = *(argv[2]);
		//Preprocessing: mesh and group creation
		int spaceDim = 2;
		
		// Prepare for the mesh
		cout << "Building mesh" << endl;
		cout << "Construction of a cartesian mesh" << endl;
		double inf = 0.0;
		double sup = 1.0;
		double discontinuity;
		int nx, ny, ncells;
		if (Direction == 'x'){
			nx=50	;
			ny=2;
			discontinuity = (inf + sup)/2.0 +  0.75/nx;
			ncells = nx;
			
		}
		else if (Direction == 'y'){
			nx=2;
			ny=50;
			discontinuity = (inf + sup)/2.0 +  0.75/ny;
			ncells = ny;
		}

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
		myProblem.setPeriodicFaces(M, Direction, ncells );
		
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
			
			double coordLeft, coordRight, coordFace; 
			if(Fj.getNumberOfCells()==2 ){ 
				myProblem.setOrientation(j,vec_normal_sigma);
				myProblem.setInteriorIndex(j);
				Cell Ctemp2 = M.getCell(idCells[1]);
				if (Direction == 'x'){
					coordLeft = Ctemp1.x();
					coordRight = Ctemp2.x();
					coordFace = Fj.x();
					
				}
				else if (Direction == 'y'){
					coordLeft = Ctemp1.y();
					coordRight = Ctemp2.y();
					coordFace = Fj.y() ;
				}
				Pressure0[idCells[0]] = initialPressure(coordLeft,discontinuity);
				Pressure0[idCells[1]] = initialPressure(coordRight,discontinuity);
				Velocity0[j] = dotprod(initialVelocity(coordFace, discontinuity, Direction),vec_normal_sigma );
			}
			else if (Fj.getNumberOfCells()==1  ){ // If boundary face and if periodic check that the boundary face is the computed (avoid passing twice ) 
				for (int idim = 0; idim <spaceDim; idim ++){
						if (vec_normal_sigma[idim] < 0)
							vec_normal_sigma[idim] = -vec_normal_sigma[idim];
				}
				myProblem.setOrientation(j,vec_normal_sigma);
				if (Direction == 'x')  coordFace = Fj.x();
				else if (Direction == 'y') coordFace = Fj.y() ;

				if  (myProblem.IsFaceBoundaryNotComputedInPeriodic(j) == false && myProblem.IsFaceBoundaryComputedInPeriodic(j) == false){
					if (wall =='l' && coordFace <1.0/(4*ncells) ){
						myProblem.setWallBoundIndex(j);
						wallVelocityMap[j] = 0;
					}
					else if (wall =='r' && abs(coordFace-1) <1.0/(4*ncells) ){
						myProblem.setWallBoundIndex(j);
						wallVelocityMap[j] = 0;
					}
					else {
						myProblem.setSteggerBoundIndex(j);	
						wallVelocityMap[j] = dotprod(initialVelocity(coordFace, discontinuity, Direction),vec_normal_sigma ) ;
						wallPressureMap[j] = initialPressure(coordFace,discontinuity);
					}
				}
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
		unsigned MaxNbOfTimeStep = 100000;
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

		// Should check if tmax, ncells, cfl and pl, pr, ul, ur are the same 

		/* cout << "Python script for exact solution" << endl;
		int result = system("python3 EulerBarotropicStaggered_1DRiemannProblem.py");  
		if (result == 0) {
			cout << "Script executed" << endl;
		} else {
			cerr << "ERROR in execution python script" << endl;
		} */
	}
		
	return EXIT_SUCCESS;
}
