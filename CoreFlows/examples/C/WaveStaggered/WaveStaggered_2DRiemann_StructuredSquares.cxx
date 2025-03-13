#include "WaveStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

double initialPressure( double z, double discontinuity){
	if (z < discontinuity)
		return 2;
	else
		return 1;
}

std::vector<double> initialVelocity(double z, double discontinuity, char Direction){
	std::vector<double> vec(2);
	if (z < discontinuity){
		if (Direction == 'x'){
			vec[0] = 1;
			vec[1] = 0;
		}
		if (Direction == 'y'){
			vec[0] = 0;
			vec[1] = 0;
		}
	}
	else{
		if (Direction == 'x'){
			vec[0] = 1;
			vec[1] = 0;
		}
		if (Direction == 'y'){
			vec[0] = 0;
			vec[1] = 1;
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


int main(int argc, char** argv){
	if (argc<2 || (*(argv[1]) != 'x' && *(argv[1]) != 'y') ){
		cout << "ERROR : you have to give a direction for the pseudo 1d Riemann problem, either 'x' or 'y' ";
	}
	else{
		char Direction = *(argv[1]);
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
			nx=200;
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
		
		double kappa = 1;
		double rho = 1;
		double c = sqrt(kappa/rho);
		WaveStaggered myProblem(spaceDim,rho, kappa);

		// Prepare for the initial condition
		// set the boundary conditions
		
		//Initial field creation
		cout << "Building initial data" << endl;
		std::map<int ,double> wallPressureMap;
		std::map<int ,double> wallVelocityMap ;
		Field Pressure0("pressure", CELLS, M, 1);
		Field Velocity0("velocity", FACES, M, 1);
		
		myProblem.setPeriodicFaces(M, Direction, ncells);
		std::map<int,int> FacePeriodicMap = myProblem.getFacePeriodicMap();
		
		for (int j=0; j< M.getNumberOfFaces(); j++ ){
			Face Fj = M.getFace(j);
			std::vector<int> idCells = Fj.getCellsId();
			std::vector<double> vec_normal_sigma(2) ; //TODO = 0!!
			Cell Ctemp1 = M.getCell(idCells[0]);
			for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
				if (j == Ctemp1.getFacesId()[l]){
					for (int idim = 0; idim < spaceDim; ++idim)
						vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
				}
			}
			
			double coordLeft, coordRight, coordFace; 
			if(Fj.getNumberOfCells()==2 ){ // myProblem.IsFaceBoundaryComputedInPeriodic(j)
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
				std::vector<double > InitialVel = initialVelocity(coordFace, discontinuity, Direction);
				double dotprod = 0;
				for (int k = 0 ; k <InitialVel.size() ; k++)
						dotprod += InitialVel[k] * vec_normal_sigma[k];
				Velocity0[j] = dotprod;
			}
			else if (Fj.getNumberOfCells()==1  ){ // If boundary face and if periodic check that the boundary face is the computed (avoid passing twice ) 
				for (int idim = 0; idim <spaceDim; idim ++){
						if (vec_normal_sigma[idim] < 0)
							vec_normal_sigma[idim] = -vec_normal_sigma[idim];
				}
				myProblem.setOrientation(j,vec_normal_sigma);
				if  (myProblem.IsFaceBoundaryNotComputedInPeriodic(j) == false && myProblem.IsFaceBoundaryComputedInPeriodic(j) == false)
					myProblem.setSteggerBoundIndex(j);	
				if (Direction == 'x'){
					coordFace = Fj.x();
					
				}
				else if (Direction == 'y'){
					coordFace = Fj.y() ;
				}							
				std::vector<double > BoundaryVel = initialVelocity(coordFace,discontinuity, Direction);
				double dotprod = 0;
				for (int k = 0 ; k <BoundaryVel.size() ; k++)
					dotprod += BoundaryVel[k] * vec_normal_sigma[k];
				wallVelocityMap[j] = dotprod;
				wallPressureMap[j] = initialPressure(coordFace,discontinuity);
			}
		}
		
		myProblem.setInitialField(Pressure0);
		myProblem.setInitialField(Velocity0);
		myProblem.setboundaryPressure(wallPressureMap);
		myProblem.setboundaryVelocity(wallVelocityMap);

		// set the numerical method
		myProblem.setTimeScheme(Explicit);
		
		// name of result file
		string fileName = "WaveStaggered_2DRiemann_StructuredSquares";

		// parameters calculation
		unsigned MaxNbOfTimeStep = 1000000;
		int freqSave = 1;
		double cfl = 0.5;
		double maxTime = 0.07;
		double precision = 1e-6;

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
	}

	return EXIT_SUCCESS;
}
