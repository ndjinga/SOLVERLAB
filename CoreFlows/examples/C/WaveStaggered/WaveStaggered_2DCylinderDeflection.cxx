#include "WaveStaggered.hxx"
#include "math.h"
#include <cassert>
#include <iomanip>	

using namespace std;

std::vector<double> ExactVelocity(double r, double theta, double r1, double r0){
	std::vector<double> vec(2);
	vec[0] = pow(r1,2)/(pow(r1,2) - pow(r0,2))*(1 - pow(r0,2)/pow(r,2) * cos(2*theta)); 
	vec[1] = pow(r1,2)/(pow(r1,2) - pow(r0,2))*(  - pow(r0,2)/pow(r,2) * sin(2*theta)); 
	return vec;
}

double initialPressure( double x, double y){ return 0; }

std::vector<double> initialVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = 0; 
	vec[1] = 0;
	return vec;
}
std::vector<double> initialBoundVelocity(double x, double y){
	std::vector<double> vec(2);
	vec[0] = 1;
	vec[1] = 0;
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
	
	// Prepare for the mesh
	cout << "Building mesh" << endl;
	double r0 = 0.8;
	double r1 = 6;

	Mesh M;
	// ./resources/AnnulusSpiderWeb5x16.med
	cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
	cout << "Loading of a mesh named "<<argv[1] << endl;
	string filename = argv[1];
	M=Mesh(filename);

	double kappa = 1;
	double rho = 1;
	WaveStaggered myProblem(spaceDim,rho, kappa);

	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	Field Pressure0("pressure", CELLS, M, 1);
	Field Velocity0("velocity", FACES, M, 1);
	std::vector<double> ExactVelocityAtFaces(M.getNumberOfFaces());
	
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		Face Fj = M.getFace(j);
		std::vector<int> idCells = Fj.getCellsId();
		std::vector<double> vec_normal_sigma(2, 0.0) ; 
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
		Velocity0[j] = dotprod( initialVelocity(Fj.x(),Fj.y()) ,vec_normal_sigma);
		Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y());
		
		if(Fj.getNumberOfCells()==2){
			myProblem.setInteriorIndex(j);
			Pressure0[idCells[1]] = initialPressure(M.getCell(idCells[1]).x(),M.getCell(idCells[1]).y());
		}
		else if (Fj.getNumberOfCells()==1){
			if (( sqrt( Fj.x()*Fj.x()+ Fj.y()*Fj.y() )  ) <= (r0 +r1)/2.0 ){// if face is on interior (wallbound condition) r_int = 1.2 ou 0.8 selon le maillage
				myProblem.setWallBoundIndex(j);
				wallVelocityMap[j] =  0;
			}
			else {// if face is on exterior (stegger condition) 			
				myProblem.setSteggerBoundIndex(j);								
				wallVelocityMap[j] = dotprod( initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma);
				wallPressureMap[j] = initialPressure(Fj.x(),Fj.y());
			}
		}
	}
	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit	);    
    // name of result file
	string fileName = "WaveStaggered_2DCylinderDeflection";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 100000000	;
	int freqSave = 4000;
	double cfl = 0.99;
	double maxTime = 200;
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

	myProblem.addRotRot(false);
	
	// evolution
	myProblem.initialize();
	myProblem.InterpolateFromFacesToCells(ExactVelocityAtFaces);
	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl; 

	cout << "------------ End of calculation !!! -----------" << endl;
	cout << "\nBoundary Velocity error = "<< std::setprecision(17) << std::fixed<< myProblem.ErrorInftyVelocityBoundary(wallVelocityMap)<<endl;
	cout << "Error L2 of velocity at the faces = "<< myProblem.ErrorVelocity(ExactVelocityAtFaces)[0] <<endl;
	//cout << "Error L2 of interpolated velocity at cells= "<< myProblem.ErrorL2VelocityInfty(ExactVelocityAtFaces, ExactVelocityInterpolate)[1] <<endl;
	myProblem.terminate();
	

	return EXIT_SUCCESS;
}
