#include "WaveStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

std::vector<double> ExactVelocity(double r, double theta, double r1, double r0){
	std::vector<double> vec(2);
	vec[0] = r1*r1/(r1*r1 -r0*r0)*(1 - r0*r0/(r*r) * math.cos(2*theta));
	vec[1] =  r1*r1/(r1*r1 -r0*r0)*(- r0*r0/(r*r) * math.sin(2*theta));
	return vec;
}
double initialPressure( double x, double y){
	return 0;
}
double initialBoundPressure( double x, double y){
	return 0;
}

std::vector<double> initialVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = 0;
	vec[1]  =  0;
	return vec;
}
std::vector<double> initialBoundVelocity(double x, double y){
	std::vector<double> vec(2);
	vec[0] = 1;
	vec[1] =  0;
	return vec;
}

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation

	int spaceDim = 2;
	// Prepare for the mesh
	cout << "Building mesh" << endl;
	std::string inputfile="/volatile/catB/esteban/Solverlab/SOLVERLAB_SRC/CoreFlows/examples/resources/AnnulusSpiderWeb20x64.med"
	double r0 = 0.8;
	double r1 = 6;

	Mesh M=svl.Mesh(inputfile);
	double kappa = 1;
	double rho = 1;
	double c = math.sqrt(kappa/rho)
	myProblem = WaveStaggered(spaceDim,rho, kappa);

	// Prepare for the initial condition
	// set the boundary conditions
	
	
	//Initial field creation
	print("Building initial data " ); 
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	Field Pressure0("pressure", svl.CELLS, M, 1);
	Field Velocity0("velocity", svl.FACES, M, 1);
	Field ExactVelocityInftyAtCells("ExactVelocityInftyAtCells", svl.CELLS, M, 3); 
	Field ExactVelocityInftyAtFaces("ExactVelocityInftyAtFaces", svl.FACES, M, 1)
	Field ExactVelocityInftyInterpolate ("ExactVelocityInftyAtInterpolate", svl.CELLS, M, 3);
	for (int l=0 ; l< M.getNumberOfCells(); l++){
		for (int k =0; k< spaceDim ; k++)
			ExactVelocityInftyInterpolate[l, k] =0;
	}


    //Initial field creation
	cout << "Building initial data " <<endl; 
	
	for (int j=0; j< M.getNumberOfFaces(); j++ ){
		std::vector<int> idCells = Fj.getCellsId();
		std::vector<double> vec_normal_sigma(2) ; //TODO = 0!!
		Cell Ctemp1 = M.getCell(idCells[0]);
		for(int l=0; l<Ctemp1.getNumberOfFaces(); l++){//we look for l the index of the face Fj for the cell Ctemp1
			if (i == Ctemp1.getFacesId()[l]){
				for (int idim = 0; idim < _Ndim; ++idim)
					vec_normal_sigma[idim] = Ctemp1.getNormalVector(l,idim);
			}
		}

		for (int v= 0; v < idCells.size(); v++){
			Ctemp = _mesh.getCell(idCells[v]); 
			double orien = getOrientation(i,Ctemp);
			double r =  sqrt(Fj.x()*Fj.x() + Fj.y()*Fj.y());
			double theta = math.arctan(Fj.y()/Fj.x());
			ExactVelocityInftyAtFaces[j] = ExactVelocity(r, theta, r1, r0) vec_normal_sigma //TODO dot product
			if ( Fj.getNumberOfCells() ==1){
				if ( r<= (r0 + r1)/2.0){ 
					myProblem.setWallBoundIndex(j);
					wallVelocityMap[j] = 0;
				}
				else {
					wallPressureMap[idCells[v]] = initialBOundPressure(Ctemp.x(),Ctemp1.y());
					wallVelocityMap[j] = initialBoundVelocity(Fj.x(),Fj.y()) vec_normal_sigma; //TODO dot product 
				}
			}
			else{
				Pressure0[idCells[v]] = initialPressure(Ctemp.x(),Ctemp1.y());
				Velocity0[j] = initialVelocity(Fj.x(),Fj.y()) vec_normal_sigma; //TODO dot product 
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
	string fileName = "WaveStaggered_2DCylinderDeflection";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 1;
	int freqSave = 1;
	double cfl = 0.6;
	double maxTime = 500;
	double precision = 1e-8;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(svl.VTK)
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
