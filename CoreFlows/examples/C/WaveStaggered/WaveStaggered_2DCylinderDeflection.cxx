#include "WaveStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

std::vector<double> ExactVelocity(double r, double theta, double r1, double r0){
	std::vector<double> vec(2);
	vec[0] = r1*r1/(r1*r1 -r0*r0)*(1 - r0*r0/(r*r) * cos(2*theta));
	vec[1] =  r1*r1/(r1*r1 -r0*r0)*(- r0*r0/(r*r) * sin(2*theta));
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
	vec[0] = 1;
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
	if(argc<2)
	{
		    cout << "- DOMAIN : SQUARE" << endl;
		    cout << "- MESH : CARTESIAN, GENERATED INTERNALLY WITH CDMATH" << endl<< endl;
		    cout << "Construction of a cartesian mesh" << endl;
	    double xinf=-0.5;
	    double xsup= 0.5;
	    double yinf=-0.5;
	    double ysup= 0.5;
	    int nx=50;
	    int ny=50;
	    M=Mesh(xinf,xsup,nx,yinf,ysup,ny);
	    double eps=1e-6;
	    M.setGroupAtPlan(xsup,0,eps,"RightEdge");
	    M.setGroupAtPlan(xinf,0,eps,"LeftEdge");
	    M.setGroupAtPlan(yinf,1,eps,"BottomEdge");
	    M.setGroupAtPlan(ysup,1,eps,"TopEdge");
	}
	else
	{
		// ./resources/AnnulusSpiderWeb5x16.med
	    cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
	    cout << "Loading of a mesh named "<<argv[1] << endl;
	    string filename = argv[1];
	    M=Mesh(filename);
	}

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

	Field ExactVelocityAtCells("ExactVelocityAtCells", CELLS, M, 3); //TODO not used ?
	Field ExactVelocityAtFaces("ExactVelocityAtFaces", FACES, M, 1);
	Field ExactVelocityInterpolate("ExactVelocityInterpolate", CELLS, M, 3);
	
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
		double r =  sqrt(Fj.x()*Fj.x() + Fj.y()*Fj.y());
		double theta = atan(Fj.y()/Fj.x());
		if (fabs(theta)<1e-10 && Fj.x() > 1e-10){ //TODO suppression ?
			for (int idim = 0; idim < spaceDim; ++idim)
				vec_normal_sigma[idim] *=-1;
		}
		myProblem.setOrientation(j,vec_normal_sigma);
		
		ExactVelocityAtFaces[j] = dotprod( ExactVelocity(r, theta, r1, r0), vec_normal_sigma); 

		if(Fj.getNumberOfCells()==2){
			Cell Ctemp2 = M.getCell(idCells[1]);
			myProblem.setInteriorIndex(j);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y());
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp2.y());
			Velocity0[j] = dotprod( initialVelocity(Fj.x(),Fj.y()) ,vec_normal_sigma);
		}
		else if (Fj.getNumberOfCells()==1){
			if (( sqrt( Fj.x()*Fj.x()+ Fj.y()*Fj.y() )  ) <= (r0 +r1)/2.0 ){// if face is on interior (wallbound condition) r_int = 1.2 ou 0.8 selon le maillage
				myProblem.setWallBoundIndex(j);
				wallVelocityMap[j] =  0;
			}
			else {// if face is on exterior (stegger condition) 			
				myProblem.setSteggerBoundIndex(j);								
				wallVelocityMap[j] = dotprod( initialBoundVelocity(Fj.x(),Fj.y()), vec_normal_sigma);
				wallPressureMap[j] = initialBoundPressure(Fj.x(),Fj.y());
			} // building exact solution at faces and its interpolation at cell	
			ExactVelocityAtFaces[j] = wallVelocityMap[j];
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
	unsigned MaxNbOfTimeStep = 1	;
	int freqSave = 300;
	double cfl = 0.5;
	double maxTime = 100000;
	double precision = 1e-12;

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
	myProblem.InterpolateFromFacesToCells(ExactVelocityAtFaces, ExactVelocityInterpolate);
	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();
	

	return EXIT_SUCCESS;
}
