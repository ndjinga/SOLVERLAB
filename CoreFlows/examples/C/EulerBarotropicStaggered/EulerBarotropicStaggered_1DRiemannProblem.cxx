#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

double initialPressure(double x){
	if (x < 0)
		return 1;//155e5;
	else
		return 2;
}

double initialVelocity(double x){
	if (x < 0)
		return 2;
	else
		return 2;
}
int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	cout << "Building Cartesian mesh" << endl;
	double xinf=-1.0;
	double xsup=1.0;
	double discontinuity = (xinf+xsup)/2.;
	int nx=100;
	Mesh M(xinf,xsup,nx);
	int spaceDim = M.getSpaceDimension();

	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	Field Pressure0("pressure", CELLS, M, 1);
	Field Velocity0("velocity", FACES, M, 1);

    //Initial field creation
	cout << "Building initial data " <<endl; 
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
		myProblem.setOrientation(j,vec_normal_sigma);

		if(Fj.getNumberOfCells()==2 ){ 
			myProblem.setInteriorIndex(j);
			Cell Ctemp2 = M.getCell(idCells[1]);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x());
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x());
			Velocity0[j] = initialVelocity(Fj.x());
		}
		else if (Fj.getNumberOfCells()==1  ){ 
			myProblem.setSteggerBoundIndex(j);								
			wallVelocityMap[j] =initialVelocity(Fj.x());
			wallPressureMap[j] = initialPressure(Fj.x());
		}
	}
	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Velocity0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Explicit);
    
    // name of result file
	string fileName = "EulerBarotropicStaggered_1DRiemannProblem";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 10000;
	int freqSave = 1;
	double cfl = 0.5;
	double maxTime = 0.03;
	double precision = 1e-13;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setSaveFileFormat(CSV);
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
