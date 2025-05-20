#include "EulerBarotropicStaggered.hxx"
#include "math.h"
#include <cassert>

using namespace std;

double initialPressure( double x, double y){
	return 432;
}

std::vector<double> initialVelocity(double x,double y){
	std::vector<double> vec(2);
	vec[0] = -667.18; 
	vec[1] = 1888.0789;
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
	int spaceDim = 2;
	Mesh M;
	double a = 1.0;
	double gamma = 2.0;
	EulerBarotropicStaggered myProblem = EulerBarotropicStaggered(GasStaggered, around1bar300K, a, gamma, spaceDim );
	if (argc>1  ){
		// ./resources/AnnulusSpiderWeb5x16.med or ./resources/AnnulusTriangles60.med
		cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
		cout << "Loading of a mesh named "<<argv[1] << endl;
		string filename = argv[1];
		M=Mesh(filename);

	}
	else{
		PetscInitialize(&argc, &argv, 0,0);
		// Prepare for the mesh
		cout << "Building mesh" << endl;
		cout << "Construction of a cartesian mesh" << endl;
		double inf = 0.0;
		double sup = 1.0;
		int ncells = 10;
		M=Mesh(inf,sup,ncells,inf,sup,ncells);

		assert(fabs(inf)<1e-11);
		assert(fabs(sup - 1.0)<1e-11);
		//myProblem.setPeriodicFaces(M, 'x', ncells); //Only works on [0,1]² -> not useful to adapt //TODO
	}

	//Initial field creation
	cout << "Building initial data" << endl;
	std::map<int ,double> wallPressureMap;
	std::map<int ,double> wallVelocityMap ;
	std::vector<double> wallVelocityVector(spaceDim);
	Field Pressure0("pressure", CELLS, M, 1);
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
		myProblem.setOrientation(j,vec_normal_sigma);
		if(Fj.getNumberOfCells()==2){
			Cell Ctemp2 = M.getCell(idCells[1]);
			myProblem.setInteriorIndex(j);
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x(),Ctemp1.y());
			Pressure0[idCells[1]] = initialPressure(Ctemp2.x(),Ctemp2.y());
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * ( initialPressure(Ctemp1.x(),Ctemp1.y()) + initialPressure(Ctemp2.x(),Ctemp2.y())  )/2;
			}
		else if (Fj.getNumberOfCells()==1){			
			Pressure0[idCells[0]] = initialPressure(Ctemp1.x() , Ctemp1.y() );
			// Interpolation must be consistent with the one done in UpdateDualDensity
			Momentum0[j] = dotprod(initialVelocity(Fj.x(),Fj.y()),vec_normal_sigma ) * (initialPressure(Ctemp1.x() , Ctemp1.y() ) + initialPressure(Fj.x(),Fj.x()))/2.0;
			if (myProblem.IsFaceBoundaryNotComputedInPeriodic(j) == false && myProblem.IsFaceBoundaryComputedInPeriodic(j) == false)
				myProblem.setSteggerBoundIndex(j);	
			// Boundary normal velocity, pressure and full velocity vector
			wallVelocityMap[j] = dotprod(initialVelocity(Fj.x(), Fj.y()),vec_normal_sigma );
			wallPressureMap[j] = initialPressure(Fj.x(), Fj.y()) ;
			for (int idm = 0; idm <spaceDim; idm ++)
				wallVelocityVector[idm] = initialVelocity(Fj.x(), Fj.y())[idm];
			myProblem.setboundaryVelocityVector(j, wallVelocityVector);
		}
	}

	myProblem.setInitialField(Pressure0);
	myProblem.setInitialField(Momentum0);
	myProblem.setboundaryPressure(wallPressureMap);
	myProblem.setboundaryVelocity(wallVelocityMap);

    // set the numerical method
	myProblem.setTimeScheme(Implicit);
    
    // name of result file
	string fileName = "EulerBarotropicStaggered_2DStatio";

    // parameters calculation
	unsigned MaxNbOfTimeStep = 1	;
	int freqSave = 1		;
	double cfl = 0.99;
	double maxTime = 50;
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
	Field ExactVelocityInterpolate("ExactVelocityInterpolate", CELLS, M, 3);
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
