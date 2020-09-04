#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group importation
	cout << "Reading a mesh with sudden cross-section change for test SinglePhase_2DWallHeatedChannel_ChangeSect()" << endl;
	Mesh M("resources/VaryingSectionDuct.med");

	// Conditions aux limites 
	//Bords externes
	double xinf=0.0;
	double xsup=0.01;
	double yinf=0.0;
	double ysup=0.01;
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"Wall");
	M.setGroupAtPlan(xinf,0,eps,"Wall");
	M.setGroupAtPlan(yinf,1,eps,"Inlet");//
	M.setGroupAtPlan(ysup,1,eps,"Outlet");//
	//Bords internes
	int nx=60, ny=60;//Nombre de cellules utilisees au depart dans Salome ou Alamos
	double dx = (xsup-xinf)/nx, dy = (ysup-yinf)/ny;//taille d'une cellule
	for(int i=0; i<ny/2;i++){
		 M.setGroupAtFaceByCoords((xsup-xinf)/4,(ysup-yinf)/4+(i+0.5)*dy,0,eps,"Wall");//Paroi verticale intérieure gauche
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4,(ysup-yinf)/4+(i+0.5)*dy,0,eps,"Wall");//Paroi verticale intérieure droitee
	}
	for(int i=0; i<nx/4;i++){
		 M.setGroupAtFaceByCoords((i+0.5)*dx,(ysup-yinf)/4,0,eps,"Wall");//paroi horizontale en bas à gauche
		 M.setGroupAtFaceByCoords((i+0.5)*dx,(ysup-yinf)*3/4,0,eps,"Wall");//paroi horizontale en haut à gauche
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4+(i+0.5)*dx,(ysup-yinf)/4,0,eps,"Wall");//paroi horizontale en bas à droite
		 M.setGroupAtFaceByCoords((xsup-xinf)*3/4+(i+0.5)*dx,(ysup-yinf)*3/4,0,eps,"Wall");//paroi horizontale en haut à droite
	}

	int spaceDim = M.getSpaceDimension();

	// set the limit field for each boundary
	LimitField limitWall;
	map<string, LimitField> boundaryFields;
	limitWall.bcType=Wall;
	limitWall.T = 623;//Temperature des parois chauffantes
	limitWall.p = 155e5;
	limitWall.v_x = vector<double>(1,0);
	limitWall.v_y = vector<double>(1,0);
	boundaryFields["Wall"]= limitWall;

	LimitField limitInlet;
	limitInlet.bcType=Inlet;
	limitInlet.T = 573;//Temperature d'entree du fluide
	limitInlet.v_x = vector<double>(1,0);
	limitInlet.v_y = vector<double>(1,2.5);//Vitesse d'entree du fluide
	boundaryFields["Inlet"]= limitInlet;

	LimitField limitOutlet;
	limitOutlet.bcType=Outlet;
	limitOutlet.p = 155e5;
	boundaryFields["Outlet"]= limitOutlet;

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	// Prepare for the initial condition
	int nVar = myProblem.getNumberOfVariables();
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 155e5;
	VV_Constant(1) = 0;
	VV_Constant(2) = 2.5;
	VV_Constant(3) = 573;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

	// physical constants
	vector<double> viscosite(1), conductivite(1);
	viscosite[0]= 8.85e-5;
	conductivite[0]=1000;//transfert de chaleur du à l'ébullition en paroi.

	//Set boundary values
	myProblem.setBoundaryFields(boundaryFields);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);

	// name result file
	string fileName = "2DWallHeatedChannel_ChangeSect";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl =.5;
	double maxTime = 500;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,50);
	myProblem.saveVelocity();

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
