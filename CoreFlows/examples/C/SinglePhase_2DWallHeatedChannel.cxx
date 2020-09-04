#include "SinglePhase.hxx"

using namespace std;

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	double xinf=-0.005;
	double xsup= 0.005;
	double yinf=0.0;
	double ysup=2.0;
	int nx=50;
	int ny=50;
	cout << "Building a regular mesh with "<<nx<<" times "<< ny<< " cells " << endl;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny);
	double eps=1.E-6;
	M.setGroupAtPlan(xsup,0,eps,"wall");
	M.setGroupAtPlan(xinf,0,eps,"wall");
	M.setGroupAtPlan(yinf,1,eps,"Neumann");//Inlet
	M.setGroupAtPlan(ysup,1,eps,"Neumann");//Outlet
	int spaceDim = M.getSpaceDimension();

	// physical constants
	vector<double> viscosite(1), conductivite(1);
	viscosite[0]= 8.85e-5;
	conductivite[0]=1000;//transfert de chaleur du à l'ébullition en paroi.

	// set the limit field for each boundary
	LimitField limitWall;
	map<string, LimitField> boundaryFields;
	limitWall.bcType=Wall;
	limitWall.T = 623;//Temperature des parois chauffantes
	limitWall.p = 155e5;
	limitWall.v_x = vector<double>(1,0);
	limitWall.v_y = vector<double>(1,0);
	boundaryFields["wall"]= limitWall;

	LimitField limitInlet;
	limitInlet.bcType=Inlet;
	limitInlet.T = 573;//Temperature d'entree du fluide
	limitInlet.v_x = vector<double>(1,0);
	limitInlet.v_y = vector<double>(1,5);//Vitesse d'entree du fluide
	boundaryFields["Inlet"]= limitInlet;

	LimitField limitOutlet;
	limitOutlet.bcType=Outlet;
	limitOutlet.p = 155e5;
	boundaryFields["Outlet"]= limitOutlet;

	LimitField limitNeumann;
	limitNeumann.bcType=Neumann;
	boundaryFields["Neumann"] = limitNeumann;

	SinglePhase  myProblem(Liquid,around155bars600K,spaceDim);
	int nVar = myProblem.getNumberOfVariables();

	//Initial field creation
	cout << "Construction de la condition initiale" << endl;
	/* First case constant initial data */
	Vector VV_Constant(nVar);
	// constant vector
	VV_Constant(0) = 155e5;
	VV_Constant(1) = 0;
	VV_Constant(2) = 5;
	VV_Constant(3) = 573;
	/* Second case restart from a previous calculation */
/*
	string inputfile = "SinglePhasePrim_2DChannelWithViscosityWithConduction";//nom du fichier (sans le .med)
	int iter=50400;//numero du pas de temps a recuperer dans le fichier
	Field VV1(inputfile,CELLS,"P,vx,vy,T",iter,0);//Chargement des valeurs du champ
	for(int i=0;i<M.getNumberOfCells();i++)
		for(int j=0;j<nVar;j++)
			VV(i,j)=VV1(i,j);//recuperation des valeurs du champ
*/
	myProblem.setInitialFieldConstant(M,VV_Constant);

	//set the boundary conditions
	myProblem.setBoundaryFields(boundaryFields);

	// physical parameters
	myProblem.setViscosity(viscosite);
	myProblem.setConductivity(conductivite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Implicit);

	// set the Petsc resolution
	myProblem.setLinearSolver(GMRES,ILU,true);

	// name result file
	string fileName = "2DHeatedWallChannel";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 10;
	double maxTime = 50;
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
