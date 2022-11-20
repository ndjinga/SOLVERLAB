#include "DiffusionEquation.hxx"

using namespace std;

#define PI 3.14159265


int main(int argc, char** argv)
{
    bool FECalculation, NeumannBC;
    string SCHEME, BC;
    
	if(argc<3)
		throw CdmathException("Provide two boolean arguments : isFECalcuation and isNeumanBC");
	else
	{
	    SCHEME = argv[1];
	    BC = argv[2];
	    
	    FECalculation = (SCHEME=="FE" or SCHEME=="EF");
	    NeumannBC = (BC=="Neumann" or BC=="NEUMANN");
	}

	cout << "-- RESOLUTION OF THE 2D DIFFUSION SYSTEM ON A SQUARE"<<endl;
	if(FECalculation )
		cout << "- Numerical scheme : implicit finite element scheme" << endl;
	else
		cout << "- Numerical scheme : implicit finite volume scheme" << endl;
	
	if(NeumannBC)
		cout << "- Boundary condition : Neumann" << endl;
	else
		cout << "- Boundary condition : Dirichlet" << endl;
	cout<<endl;
	
	//Preprocessing: mesh creation
	double xinf=0.0;
	double xsup=1.;
	int nx=2;
	double yinf=0.0;
	double ysup=1.;
	int ny=2;
	cout << "Building of a 2D mesh with "<<nx<<"x"<<ny<<" cells" << endl;
	Mesh square_mesh;
	if( FECalculation )
		square_mesh=Mesh(xinf,xsup,nx,yinf,ysup,ny,0);//right triangle mesh build from spliting the cells of a nx*ny cartesian mesh
	else
		square_mesh=Mesh(xinf,xsup,nx,yinf,ysup,ny);//nx*ny cartesian finite volume mesh
		
	if( square_mesh.getSpaceDimension()!=2 or square_mesh.getMeshDimension()!=2) 
		throw CdmathException("Wrong space or mesh dimension : space and mesh dimensions should be 2");
	if(FECalculation and not square_mesh.isTriangular()) 
		throw CdmathException("Wrong cell types : mesh is not made of triangles");
	
	int nbNodes = square_mesh.getNumberOfNodes();
	int nbCells = square_mesh.getNumberOfCells();
	int spaceDim = square_mesh.getSpaceDimension();

	//Solid parameters
	double cp_ur=300;//Uranium specific heat
	double rho_ur=10000;//Uranium density
	double lambda_ur=5;
 
	DiffusionEquation  myProblem(spaceDim,FECalculation,rho_ur,cp_ur,lambda_ur);

	//Setting initial field
	Vector VV_Constant(1);
	VV_Constant(0) = 20;//Solid temperature

	cout << "Building initial data" << endl;
	if(FECalculation )
		myProblem.setInitialFieldConstant(square_mesh,VV_Constant, NODES);
	else
		myProblem.setInitialFieldConstant(square_mesh,VV_Constant, CELLS);

	//Set heat source
	Field Phi;
	double x,y;
	if(FECalculation)
	{
		Phi=Field("Heat power field", NODES, square_mesh, 1);
		for (int j = 0; j < nbNodes; j++) {
			x=square_mesh.getNode(j).x();
			y=square_mesh.getNode(j).y();
			Phi(j) = sin(PI*x)*sin(PI*y);
		}
	}
	else
	{
		Phi=Field("Heat power field", CELLS, square_mesh, 1);
		for (int j = 0; j < nbCells; j++) {
			x=square_mesh.getCell(j).x();
			y=square_mesh.getCell(j).y();
			Phi(j) = sin(PI*x)*sin(PI*y);
		}
	}
	myProblem.setHeatPowerField(Phi);
	Phi.writeVTK("2DheatPowerField");

	//set the boundary conditions
	if( NeumannBC )
		myProblem.setNeumannBoundaryCondition("Boundary");
	else
		myProblem.setDirichletBoundaryCondition("Boundary");

	// set the numerical method
	myProblem.setTimeScheme( Implicit);

	// name result file
	string fileName;
	if( NeumannBC )
		if(FECalculation)
			fileName = "2D_SQUARE_FE_Implicit_Neumann";
		else
			fileName = "2D_SQUARE_FV_Implicit_Neumann";
	else
		if(FECalculation)
			fileName = "2D_SQUARE_FE_Implicit_Dirichlet";
		else
			fileName = "2D_SQUARE_FV_Implicit_Dirichlet";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 1;
	double cfl = 100;
	double maxTime = 1000000;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);

	//set file saving format
	myProblem.setSaveFileFormat(CSV);

	// evolution
	myProblem.initialize();
	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of simulation -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
