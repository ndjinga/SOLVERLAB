#include "DriftModel.hxx"

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 3;
	// Prepare for the mesh
	cout << "Building cartesian mesh" << endl;
	double xinf = 0 ;
	double xsup=2.0;
	double yinf=0.0;
	double ysup=2.0;
	double zinf=0.0;
	double zsup=4.0;
	int nx=10;
	int ny=nx;
	int nz=20;
	double xcloison=(xinf+xsup)/2;
	double ycloison=(yinf+ysup)/2;
	double zcloisonmin=1;
	double zcloisonmax=3;
	Mesh M(xinf,xsup,nx,yinf,ysup,ny,zinf,zsup,nz);
	// set the limit field for each boundary
	double eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"wall");
	M.setGroupAtPlan(xinf,0,eps,"wall");
	M.setGroupAtPlan(ysup,1,eps,"wall");
	M.setGroupAtPlan(yinf,1,eps,"wall");
	M.setGroupAtPlan(zsup,2,eps,"outlet");
	M.setGroupAtPlan(zinf,2,eps,"inlet");
	double dx=(xsup-xinf)/nx;
	double dy=(ysup-yinf)/ny;
	double dz=(zsup-zinf)/nz;
	int ncloison=nz*(zcloisonmax-zcloisonmin)/(zsup-zinf);
	int i=0	;
	int j=0;
	while( i< ncloison){
		while(j< ny){
			M.setGroupAtFaceByCoords(xcloison,(j+0.5)*dy,zcloisonmin+(i+0.5)*dz,eps,"wall");
			M.setGroupAtFaceByCoords((j+0.5)*dx,ycloison,zcloisonmin+(i+0.5)*dz,eps,"wall");
			j=j+1;
		}
		i=i+1;
	}

    // set the limit field for each boundary
	double wallVelocityX=0;
	double wallVelocityY=0;
	double wallVelocityZ=0;
	double wallTemperature=573;
	double inletConc=0;
	double inletVelocityX=0;
	double inletVelocityY=0;
	double inletVelocityZ=1;
	double inletTemperature=563;
	double outletPressure=155e5;

    // physical constants
	vector<double> gravite = vector<double>(spaceDim,0);

	gravite[0]=0;
	gravite[1]=0;
	gravite[2]=-10;

	double heatPower1=0;
	double heatPower2=0.25e8;
	double heatPower3=0.5e8;
	double heatPower4=1e8;

	DriftModel myProblem = DriftModel(around155bars600K,spaceDim);
	int nVar =myProblem.getNumberOfVariables();
	Field heatPowerField=Field("heatPowerField", CELLS, M, 1);

	int nbCells=M.getNumberOfCells();

	for (int i=0;i<nbCells;i++){
		double x=M.getCell(i).x();
		double y=M.getCell(i).y();
		double z=M.getCell(i).z();
		if (z> zcloisonmin && z< zcloisonmax)
			if (y<ycloison && x<xcloison)
				heatPowerField[i]=heatPower1;
			if (y<ycloison && x>xcloison)
				heatPowerField[i]=heatPower2;
			if (y>ycloison && x<xcloison)
				heatPowerField[i]=heatPower3;
			if (y>ycloison && x>xcloison)
				heatPowerField[i]=heatPower4;
		else
			heatPowerField[i]=0;
	}
	heatPowerField.writeVTK("heatPowerField",true);

    //Prepare for the initial condition
	Vector VV_Constant =Vector(nVar);

	// constant vector
	VV_Constant[0] = inletConc ;
	VV_Constant[1] = outletPressure ;
	VV_Constant[2] = inletVelocityX;
	VV_Constant[3] = inletVelocityY;
	VV_Constant[4] = inletVelocityZ;
	VV_Constant[5] = inletTemperature ;

    //Initial field creation
	cout<<"Building initial data " <<endl;
	myProblem.setInitialFieldConstant(M,VV_Constant);

    // the boundary conditions
	myProblem.setOutletBoundaryCondition("outlet", outletPressure);
	myProblem.setInletBoundaryCondition("inlet", inletTemperature, inletConc, inletVelocityX, inletVelocityY, inletVelocityZ);
	myProblem.setWallBoundaryCondition("wall", wallTemperature, wallVelocityX, wallVelocityY,wallVelocityZ);

	// set physical parameters
	myProblem.setHeatPowerField(heatPowerField);
	myProblem.setGravity(gravite);

	// set the numerical method
	myProblem.setNumericalScheme(upwind, Explicit);
	myProblem.setWellBalancedCorrection(false);

	// name file save
	string fileName = "3DCanalCloison";

	// parameters calculation
	unsigned MaxNbOfTimeStep = 3;
	int freqSave = 1;
	double cfl = 0.3;
	double maxTime = 5;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
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
