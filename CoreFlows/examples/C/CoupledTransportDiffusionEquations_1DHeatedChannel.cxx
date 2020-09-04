#include "TransportEquation.hxx"
#include "DiffusionEquation.hxx"

using namespace std;

#define PI 3.14159265

void power_field_CoupledTransportDiffusionTest(Field & Phi){
	double L=4.2;
	double lambda=0.2;
	double phi=1e5;
	double x;
	Mesh M = Phi.getMesh();
	int nbCells = M.getNumberOfCells();
	for (int j = 0; j < nbCells; j++) {
		x=M.getCell(j).x();
		Phi(j) = phi*cos(PI*(x-L/2)/(L+lambda));
	}
}

int main(int argc, char** argv)
{
	//Preprocessing: mesh and group creation
	double xinf=0.0;
	double xsup=4.2;
	int nx=100;
	double eps=1.E-6;
	cout << "Building of the diffusion mesh with "<<nx<<" cells" << endl;
	Mesh diffusionMesh(xinf,xsup,nx);
	diffusionMesh.setGroupAtPlan(xsup,0,eps,"Neumann");
	diffusionMesh.setGroupAtPlan(xinf,0,eps,"Neumann");

	cout << "Building of the transport mesh with "<<nx<<" cells" << endl;
	Mesh transportMesh(xinf,xsup,nx);
	transportMesh.setGroupAtPlan(xsup,0,eps,"Neumann");
	transportMesh.setGroupAtPlan(xinf,0,eps,"Inlet");
	int spaceDim = 1;

	// Boundary conditions 
	map<string, LimitFieldDiffusion> boundaryFieldsDiffusion;
	map<string, LimitFieldTransport> boundaryFieldsTransport;

	// Boundary conditions for the solid
	LimitFieldDiffusion limitNeumann;
	limitNeumann.bcType=NeumannDiffusion;
	boundaryFieldsDiffusion["Neumann"] = limitNeumann;

	// Boundary conditions for the fluid
	LimitFieldTransport limitInlet;
	limitInlet.bcType=InletTransport;
	limitInlet.h =1.3e6;//Inlet water enthalpy
	boundaryFieldsTransport["Inlet"] = limitInlet;

	//Set the fluid transport velocity
	vector<double> transportVelocity(1,5);//Vitesse du fluide

	//Solid parameters
	double cp_ur=300;//Uranium specific heat
	double rho_ur=10000;//Uranium density
	double lambda_ur=5;

	TransportEquation  myTransportEquation(LiquidPhase, around155bars600KTransport,transportVelocity);
	Field fluidEnthalpy("Enthalpie", CELLS, transportMesh, 1);
	bool FECalculation=false;
    DiffusionEquation  myDiffusionEquation(spaceDim,FECalculation,rho_ur, cp_ur, lambda_ur);

	Field solidTemp("Solid temperature", CELLS, diffusionMesh, 1);
	Field fluidTemp("Fluid temperature", CELLS, transportMesh, 1);

	double heatTransfertCoeff=10000;//fluid/solid heat exchange coefficient
	myTransportEquation.setHeatTransfertCoeff(heatTransfertCoeff);
	myDiffusionEquation.setHeatTransfertCoeff(heatTransfertCoeff);

	//Set heat source in the solid
	Field Phi("Heat power field", CELLS, diffusionMesh, 1);
	power_field_CoupledTransportDiffusionTest(Phi);
	myDiffusionEquation.setHeatPowerField(Phi);
	Phi.writeVTK("1DheatPowerField");

	//Initial field creation
	Vector VV_Constant(1);
	VV_Constant(0) = 623;//Rod clad temperature nucleaire

	cout << "Construction de la condition initiale " << endl;
	// generate initial condition
	myDiffusionEquation.setInitialFieldConstant(diffusionMesh,VV_Constant);


	VV_Constant(0) = 1.3e6;
	myTransportEquation.setInitialFieldConstant(transportMesh,VV_Constant);

	//set the boundary conditions
	myTransportEquation.setBoundaryFields(boundaryFieldsTransport);//Neumann and Inlet BC will be used
	myDiffusionEquation.setBoundaryFields(boundaryFieldsDiffusion);//Only Neumann BC will be used

	// set the numerical method
	myDiffusionEquation.setTimeScheme( Explicit);
	myTransportEquation.setTimeScheme( Explicit);

	// name result file
	string fluidFileName = "1DFluidEnthalpy";
	string solidFileName = "1DSolidTemperature";

	// parameters calculation
	unsigned MaxNbOfTimeStep =3;
	int freqSave = 10;
	double cfl = 0.5;
	double maxTime = 1000000;
	double precision = 1e-6;

	myDiffusionEquation.setCFL(cfl);
	myDiffusionEquation.setPrecision(precision);
	myDiffusionEquation.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myDiffusionEquation.setTimeMax(maxTime);
	myDiffusionEquation.setFreqSave(freqSave);
	myDiffusionEquation.setFileName(solidFileName);

	myTransportEquation.setCFL(cfl);
	myTransportEquation.setPrecision(precision);
	myTransportEquation.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myTransportEquation.setTimeMax(maxTime);
	myTransportEquation.setFreqSave(freqSave);
	myTransportEquation.setFileName(fluidFileName);

	// loop on time steps
	myDiffusionEquation.initialize();
	myTransportEquation.initialize();

	double time=0,dt=0;
	int nbTimeStep=0;
	bool stop=false, stop_transport=false, stop_diffusion=false; // Does the Problem want to stop (error) ?
	bool ok; // Is the time interval successfully solved ?

	// Time step loop
	while(!stop && !(myDiffusionEquation.isStationary() && myTransportEquation.isStationary()) &&time<maxTime && nbTimeStep<MaxNbOfTimeStep)
	{
		ok=false; // Is the time interval successfully solved ?
		fluidTemp=myTransportEquation.getFluidTemperatureField();
		solidTemp=myDiffusionEquation.getRodTemperatureField();
		myDiffusionEquation.setFluidTemperatureField(fluidTemp);
		myTransportEquation.setRodTemperatureField(solidTemp);
		// Guess the next time step length
		dt=min(myDiffusionEquation.computeTimeStep(stop),myTransportEquation.computeTimeStep(stop));
		if (stop){
			cout << "Failed computing time step "<<nbTimeStep<<", time = " << time <<", dt= "<<dt<<", stopping calculation"<< endl;
		break;
		}
		// Loop on the time interval tries
		while (!ok && !stop )
		{
			stop_transport=!myTransportEquation.initTimeStep(dt);
			stop_diffusion=!myDiffusionEquation.initTimeStep(dt);
			stop=stop_diffusion && stop_transport;

			// Prepare the next time step
			if (stop){
				cout << "Failed initializing time step "<<nbTimeStep<<", time = " << time <<", dt= "<<dt<<", stopping calculation"<< endl;
			break;
			}
			// Solve the next time step
			ok=myDiffusionEquation.solveTimeStep()&& myTransportEquation.solveTimeStep();

			if (!ok)   // The resolution failed, try with a new time interval.
			{
				myDiffusionEquation.abortTimeStep();
				myTransportEquation.abortTimeStep();
					cout << "Failed solving time step "<<nbTimeStep<<", time = " << time<<" dt= "<<dt<<", cfl= "<<cfl <<", stopping calculation"<< endl;
					stop=true; // Impossible to solve the next time step, the Problem has given up
					break;
			}
			else // The resolution was successful, validate and go to the next time step.
			{
				cout << "Time step = "<< nbTimeStep << ", dt = "<< dt <<", time = "<<time << endl;
				myDiffusionEquation.validateTimeStep();
				myTransportEquation.validateTimeStep();
				time=myDiffusionEquation.presentTime();
				nbTimeStep++;
			}
		}
	}
	if(myDiffusionEquation.isStationary() && myTransportEquation.isStationary())
		cout << "Stationary state reached" <<endl;
	else if(time>=maxTime)
		cout<<"Maximum time "<<maxTime<<" reached"<<endl;
	else if(nbTimeStep>=MaxNbOfTimeStep)
		cout<<"Maximum number of time steps "<<MaxNbOfTimeStep<<" reached"<<endl;
	else
		cout<<"Error problem wants to stop!"<<endl;

	cout << "End of calculation time t= " << time << " at time step number "<< nbTimeStep << endl;
	if (ok)
		cout << "Coupled simulation "<<fluidFileName<<" and "<<solidFileName<<" was successful !" << endl;
	else
		cout << "Coupled simulation "<<fluidFileName<<" and "<<solidFileName<<"  failed ! " << endl;

	cout << "------------ End of calculation -----------" << endl;
	myDiffusionEquation.terminate();
	myTransportEquation.terminate();

	return EXIT_SUCCESS;
}
