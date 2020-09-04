#include "StationaryDiffusionEquation.hxx"
#include "Node.hxx"
#include "math.h"
#include <assert.h>

double pi = M_PI;
using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

	/* Mesh data */
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=20;
	int ny=20;

    /* Mesh construction */
	Mesh M(xinf,xsup,nx,yinf,ysup,ny,0); //Regular triangular mesh

	/* set the limit field for each boundary */
	double eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"Bord1");
	M.setGroupAtPlan(xinf,0,eps,"Bord2");
	M.setGroupAtPlan(ysup,1,eps,"Bord3");
	M.setGroupAtPlan(yinf,1,eps,"Bord4");

    /* set the boundary values for each boundary */
	double T1=0;
	double T2=0;
	double T3=0;
	double T4=0;

	cout<< "Building of a regular triangular 2D mesh from a square mesh with "<< nx<<"x" <<ny<< " cells"<<endl;

    /* Create the problem */
    bool FEComputation=false;
	StationaryDiffusionEquation myProblem(spaceDim,FEComputation);
	myProblem.setMesh(M);

    /* set the boundary conditions */
	myProblem.setNeumannBoundaryCondition("Bord1");
	myProblem.setNeumannBoundaryCondition("Bord2");
	myProblem.setNeumannBoundaryCondition("Bord3");
	myProblem.setNeumannBoundaryCondition("Bord4");

	/* Set the right hand side function*/
	Field my_RHSfield("RHS_field", CELLS, M, 1);
    Cell Ci; 
    double x, y;
	for(int i=0; i< M.getNumberOfCells(); i++)
    {
		Ci= M.getCell(i);
		x = Ci.x();
		y = Ci.y();

		my_RHSfield[i]=2*pi*pi*cos(pi*x)*cos(pi*y);//mettre la fonction definie au second membre de l'edp
	}
	myProblem.setHeatPowerField(my_RHSfield);
	myProblem.setLinearSolver(GMRES,ILU);

    /* name the result file */
	string fileName = "StationnaryDiffusion_2DFV_RegularTriangles_Neumann";
	myProblem.setFileName(fileName);

	/* Run the computation */
	myProblem.initialize();
	bool ok = myProblem.solveStationaryProblem();
	if (!ok)
		cout << "Simulation of "<<fileName<<" failed !" << endl;
	else
	{
		/********************** Postprocessing and measure od numerical error******************************/
		Field my_ResultField = myProblem.getOutputTemperatureField();
		/* The following formulas use the fact that the exact solution is equal the right hand side divided by 2*pi*pi */
		double max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(2*pi*pi);
		double max_sol_num=my_ResultField.max();
		double min_sol_num=my_ResultField.min();
		double erreur_abs=0;
		for(int i=0; i< M.getNumberOfCells() ; i++)
			if( erreur_abs < abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]) )
				erreur_abs = abs(my_RHSfield[i]/(2*pi*pi) - my_ResultField[i]);
		
		cout<<"Absolute error = max(| exact solution - numerical solution |) = "<<erreur_abs <<endl;
		cout<<"Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = "<<erreur_abs/max_abs_sol_exacte<<endl;
		cout<<"Maximum numerical solution = "<< max_sol_num<< " Minimum numerical solution = "<< min_sol_num << endl;
		
		assert( erreur_abs/max_abs_sol_exacte <1.);

        cout << "Simulation of "<<fileName<<" is successful ! " << endl;
    }

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
