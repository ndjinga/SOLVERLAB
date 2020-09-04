#include "StationaryDiffusionEquation.hxx"
#include "Node.hxx"
#include "math.h"
#include <assert.h>

double pi = M_PI;
using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 3;

	/* Mesh data */
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	double zinf=0.0;
	double zsup=1.0;
	int nx=2;
	int ny=2;
	int nz=2;

    /* Mesh construction : splitting polity to 0 yield all nodes considered boundary nodes */
	Mesh M(xinf,xsup,nx,yinf,ysup,ny,zinf,zsup,nz,0); //Regular tetrahadral mesh

	/* set the limit field for each boundary */
	double eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"Bord1");
	M.setGroupAtPlan(xinf,0,eps,"Bord2");
	M.setGroupAtPlan(ysup,1,eps,"Bord3");
	M.setGroupAtPlan(yinf,1,eps,"Bord4");
	M.setGroupAtPlan(zsup,2,eps,"Bord5");
	M.setGroupAtPlan(zinf,2,eps,"Bord6");

    /* set the boundary values for each boundary */
	double T1=0;
	double T2=0;
	double T3=0;
	double T4=0;
	double T5=0;
	double T6=0;

	cout<< "Built a regular tetrahedral 3D mesh from a cube mesh with "<< nx<<"x" <<ny<<"x" <<nz<< " cells"<<endl;

    /* Create the problem */
    bool FEComputation=true;
	StationaryDiffusionEquation myProblem(spaceDim,FEComputation);
	myProblem.setMesh(M);

    /* set the boundary conditions */
	myProblem.setDirichletBoundaryCondition("Bord1",T1);
	myProblem.setDirichletBoundaryCondition("Bord2",T2);
	myProblem.setDirichletBoundaryCondition("Bord3",T3);
	myProblem.setDirichletBoundaryCondition("Bord4",T4);
	myProblem.setDirichletBoundaryCondition("Bord5",T5);
	myProblem.setDirichletBoundaryCondition("Bord6",T6);

	/* Set the right hand side function*/
	Field my_RHSfield("RHS_field", NODES, M, 1);
    Node Ni; 
    double x, y, z;
	for(int i=0; i< M.getNumberOfNodes(); i++)
    {
		Ni= M.getNode(i);
		x = Ni.x();
		y = Ni.y();
		z = Ni.z();

		my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);//mettre la fonction definie au second membre de l'edp
	}
	myProblem.setHeatPowerField(my_RHSfield);
	myProblem.setLinearSolver(GMRES,ILU);

    /* name the result file */
	string fileName = "StationaryDiffusion_3DFE_StructuredTetrahedra";
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
		/* The following formulas use the fact that the exact solution is equal the right hand side divided by 3*pi*pi */
		double max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(3*pi*pi);
		double max_sol_num=my_ResultField.max();
		double min_sol_num=my_ResultField.min();
		double erreur_abs=0;
		for(int i=0; i< M.getNumberOfNodes() ; i++)
			if( erreur_abs < abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i]) )
				erreur_abs = abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i]);
		
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
