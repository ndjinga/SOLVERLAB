//============================================================================
// Author      : Michael NDJINGA
// Date        : November 2020
// Description : 2D linear wave system
//============================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Field.hxx"
#include "CdmathException.hxx"

#include <petscksp.h>

using namespace std;

double p0  =155e5;   //reference pressure in a pressurised nuclear vessel
double c0  =700.;    //reference sound speed for water at 155 bars, 600K
double rho0=p0/c0*c0;//reference density
double precision=1e-5;

void initial_conditions_shock(Mesh my_mesh,Field& pressure_field,Field& velocity_field)
{
    double rayon=0.35;
    double xcentre=0.;
    double ycentre=0;

    int dim    =my_mesh.getMeshDimension();
    int nbCells=my_mesh.getNumberOfCells();

    for (int j=0 ; j<nbCells ; j++)
    {
        double x = my_mesh.getCell(j).x() ;
        double y = my_mesh.getCell(j).y() ;
        double valX=(x-xcentre)*(x-xcentre);
        double valY=(y-ycentre)*(y-ycentre);
        double val=sqrt(valX+valY);

		for(int idim=0; idim<dim; idim++)
			velocity_field[j,idim]=0;
			
        if (val<rayon)
            pressure_field(j) = 155e5;
        else
            pressure_field(j) = 70e5;
    }
}

void addValue( int i, int j, Matrix M, Mat * mat )
{
    int I,J;
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
        {
            I=i+k;
            J=j+l;
            MatSetValues( *mat,1, &I, 1, &J, &M(k,l), ADD_VALUES);
        }
}

Matrix jacobianMatrices(Vector normal, double coeff)
{
    int dim=normal.size();
    Matrix    A(dim+1,dim+1);
    Matrix absA(dim+1,dim+1);

    absA(0,0)=c0*coeff;
    for(int i=0 ; i<dim; i++)
    {
        A(i+1,0)=      normal[i]*coeff;
        A(0,i+1)=c0*c0*normal[i]*coeff;
        for( int j=0 ; j<dim; j++)
            absA(i+1,j+1)=c0*normal[i]*normal[j]*coeff;
    }
    return (A - absA)*(1./2);
}
    
void computeDivergenceMatrix(Mesh my_mesh, Mat * implMat, double dt)
{
    int nbCells = my_mesh.getNumberOfCells();
    int dim=my_mesh.getMeshDimension();
    int nbComp=dim+1;
    Vector normal(dim);

    Matrix idMoinsJacCL(nbComp);
    
    for(int j=0; j<nbCells; j++)//On parcourt les cellules
    {
        Cell Cj = my_mesh.getCell(j);
        int nbFaces = Cj.getNumberOfFaces();

        for(int k=0; k<nbFaces; k++)
        {
            int indexFace = Cj.getFacesId()[k];
            Face Fk = my_mesh.getFace(indexFace);
            for(int i =0; i<dim ; i++)
                normal[i] = Cj.getNormalVector(k, i);//normale sortante

            Matrix Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure());

            int cellAutre =-1;
            if ( not Fk.isBorder())
            {
                /* hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if (Fk.getCellsId()[0] == j) 
                    // hypothese verifiée 
                    cellAutre = Fk.getCellsId()[1];
                else if(Fk.getCellsId()[1] == j) 
                    // hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                else
                    throw CdmathException("computeDivergenceMatrix: problem with mesh, unknown cell number");
                    
                addValue(j*nbComp,cellAutre*nbComp,Am      ,implMat);
                addValue(j*nbComp,        j*nbComp,Am*(-1.),implMat);
            }
            else 
            {    
				if( Fk.getGroupName() != "Periodic" && Fk.getGroupName() != "Neumann")//Wall boundary condition unless Periodic/Neumann specified explicitly
                {
					Vector v(dim+1);
                    for(int i=0; i<dim; i++)
                        v[i+1]=normal[i];
                    Matrix idMoinsJacCL=v.tensProduct(v)*2;
                    
                    addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL,implMat);
				}
                else if( Fk.getGroupName() == "Periodic")//Periodic boundary condition
                {
                    int indexFP=my_mesh.getIndexFacePeriodic(indexFace);
                    Face Fp = my_mesh.getFace(indexFP);
                    cellAutre = Fp.getCellsId()[0];
                    
                    addValue(j*nbComp,cellAutre*nbComp,Am      ,implMat);
                    addValue(j*nbComp,        j*nbComp,Am*(-1.),implMat);
				}
                else if(Fk.getGroupName() != "Neumann")//Nothing to do for Neumann boundary condition
                {
                    cout<< Fk.getGroupName() <<endl;
                    throw CdmathException("computeDivergenceMatrix: Unknown boundary condition name");
				}
            }
        }   
    }     
	MatAssemblyBegin(*implMat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  *implMat, MAT_FINAL_ASSEMBLY);
}

void WaveSystem2D(double tmax, double ntmax, double cfl, int output_freq, const Mesh& my_mesh, const string file, int rank, int size)
{
	int globalNbUnknowns;
	
	if(rank != 0)
	{
		MPI_Bcast(&globalNbUnknowns, 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout<<"process "<< rank << " just received globalNbUnknowns= "<< globalNbUnknowns<<endl;
	}
	else
	{
		/* Retrieve mesh data */
	    int nbCells = my_mesh.getNumberOfCells();
	    int dim=my_mesh.getMeshDimension();
	    int nbComp=dim+1;
	    globalNbUnknowns=nbCells*nbComp;
		int buffer[1];
		MPI_Bcast(&globalNbUnknowns, 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout<<"process "<< rank << " just sent globalNbUnknowns= "<< globalNbUnknowns<<endl;
	    
	    std::string meshName=my_mesh.getName();
	    int nbVoisinsMax=my_mesh.getMaxNbNeighbours(CELLS);
	    double dx_min=my_mesh.minRatioVolSurf();
	
		/* time iteration variables */
	    int it=0;
	    bool isStationary=false;
	    double time=0.;
	    double dt = cfl * dx_min / c0;
	    double norm;
	    
	    /* Initial conditions */
	    cout << "Construction of the initial condition …" << endl;
	    
	    Field pressure_field("Pressure",CELLS,my_mesh,1) ;
	    Field velocity_field("Velocity",CELLS,my_mesh,1) ;
	    initial_conditions_shock(my_mesh,pressure_field, velocity_field);
	
	    cout << "Saving the solution at T=" << time << "…" << endl;
	    pressure_field.setTime(time,it);
	    pressure_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_pressure");
	    velocity_field.setTime(time,it);
	    velocity_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_velocity");
	    /* --------------------------------------------- */
	
	    /* iteration vectors */
		Vec Un, dUn;
		VecCreate(PETSC_COMM_WORLD,&Un);
		VecSetSizes(Un,PETSC_DECIDE,nbCells*nbComp);
		VecSetFromOptions(Un);
		VecDuplicate (Un,&dUn);
		int idx;//Index where to add the block of values
		double value;//value to add in the vector
	
	    for(int k =0; k<nbCells; k++)
	    {
			idx = k*nbComp;
			value=pressure_field[k];//vale to add in the vector
			VecSetValues(Un,1,&idx,&value,INSERT_VALUES);
			for(int idim =0; idim<dim; idim++)
			{
				idx = k*nbComp+1+idim;
				value =rho0*velocity_field[k,idim];
				VecSetValues(Un,1,&idx,&value,INSERT_VALUES);
			}
		}
	
		VecAssemblyBegin(Un);
		VecAssemblyEnd(Un);
	
	    Mat divMat;
	   	MatCreateSeqAIJ(PETSC_COMM_WORLD,nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp,NULL,&divMat);
	    computeDivergenceMatrix(my_mesh,&divMat,dt);

    /* Time loop */
    cout<< "Starting computation of the linear wave system with an explicit UPWIND scheme …" << endl;
    while (it<ntmax && time <= tmax )
    {
		MatMult(divMat,Un,dUn);
        VecAXPY(Un,-1,dUn);
        
        time=time+dt;
        it=it+1;
 
		VecNorm(dUn,NORM_2,&norm);
         isStationary = norm<precision;
        /* Sauvegardes */
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax)
        {
            cout<<"-- Iteration: " << it << ", Time: " << time << ", dt: " << dt<<endl;

            for(int k=0; k<nbCells; k++)
            {
				idx = k*(dim+1)+0;
				VecGetValues(Un,1,&idx,&value);
                pressure_field[k]  =value;
				for(int idim =0; idim<dim; idim++)
				{
					idx = k*nbComp+1+idim;
					VecGetValues(Un,1,&idx,&value);
					velocity_field[k,idim] = value/rho0;
				}
			}
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_pressure",false);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_velocity",false);
		}
    }
    cout<<"End of calculation -- Iteration: " << it << ", Time: "<< time<< ", dt: " << dt<<endl;

    if(it>=ntmax)
        cout<< "Nombre de pas de temps maximum ntmax= "<< ntmax<< " atteint"<<endl;
    else if(isStationary)
        cout<< "Régime stationnaire atteint au pas de temps "<< it<< ", t= "<< time<<endl;       
    else
        cout<< "Temps maximum Tmax= "<< tmax<< " atteint"<<endl;

	MatDestroy(&divMat);
	}
}
 
int main(int argc, char *argv[])
{
	/* PETSc initialisation */
	PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
	PetscMPIInt    size;        /* size of communicator */
	PetscMPIInt    rank;        /* processor rank */
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	
    // Problem data
    double cfl=0.49;
    double tmax=1.;
    double ntmax=3;//20000;
    int freqSortie=100;
    string fileOutPut="SphericalWave";
    Mesh myMesh;

	if(size>1)
		PetscPrintf(PETSC_COMM_WORLD,"---- More than one processor detected, running parallel simulation -----\n");
		
	if(rank == 0)
	{
	    cout << "RESOLUTION OF THE 2D WAVE SYSTEM: Upwind explicit scheme" << endl;
	    cout << "- WALL BC" << endl;
	
		if(argc<2)
		{
		    cout << "- DOMAIN: SQUARE" << endl;
		    cout << "- MESH: CARTESIAN, GENERATED INTERNALLY WITH CDMATH" << endl;
		    cout << "Construction of a cartesian mesh …" << endl;
		    double xinf=0.0;
		    double xsup=1.0;
		    double yinf=0.0;
		    double ysup=1.0;
		    int nx=50;
		    int ny=50;
		    myMesh=Mesh(xinf,xsup,nx,yinf,ysup,ny);
		    
		    double eps=1.E-10;
		    myMesh.setGroupAtPlan(xsup,0,eps,"RightEdge");
		    myMesh.setGroupAtPlan(xinf,0,eps,"LeftEdge");
		    myMesh.setGroupAtPlan(yinf,1,eps,"BottomEdge");
		    myMesh.setGroupAtPlan(ysup,1,eps,"TopEdge");
		}
		else
		{
		    cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl;
		    cout << "Loading of a mesh …" << endl;
		    string filename = argv[1];
		    myMesh=Mesh(filename);
		}
	}
	WaveSystem2D(tmax,ntmax,cfl,freqSortie,myMesh,fileOutPut, rank, size);

	if(rank == 0)
		cout << "Simulation complete." << endl;

	PetscFinalize();
    return 0;
}