//============================================================================
// Author      : Michael NDJINGA
// Date        : November 2020
// Description : multiD linear wave system run in parallel with calls to MPI
//               Test used in the MPI version of the Salome platform
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
    double rayon=0.15;
    double xcentre=0.;
    double ycentre=0;
    double zcentre=0;
	
	double x, y, z;
	double val, valX, valY, valZ;
	
    int dim    =my_mesh.getMeshDimension();
    int nbCells=my_mesh.getNumberOfCells();

    for (int j=0 ; j<nbCells ; j++)
    {
        x = my_mesh.getCell(j).x() ;
		if(dim>1)
		{
			y = my_mesh.getCell(j).y() ;
			if(dim==3)
				z = my_mesh.getCell(j).z() ;
		}

        valX=(x-xcentre)*(x-xcentre);
		if(dim==1)
			val=sqrt(valX);
		else if(dim==2)
		{
			valY=(y-ycentre)*(y-ycentre);
			val=sqrt(valX+valY);		
		}
		else if(dim==3)
		{
			valY=(y-ycentre)*(y-ycentre);
			valZ=(z-zcentre)*(z-zcentre);
			val=sqrt(valX+valY+valZ);		
		}
		
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
}

void WaveSystem(double tmax, int ntmax, double cfl, int output_freq, const Mesh& my_mesh, const string file, int rank, int size, string resultDirectory)
{
	/* Time iteration variables */
    int it=0;
    bool isStationary=false;
    double time=0.;
    double dt;
    double norm;
	    
	/* PETSc variables */
	int globalNbUnknowns;
	int  localNbUnknowns;
	int d_nnz, o_nnz;
	Vec Un, dUn, Un_seq;
    Mat divMat;
	VecScatter scat;

	
	/* Mesh parameters managed only by proc 0 */
    int nbCells;
    int dim;
    int nbComp;	    
	int idx;//Index where to add the vector values
	double value;//value to add in the vector	
    Field pressure_field, velocity_field;
    std::string meshName;

	if(rank == 0)
	    globalNbUnknowns=my_mesh.getNumberOfCells()*(my_mesh.getMeshDimension()+1);//nbCells*nbComp
	
	MPI_Bcast(&globalNbUnknowns, 1, MPI_INT, 0, PETSC_COMM_WORLD);
 
    /* iteration vectors */
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE    ,globalNbUnknowns,&Un);
	VecDuplicate (Un,&dUn);
	if(rank == 0)
		VecCreateSeq(PETSC_COMM_SELF,globalNbUnknowns,&Un_seq);//For saving results on proc 0

	VecScatterCreateToZero(Un,&scat,&Un_seq);

	/* System matrix */
	VecGetLocalSize(Un, &localNbUnknowns);

	if(rank == 0)
	{
		int nbVoisinsMax = my_mesh.getMaxNbNeighbours(CELLS);
	    d_nnz=(nbVoisinsMax+1)*(my_mesh.getMeshDimension()+1);//(nbVoisinsMax+1)*nbComp
	    o_nnz= nbVoisinsMax   *(my_mesh.getMeshDimension()+1);//                 nbComp
	}
	MPI_Bcast(&d_nnz, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&o_nnz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

   	MatCreateAIJ(PETSC_COMM_WORLD,localNbUnknowns,localNbUnknowns,globalNbUnknowns,globalNbUnknowns,d_nnz,NULL,o_nnz,NULL,&divMat);

	if(rank == 0)
		{
		/* Retrieve mesh data */
	    nbCells = my_mesh.getNumberOfCells();
	    dim=my_mesh.getMeshDimension();
	    nbComp=dim+1;	    
	    meshName=my_mesh.getName();
	    double dx_min=my_mesh.minRatioVolSurf();
		dt = cfl * dx_min / c0;
		
	    /* Initial conditions */
	    cout<<"Building the initial condition on processor 0" << endl;
	    
	    pressure_field=Field("Pressure",CELLS,my_mesh,1) ;
	    velocity_field=Field("Velocity",CELLS,my_mesh,dim) ;
	    initial_conditions_shock(my_mesh,pressure_field, velocity_field);
	
	    cout << "Saving the solution at T=" << time <<"  on processor 0"<<endl;
	    pressure_field.setTime(time,it);
	    pressure_field.writeMED(resultDirectory+"/WaveSystem"+to_string(dim)+"DUpwind_"+to_string(size)+"Procs_"+meshName+"_pressure");
	    velocity_field.setTime(time,it);
	    velocity_field.writeMED(resultDirectory+"/WaveSystem"+to_string(dim)+"DUpwind_"+to_string(size)+"Procs_"+meshName+"_velocity");
	    /* --------------------------------------------- */
	

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
	    computeDivergenceMatrix(my_mesh,&divMat,dt);
	}		

	VecAssemblyBegin(Un);
	VecAssemblyEnd(Un);
	
	//VecView(Un, PETSC_VIEWER_STDOUT_WORLD );

	MatAssemblyBegin(divMat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  divMat, MAT_FINAL_ASSEMBLY);

	//MatView(divMat,	PETSC_VIEWER_STDOUT_WORLD );

	MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    /* Time loop */
	PetscPrintf(PETSC_COMM_WORLD,"Starting computation of the linear wave system on all processors : \n\n");
    while (it<ntmax && time <= tmax && !isStationary)
    {
		MatMult(divMat,Un,dUn);
        VecAXPY(Un,-1,dUn);
        
        time=time+dt;
        it=it+1;
 
		VecNorm(dUn,NORM_2,&norm);
        isStationary = norm<precision;
        /* Sauvegardes */
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax )
        {
			PetscPrintf(PETSC_COMM_WORLD,"-- Iteration: %d, Time: %f, dt: %f, saving results on processor 0 \n", it, time, dt);
			VecScatterBegin(scat,Un,Un_seq,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd(  scat,Un,Un_seq,INSERT_VALUES,SCATTER_FORWARD);

			if(rank == 0)
			{
	            for(int k=0; k<nbCells; k++)
	            {
					idx = k*(dim+1)+0;
					VecGetValues(Un_seq,1,&idx,&value);
	                pressure_field[k]  =value;
					for(int idim =0; idim<dim; idim++)
					{
						idx = k*nbComp+1+idim;
						VecGetValues(Un_seq,1,&idx,&value);
						velocity_field[k,idim] = value/rho0;
					}
				}
	            pressure_field.setTime(time,it);
	            pressure_field.writeMED(resultDirectory+"/WaveSystem"+to_string(dim)+"DUpwind_"+to_string(size)+"Procs_"+meshName+"_pressure",false);
	            velocity_field.setTime(time,it);
	            velocity_field.writeMED(resultDirectory+"/WaveSystem"+to_string(dim)+"DUpwind_"+to_string(size)+"Procs_"+meshName+"_velocity",false);
			}
		}
    }

	PetscPrintf(PETSC_COMM_WORLD,"\n End of calculations at iteration: %d, and time: %f\n", it, time);
    if(it>=ntmax)
        PetscPrintf(PETSC_COMM_WORLD, "Nombre de pas de temps maximum ntmax= %d atteint\n", ntmax);
    else if(isStationary)
        PetscPrintf(PETSC_COMM_WORLD, "Régime stationnaire atteint au pas de temps %d, t= %f\n", it, time);       
    else
        PetscPrintf(PETSC_COMM_WORLD, "Temps maximum tmax= %f atteint\n", tmax);

	VecDestroy(&Un);
	VecDestroy(&Un_seq);
	VecDestroy(&dUn);
	MatDestroy(&divMat);
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
    double tmax=1.;
    int ntmax=2000;
    int freqSortie=10;
    string fileOutPut="SphericalWave";
    Mesh myMesh;
	string resultDirectory="./";
	
	if(size>1)
		PetscPrintf(PETSC_COMM_WORLD,"---- More than one processor detected : running a parallel simulation ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Limited parallelism : input and output remain sequential ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Only the matrix-vector products are done in parallel ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"---- Processor 0 is in charge of building the mesh, saving the results, filling and then distributing the matrix to other processors.\n\n");
		
	if(rank == 0)
	{
	    cout << "-- Starting the RESOLUTION OF THE 2D WAVE SYSTEM on "<< size <<" processors"<<endl;
	    cout << "- Numerical scheme : Upwind explicit scheme" << endl;
	    cout << "- Boundary conditions : WALL" << endl;
	
		/* Read or create mesh */
		if(argc<2)
		{
		    cout << "- DOMAIN : SQUARE" << endl;
		    cout << "- MESH : CARTESIAN, GENERATED INTERNALLY WITH CDMATH" << endl<< endl;
		    cout << "Construction of a cartesian mesh on processor 0" << endl;
		    double xinf=-0.5;
		    double xsup= 0.5;
		    double yinf=-0.5;
		    double ysup= 0.5;
		    int nx=50;
		    int ny=50;
		    myMesh=Mesh(xinf,xsup,nx,yinf,ysup,ny);
		    
		    double eps=precision;
		    myMesh.setGroupAtPlan(xsup,0,eps,"RightEdge");
		    myMesh.setGroupAtPlan(xinf,0,eps,"LeftEdge");
		    myMesh.setGroupAtPlan(yinf,1,eps,"BottomEdge");
		    myMesh.setGroupAtPlan(ysup,1,eps,"TopEdge");
		}
		else
		{
		    cout << "- MESH:  GENERATED EXTERNALLY WITH SALOME" << endl<< endl;
		    cout << "Loading of a mesh named "<<argv[1]<<" on processor 0" << endl;
		    string filename = argv[1];
		    myMesh=Mesh(filename);
		}

		/*Detect directory where to same results*/
		if(argc>2)
		    resultDirectory = argv[2];
	}
	
    double cfl=1./myMesh.getSpaceDimension();
	WaveSystem(tmax,ntmax,cfl,freqSortie,myMesh,fileOutPut, rank, size, resultDirectory);

	if(rank == 0)
		cout << "Simulation complete." << endl;

	PetscFinalize();
    return 0;
}
