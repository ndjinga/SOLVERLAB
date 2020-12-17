//============================================================================
// Author      : Michael NDJINGA
// Date        : November 2020
// Description : multiD linear wave system
//============================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Field.hxx"
#include "SparseMatrixPetsc.hxx"
#include "CdmathException.hxx"

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
    
SparseMatrixPetsc computeDivergenceMatrix(Mesh my_mesh, int nbVoisinsMax, double dt)
{
    int nbCells = my_mesh.getNumberOfCells();
    int dim=my_mesh.getMeshDimension();
    int nbComp=dim+1;
    Vector normal(dim);

    SparseMatrixPetsc implMat(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp);

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
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am);
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.));
            }
            else 
            {    
				if( Fk.getGroupName() != "Periodic" && Fk.getGroupName() != "Neumann")//Wall boundary condition unless Periodic/Neumann specified explicitly
                {
					Vector v(dim+1);
                    for(int i=0; i<dim; i++)
                        v[i+1]=normal[i];
                    Matrix idMoinsJacCL=v.tensProduct(v)*2;
                    
                    implMat.addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL);
				}
                else if( Fk.getGroupName() == "Periodic")//Periodic boundary condition
                {
                    int indexFP=my_mesh.getIndexFacePeriodic(indexFace);
                    Face Fp = my_mesh.getFace(indexFP);
                    cellAutre = Fp.getCellsId()[0];
                    
                    implMat.addValue(j*nbComp,cellAutre*nbComp,Am);
                    implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.));
				}
                else if(Fk.getGroupName() != "Neumann")//Nothing to do for Neumann boundary condition
                {
                    cout<< Fk.getGroupName() <<endl;
                    throw CdmathException("computeDivergenceMatrix: Unknown boundary condition name");
				}
            }
        }   
    }     
    return implMat;
}

void WaveSystem(double tmax, int ntmax, double cfl, int output_freq, const Mesh& my_mesh, const string file)
{
	/* Retrieve mesh data */
    int dim=my_mesh.getMeshDimension();
    int nbCells = my_mesh.getNumberOfCells();
    std::string meshName=my_mesh.getName();
    int nbVoisinsMax=my_mesh.getMaxNbNeighbours(CELLS);
    double dx_min=my_mesh.minRatioVolSurf();


    /* Initial conditions */
    cout << "Construction of the initial condition" << endl;
    
    Field pressure_field("Pressure",CELLS,my_mesh,1) ;
    Field velocity_field("Velocity",CELLS,my_mesh,1) ;
    initial_conditions_shock(my_mesh,pressure_field, velocity_field);

    /* iteration vectors */
    Vector  Un(nbCells*(dim+1));
    Vector dUn(nbCells*(dim+1));
    
    for(int k =0; k<nbCells; k++){
        Un[k*(dim+1)+0] =     pressure_field[k];
        Un[k*(dim+1)+1] =rho0*velocity_field[k,0];
        if(dim>=2)
            Un[k + 2*nbCells] = rho0*velocity_field[k,1] ;
        if(dim==3)
            Un[k + 3*nbCells] = rho0*velocity_field[k,2];
		}

    /*
     * MED output of the initial condition at t=0 and iter = 0
     */
    int it=0;
    bool isStationary=false;
    double time=0.;
    double dt = cfl * dx_min / c0;
    
    cout << "Saving the solution at T=" << time << endl;
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+to_string(dim)+"DUpwind"+meshName+"_velocity");
    /* --------------------------------------------- */

    SparseMatrixPetsc divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt);

    /* Time loop */
    cout<< "Starting computation of the linear wave system with an explicit UPWIND scheme" << endl;
    while (it<ntmax && time <= tmax && ! isStationary)
    {
        dUn=divMat*Un;
        Un-=dUn;
        
        time=time+dt;
        it=it+1;
 
         isStationary = dUn.norm()<precision;
        /* Sauvegardes */
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax)
        {
            cout<<"-- Iteration: " << it << ", Time: " << time << ", dt: " << dt<<endl;

            for(int k=0; k<nbCells; k++)
            {
                pressure_field[k]  =Un[k*(dim+1)+0];
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0;
                if(dim>1)
                {
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0;
                    if(dim>2)
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0;
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
}
 
int main(int argc, char *argv[])
{
	cout << "-- Starting the RESOLUTION OF THE 2D WAVE SYSTEM"<<endl;
	cout << "- Numerical scheme : Upwind explicit scheme" << endl;
	cout << "- Boundary conditions : WALL" << endl;

    // Problem data
    double cfl=0.49;
    double tmax=1.;
    double ntmax=3;//20000;
    int freqSortie=100;
    string fileOutPut="SphericalWave";

	Mesh myMesh;
	
	if(argc<2)
	{
		    cout << "- DOMAIN : SQUARE" << endl;
		    cout << "- MESH : CARTESIAN, GENERATED INTERNALLY WITH CDMATH" << endl<< endl;
		    cout << "Construction of a cartesian mesh" << endl;
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
	    cout << "Loading of a mesh named "<<argv[1] << endl;
	    string filename = argv[1];
	    myMesh=Mesh(filename);
	}

	WaveSystem(tmax,ntmax,cfl,freqSortie,myMesh,fileOutPut);
	
    cout << "Simulation complete." << endl;

    return 0;
}
