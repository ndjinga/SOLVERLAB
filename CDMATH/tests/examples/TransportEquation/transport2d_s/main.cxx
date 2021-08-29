//============================================================================
// Author      : Anouar MEKKAS
// Version     :
// Description : 2D linear transport equation on cartesian grid
//============================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Field.hxx"

using namespace std;


void initial_conditions(Field& yField)
{
    double rayon=0.15;
    double xcentre=0.25;
    double ycentre=0.25;
    Mesh myMesh=yField.getMesh();
    int nbCells=myMesh.getNumberOfCells();
    for (int j=0 ; j<nbCells ; j++)
    {
        double x = myMesh.getCell(j).x() ;
        double y = myMesh.getCell(j).y() ;
        double valX=(x-xcentre)*(x-xcentre);
        double valY=(y-ycentre)*(y-ycentre);
        double val=sqrt(valX+valY);
        if (val<rayon)
            yField(j) = 1.0;
        else
            yField(j) = 0.0;
    }
}

void sigma_flux(double VitesseX, double VitesseY, double cfl, const Field& yField, const IntTab indexFacesPerio, double& dt, Field& SumFlux)
{
    /* Fluxes calculation */
    Mesh myMesh=yField.getMesh();
    int nbCells=myMesh.getNumberOfCells();
    double normU=sqrt(VitesseX*VitesseX+VitesseY*VitesseY);
    for (int j=0 ; j<nbCells ; j++)
    {
        Cell Cj=myMesh.getCell(j);
        int nbFace=Cj.getNumberOfFaces();
        double SumF=0.0;
        double minlengthFk=1.E30;

        int cellCourante,cellAutre;
        for (int k=0 ; k<nbFace ; k++)
        {
            int indexFace=Cj.getFacesId()[k];
            Face Fk=myMesh.getFace(indexFace);
            double NormalX=Cj.getNormalVector(k,0);
            double NormalY=Cj.getNormalVector(k,1);
            double LengthFk = Fk.getMeasure();
            double UN=VitesseX*NormalX+VitesseY*NormalY;

            minlengthFk=min(minlengthFk,LengthFk/fabs(UN));
            minlengthFk=min(minlengthFk,LengthFk/fabs(VitesseX));
            minlengthFk=min(minlengthFk,LengthFk/fabs(VitesseY));

            double conc=0.0;
            cellCourante=j;
            cellAutre=-1;

            if (!Fk.isBorder())
            {
                int indexC1=Fk.getCellsId()[0];
                int indexC2=Fk.getCellsId()[1];
                /* Hypothesis: the cell of index indexC1 is the current cell index j */
                if ( indexC1 == j )
                {
                    /* Hypothesis verified */
                    cellCourante=indexC1;
                    cellAutre=indexC2;
                } else if ( indexC2 == j )
                {
                    /* Hypothesis not verified */
                    cellCourante=indexC2;
                    cellAutre=indexC1;
                }
                // Define left and right cells with the product velocity * outward normal
                // If u*n>0, then nothing to do, else invert left and right
                if (UN>1.E-15)
                    conc=yField(cellCourante);
                else
                    conc=yField(cellAutre);
            }else
            {
                /* Homogeneous Neumann boundary conditions */
                if (Fk.getGroupName().compare("LeftEdge")==0 || Fk.getGroupName().compare("RightEdge")==0)
                {
                    if (UN>1.E-15)
                        conc=yField(cellCourante);
                    else
                        conc=0.0;
                }
                /* Periodic boundary conditions */
                if (Fk.getGroupName().compare("BottomEdge")==0 || Fk.getGroupName().compare("TopEdge")==0)
                {
                        int indexFP=indexFacesPerio[indexFace];
                        /* une autre manière de recuperer l'index de la face periodique */
                        //int indexFP=myMesh.getIndexFacePeriodic(indexFace);
                        Face Fp=myMesh.getFace(indexFP);
                        int indexCp=Fp.getCellsId()[0];
                        if (UN>1.E-15)
                            conc=yField(cellCourante);
                        else
                            conc=yField(indexCp);
                }
            }
            SumF=SumF+UN*LengthFk*conc;
          }
        dt=cfl*minlengthFk/normU;
        SumFlux(j)=dt*SumF/Cj.getMeasure();
    }
}

void EquationTransport2D(double tmax, double VitesseX, double VitesseY, double cfl, int freqSortie, const Mesh& myMesh, const string file)
{
    /* Initial conditions */
    cout << "Construction of the initial condition …" << endl;
    Field yField("Y field",CELLS,myMesh,1) ;
    initial_conditions(yField);

    /*
     * MED output of the initial condition at t=0 and iter = 0
     */
    int iter=0;
    double time=0.;
    cout << "Saving the solution at T=" << time << "…" << endl;
    yField.setTime(time,iter);
    yField.writeMED(file);
    yField.writeVTK(file);
    yField.writeCSV(file);
    /* --------------------------------------------- */

    /* Time loop */
    cout << "Resolution of the transport equation with an UPWIND scheme …" << endl;
    int ntmax=3;
    double dt;
    IntTab indexFacesPerio=myMesh.getIndexFacePeriodic();
    while (iter<ntmax && time <= tmax )
    {
        Field SumFlux("Sum Flux",CELLS,myMesh,1) ;
        sigma_flux(VitesseX,VitesseY,cfl,yField,indexFacesPerio,dt,SumFlux);
        cout << "-- Iter: " << iter << ", Time: " << time << ", dt: " << dt << endl;

        /* Advancing time step */
        yField-=SumFlux;

        time+=dt;
        iter+=1;
        // Ouput every freq iterations
        if (iter%freqSortie==0)
        {
            yField.setTime(time,iter);
            yField.writeMED(file,false);
            yField.writeVTK(file,false);
            yField.writeCSV(file);
        }
    }
}

int main()
{
    cout << "RESOLUTION OF 2D TRANSPORT EQUATION:" << endl;
    cout << "- DOMAIN: SQUARE" << endl;
    cout << "- MESH: CARTESIAN, GENERATED INTERNALLY WITH CDMATH" << endl;
    cout << "- PERIODIC BC ON TOP AND BOTTOM" << endl;
    cout << "- HOMOGENEOUS NEUMANN BC ON LEFT AND RIGHT" << endl;

    // Problem data
    double cfl=0.4;
    double VitesseX=1.0;
    double VitesseY=1.0;
    double tmax=1.;
    int freqSortie=10;

    cout << "Construction of a cartesian mesh …" << endl;
    double xinf=0.0;
    double xsup=1.0;
    double yinf=0.0;
    double ysup=1.0;
    int nx=100;
    int ny=100;
    Mesh myMesh(xinf,xsup,nx,yinf,ysup,ny);
    double eps=1.E-10;
    myMesh.setGroupAtPlan(xsup,0,eps,"RightEdge");
    myMesh.setGroupAtPlan(xinf,0,eps,"LeftEdge");
    myMesh.setGroupAtPlan(yinf,1,eps,"BottomEdge");
    myMesh.setGroupAtPlan(ysup,1,eps,"TopEdge");
    string fileOutPutCart="Exercice1";
    EquationTransport2D(tmax,VitesseX,VitesseY,cfl,freqSortie,myMesh,fileOutPutCart);
    cout << "CDMATH calculation done." << endl;
    return 0;
}
