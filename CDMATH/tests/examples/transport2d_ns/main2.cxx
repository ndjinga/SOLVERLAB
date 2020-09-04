//============================================================================
// Author      : Anouar MEKKAS
// Version     :
// Description : Equation de transport lineaire 2D non structure
//============================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "Mesh.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Field.hxx"

using namespace std;

struct parameters {
    double cfl;
    double VitesseX;
    double VitesseY;
    int freqSortie;
    std::string outFile;
};

void conditions_initiales(Field& yField)
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

void sigma_flux(const parameters & p, const Field& yField, const IntTab indexFacesPerio, double& dt, Field& SumFlux)
{
    /* Calcul des flux */
    Mesh myMesh=yField.getMesh();
    int nbCells=myMesh.getNumberOfCells();
    double normU=sqrt(p.VitesseX*p.VitesseX+p.VitesseY*p.VitesseY);
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
            double UN=p.VitesseX*NormalX+p.VitesseY*NormalY;

            minlengthFk=min(minlengthFk,LengthFk/fabs(UN));
            minlengthFk=min(minlengthFk,LengthFk/fabs(p.VitesseX));
            minlengthFk=min(minlengthFk,LengthFk/fabs(p.VitesseY));

            double conc=0.0;
            cellCourante=j;
            cellAutre=-1;

            if (!Fk.isBorder())
            {
                int indexC1=Fk.getCellsId()[0];
                int indexC2=Fk.getCellsId()[1];
                /* hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if ( indexC1 == j )
                {
                    /* hypothese verifie */
                    cellCourante=indexC1;
                    cellAutre=indexC2;
                } else if ( indexC2 == j )
                {
                    /* hypothese non verifie */
                    cellCourante=indexC2;
                    cellAutre=indexC1;
                }
                // definir la cellule gauche et droite par le prduit vitesse * normale sortante
                // si u*n>0 : rien a faire sinon inverser la gauche et la droite
                if (UN>1.E-15)
                    conc=yField(cellCourante);
                else
                    conc=yField(cellAutre);
            }else
            {
                /* conditions aux limites neumann homogene */
                if (Fk.getGroupName().compare("GAUCHE")==0 || Fk.getGroupName().compare("DROITE")==0)
                {
                    if (UN>1.E-15)
                        conc=yField(cellCourante);
                    else
                        conc=0.0;
                }
                /* conditions aux limites periodiques */
                if (Fk.getGroupName().compare("BAS")==0 || Fk.getGroupName().compare("HAUT")==0)
                {
                        int indexFP=indexFacesPerio[indexFace];
                        /* une autre manière de recuperer l'index de la face periodique */
                        //int indexFP=M.getIndexFacePeriodic(indexFace);
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
        dt=p.cfl*minlengthFk/normU;
        SumFlux(j)=dt*SumF/Cj.getMeasure();
    }
}

void EquationTransport2D(int &iter, double tmin, double & tmax,
                         const parameters & p, Field & yField)
{
  
    /*
     * MED output of the initial condition at t=0 and iter = 0
     */
    double time=tmin;
    const Mesh& myMesh = yField.getMesh();
    
    yField.setTime(time,iter);
    yField.writeMED(p.outFile,true);
    yField.writeVTK(p.outFile);
    yField.writeCSV(p.outFile);
    /* --------------------------------------------- */

    /* Time loop */
    cout << "Resolution of the transport equation with an UPWIND scheme…" << endl;
    int ntmax=3 + iter;
    double dt;
    IntTab indexFacesPerio=myMesh.getIndexFacePeriodic();
    while (iter<ntmax && time <= tmax )
    {
        Field SumFlux("Fluxes sum",CELLS,myMesh,1) ;
        sigma_flux(p,yField,indexFacesPerio,dt,SumFlux);
        cout << "-- Iter: " << iter << ", Time: " << time << ", dt: " << dt << endl;

        /* Advancing time step */
        yField-=SumFlux;

        time+=dt;
        iter+=1;
        // Output every freq iterations
        if (iter % p.freqSortie==0)
        {
            yField.setTime(time,iter);
            yField.writeVTK(p.outFile,false);
            yField.writeCSV(p.outFile);
        }
    }

    cout << "MED post-treatment of the solution at T=" << time << endl;
    yField.setTime(time,iter);
    yField.writeMED(p.outFile,false);
    yField.writeCSV(p.outFile);

    tmax = time;
    cout << "End of EquationTransport2D." << endl;
}

void compute_onepass(double tmin, double tmax, const Mesh & M, const parameters & p)
{
	int iter = 0;

    /* Initial conditions */
    cout << "Construction of initial condition…" << endl;
    Field f("Y field",CELLS,M,1) ;
    conditions_initiales(f);

    EquationTransport2D(iter, tmin, tmax, p, f);
}

void compute_twopass(double tmin, double tmax, const Mesh & M, const parameters & p)
{
	int iter = 0;
	parameters q(p);
	double tstep = 0.5 * (tmin + tmax);

    /* Initial conditions */
    cout << "Construction of initial condition…" << endl;
    Field f("Y field",CELLS,M,1) ;
    conditions_initiales(f);

	q.outFile = p.outFile + "_A";
    EquationTransport2D(iter, tmin, tstep, q, f);

    Field f2(q.outFile, CELLS, "Y field", iter, 0);

	q.outFile = p.outFile + "_B";
    EquationTransport2D(iter, tstep, tmax, q, f2);

}

int main()
{
    cout << "Resolution of the 2D transport equation:" << endl;
    cout << "- DOMAIN: SQUARE" << endl;
    cout << "- MESH: TRIANGULAR, GENERATED WITH SALOME" << endl;
    cout << "- PERIODIC BC ON TOP AND BOTTOM" << endl;
    cout << "- HOMOGENEOUS NEUMANN BC ON LEFT AND RIGHT" << endl;

    // donnees du probleme
    parameters p;
    int iter;

    p.cfl=0.4;
    p.VitesseX=1.0;
    p.VitesseY=1.0;
    p.freqSortie=50;
    p.outFile="res2D";

    cout << "Construction of Cartesian mesh…" << endl;
    Mesh M("../../tests/ressources/meshSquare.med");

    double tmin = 0;
    double tmax = 0.1;
    
    compute_onepass(tmin, tmax, M, p);
    compute_twopass(tmin, tmax, M, p);
    cout << "CDMATH calculation done." << endl;

    return 0;
}
