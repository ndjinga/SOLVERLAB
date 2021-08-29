//============================================================================
// Name        : Thermique 1D Entrée sortie
// Author      : M. Ndjinga
// Version     :
// Copyright   : CEA Saclay 2014
// Description : Modélisation 1D d'un cœur de réacteur
//============================================================================

#include <iostream>
#include <cmath>
//#include<complex>

//#include "Matrix.hxx"
#include "Vector.hxx"
#include "LinearSolver.hxx"
#include "Mesh.hxx"
#include "Field.hxx"
#include "Cell.hxx"
//#include "Face.hxx"
//#include "Node.hxx"
#include "CdmathException.hxx"

using namespace std;


void champ_heaviside(Field& Temp, double milieu, double TempLeft, double TempRight)
{
    Mesh M=Temp.getMesh();
    int nbCells=M.getNumberOfCells();
    double x;
    for (int j=0 ; j<nbCells ; j++)
    {
        x = M.getCell(j).x() ;
        //cout<<"cellule "<< j<<" x= "<< x<<" milieu= "<< milieu <<endl;
        if (x<milieu)
            Temp(j) = TempLeft;
        else
            Temp(j) = TempRight;
    }
}


double sign(double x)
{
    if(x>0)
        return 1.;
    else if(x<0)
            return -1.;
        else
            return 0.;
}


double theta_upwind(double Tim1, double Ti, double Tip1, double Tip2, double VitesseX, double dt, double Tin, double Tip1n)
{
    if(abs(Tip1-Ti)<1.e-5)
        return 1.;

    double denominateur = abs(VitesseX)*(2.-sign(Tip1-Ti)*(sign(Tip2-Tip1)+sign(Ti-Tim1)));

    if( abs(denominateur) <=1.e-5 )
                return 0.5;
    else
    {
        double result = (1-sign(Tip1-Ti)*sign(Ti-Tim1))/denominateur;
        if(result <0.5)
            return 0.5;
        else
            return result;
    }
}


void solveTriDiag(vector<double>& a,vector<double>& b,vector<double>& c,Field &d, int nbCells)
{
    int n= nbCells;

    if(a.size()!=c.size() || a.size()+1 != b.size() || n != (int) b.size() || n != (int) (d.getMesh()).getNumberOfCells())
        throw CdmathException("solveTriDiag(): vectors with inapropriate size");

/*  cout<<"avant solveTri"<<endl;
    for(int i=0;i<nbCells-2;i++)
        cout<< " a= "<< a[i]<<" b= "<< b[i]<<" c= "<< c[i]<< " d= "<< d[i]<<endl;*/

    /* Algorithme de résolution d'un Système Tridiagonal*/
         n--; // since we start from x0 (not x1)
         c[0] /= b[0];
         d[0] /= b[0];
         for (int i = 1; i < n; i++) {
            c[i] /= b[i] - a[i-1]*c[i-1];
            d[i] = (d[i] - a[i-1]*d[i-1]) / (b[i] - a[i-1]*c[i-1]);
            }

            d[n] = (d[n] - a[n-1]*d[n-1]) / (b[n] - a[n-1]*c[n-1]);

            for (int i = n; i-- > 0;) {
                d[i] -= c[i]*d[i+1];
            }

        /*cout<<"apres solveTri"<<endl;
        for(int i=0;i<nbCells-2;i++)
            cout<< " a= "<< a[i]<<" b= "<< b[i]<<" c= "<< c[i]<< " d= "<< d[i]<<endl;*/
}


void solveEntreeImplicit(Field& phi, Field& Temp, int nbCells, double VitesseX, double cfl, double Tentree, double dt) {

    /* Coefficients de la matrice tridiagonalee liee au systeme*/
    /*b=coeff diag, c= coeff au dessus diag, a= coeff sous diag */
    vector<double> a(nbCells-1),b(nbCells),c(nbCells-1);
    double theta, epsilon=1.e-2, erreur=1/epsilon;
    int iterMax=100, k=0;
    Field d, ScndMbre=phi, Tempn=Temp;
    ScndMbre*=dt;
    ScndMbre+=Tempn;//On obtient Tn+dt*phi
    //On décentre le terme source
    for (int i = 1; i < nbCells-1; i++)
        {
        ScndMbre[i]-=cfl*(phi[i]-phi[i-1])*dt/2;
        }
    ScndMbre[0]-=cfl*(phi[0])*dt/2;
        //d[0]-=cfl*(phi[n-1]-phi[n-2]);

    if(nbCells<3)
        throw CdmathException("solveEntreeImplicit(): mesh should have at least three cells");

    while( k<iterMax && erreur>epsilon)
    {
    d= ScndMbre;//contiendra le second membre du système
    /*cout<< " Debut itération d= "<<endl;
    for(int i=0;i<nbCells-1;i++)
            cout<< " , "<< d(i);
    cout<<endl;*/
    cout<<"theta="<<endl;
    /* On traite la première face (entrée) */
    theta=theta_upwind(Tentree,Tentree,Temp(0),Temp(1),VitesseX, dt, Tentree, Tempn(0));
    b[0]=1+ cfl*(2*theta-1);
    c[0]=(1-theta)*cfl;
    d(0)+=theta*cfl*Tentree;

    cout<< theta << " , ";

    /* On traite la deuxième face interne */
    theta=theta_upwind(Tentree,Temp(0),Temp(1),Temp(2),VitesseX, dt, Tempn(0), Tempn(1));
    b[1]=1+ cfl*(2*theta-1);
    c[1]=(1-theta)*cfl;
    a[0]=-theta*cfl;

    cout<< theta << " , ";

    //On traite les faces internes
    for(int i=2; i<nbCells-2; i++)
    {
        theta=theta_upwind(Temp(i-2),Temp(i-1),Temp(i),Temp(i+1),VitesseX, dt, Tempn(i-1), Tempn(i));
        b[i]=1+ cfl*(2*theta-1);
        c[i]=(1-theta)*cfl;
        a[i-1]=-theta*cfl;
        cout<< theta << " , ";
    }

    /* On traite l'avant dernière face interne */
    theta=theta_upwind(Temp(nbCells-3),Temp(nbCells-2),Temp(nbCells-1),Temp(nbCells-1),VitesseX, dt, Tempn(nbCells-2), Tempn(nbCells-1));
    b[nbCells-2]=1+ cfl*(2*theta-1);
    c[nbCells-2]=(1-theta)*cfl;
    a[nbCells-3]=-theta*cfl;

    cout<< theta << " , ";

    /* On traite la dernière face (sortie) */
    theta=theta_upwind(Temp(nbCells-2),Temp(nbCells-1),Temp(nbCells-1),Temp(nbCells-1),VitesseX, dt, Tempn(nbCells-1), Tempn(nbCells));
    b[nbCells-1]=1+ cfl*theta;
    a[nbCells-2]=-theta*cfl;
    //d(nbCells-1)+=-(1-theta_sortie)*cfl*Tempn(nbCells-1);

    cout<< theta << endl;

    cout<<" k= " << k<<endl;
    /*cout<<"avant solveTri"<<endl;
    for(int i=0;i<nbCells-2;i++)
        cout<< " a= "<< a[i]<<" b= "<< b[i]<<" c= "<< c[i]<< " d= "<< d[i]<<endl;*/
    solveTriDiag(a,b,c,d,nbCells);
    /*cout<<"après solveTri"<<endl;
    for(int i=0;i<nbCells-2;i++)
        cout<< " a= "<< a[i]<<" b= "<< b[i]<<" c= "<< c[i]<< " d= "<< d[i]<<endl;*/
    k++;
    erreur = 0;
    for(int i=0;i<nbCells; i++)
        erreur=max(erreur, abs(Temp(i)-d[i])/300);
    cout<< " erreur= "<<erreur<<endl;
    /*cout<< " Fin itération d= "<<endl;
    for(int i=0;i<nbCells-1;i++)
            cout<< " , "<< d(i);
    cout<<endl;*/
    Temp=d;
//  for(int i=0;i<nbCells-1;i++)
//      cout<< " Temp= "<< Temp(i)<<endl;
    }

    if(k>=iterMax)
        throw CdmathException("solveEntreeImplicit: Newton scheme not convergent");
}


void EquationTransport1D_entree(double tmax, double VitesseX, int ntmax, double dt, double cfl,int freqSortie, const Mesh& M, const string file, double milieu, double TempLeft, double TempRight)
{
    /* --------------------------------------------- */
    /* Condition initiale */
    cout << "Construction de la condition initiale ... " << endl;
    Field Temp("Temperature",CELLS,M,1) ;
    champ_heaviside(Temp,milieu,TempLeft, TempRight);

    /*terme source */
    cout << "Construction du terme source ... " << endl;
    Field phi("Flux thermique",CELLS,M,1) ;
    champ_heaviside(phi,milieu,-10, 10);
    /*
     * Sortie MED de la condition initiale et du flux thermique à t=0 et iter = 0
     */
    int iter=0;
    double time=0.;
    int nbCells=M.getNumberOfCells();
    cout << "Post-traitement MED du flux thermique constant en temps"<< " ..." << endl;
    //phi.setTime(time,iter);
    //phi.writeCSV(file);
    cout << "Post-traitement MED de la solution à t=" << time << " ..." << endl;
    Temp.setTime(time,iter);
    //Temp.writeCSV(file);
    Temp.writeVTK(file);
    /* boucle de temps */
    cout << " Resolution de l'equation de transport 1D avec entree par un schema Implicite TVD" << endl;

    while (iter<ntmax && time <= tmax )
    {
        cout << "-- Iter : " << iter << " Time : " << time << " dt : " << dt << endl;

        /* Avancement en temps */
            solveEntreeImplicit(phi, Temp,nbCells, VitesseX, cfl, TempLeft, dt);
        time+=dt;
        iter+=1;
        // sortie visu tous les freq iterations
        if (iter%freqSortie==0)
        {
            Temp.setTime(time,iter);
            Temp.writeVTK(file,false);
            //Temp.writeCSV(file);
        }
    }
}


int Equation_Transport()
{
    cout << "RESOLUTION EQUATION DE TRANSPORT 1D :" << endl;
    cout << "- DOMAINE : Segment 0,10" << endl;
    cout << "- MAILLAGE CARTESIEN : GENERATION INTERNE CDMATH " << endl;

    // Problem data
    double VitesseX=1.0;
    double tmax=100.;
    int freqSortie=1;
    int ntmax=3;
    int nx=50;
    double xinf=0.0;
    double xsup=10.;
    double dx=(xsup-xinf)/nx;
    double cfl=1;
    double dt=cfl*dx/abs(VitesseX);

    double TempLeft=300;//285.;
    double TempRight=300;//315.;

    cout << "Début Construction du maillage Cartesien…" << endl;
    Mesh M(xinf,xsup,nx);
    cout << "Fin Construction du maillage Cartesien…" << endl;

    //  theta=1;
    EquationTransport1D_entree(tmax, VitesseX, ntmax, dt, cfl, freqSortie, M, "TemperatureEntreeImplicitTVD50CellsSourceUpwind", (xsup+xinf)/2, TempLeft, TempRight);
    //EquationTransport1D_entree(tmax, VitesseX, ntmax, dt, cfl, freqSortie, M, "TemperatureEntreeImplicitUpwindCreneau50Cells", (xsup+xinf)/4, TempLeft, TempRight);

    cout << "CDMATH calculation done." << endl;
    return 0;
}


int main() {
/* Test solveur tridiagonal
    int n=6;
    Matrix A2(n,n,n+2*(n-1));

       A2(0,0)=2.;
       vector<double> e(n-1);
       for(int i=0;i<n-1;i++)
           {
           e[i]=-1/(i+1);
           A2(i,i+1)=-1/(i+1);
           A2(i+1,i)=-1.;
           A2(i+1,i+1)=2.;
           }
       Vector Xana2(n);
       for(int i=0;i<n;i++)
            Xana2(i)=i;

       Vector B2=A2*Xana2;

       LinearSolver LS11(A2,B2,500,1.E-10,"GMRES","ILU");
       Vector X11=LS11.solve();
       bool correct=true;
       for (int i=0;i<X11.getNumberOfRows();i++)
        correct=correct && (abs(X11(i)-Xana2(i))<1.E-10);

       if(correct)
           cout<<"OK"<<endl;
       else
           cout<<"KO"<<endl;
       vector<double> a(n-1,-1),b(n,2),c(n-1,-0.5);
       Mesh M(0,1,n);
       Field d("Test",CELLS,M,1);
       int nbCells=M.getNumberOfCells();
       for(int i=0;i<nbCells;i++)
           d[i]=B2(i);
       solveTriDiag(a,b,e,d,nbCells);
       correct=true;
       for (int i=0;i<X11.getNumberOfRows();i++)
        correct=correct && (abs(d(i)-Xana2(i))<1.E-10);

       if(correct)
           cout<<"OK"<<endl;
       else
           cout<<"KO"<<endl;

        Fin Test solveur tridiagonal */

       int ret=Equation_Transport();
       return ret;
}
