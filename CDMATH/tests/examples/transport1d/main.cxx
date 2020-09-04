//============================================================================
// Author      : Anouar MEKKAS
// Version     :
// Description : 1D linear transport equation
//============================================================================

#include <iostream>
#include <cmath>

#include "Cell.hxx"
#include "Mesh.hxx"
#include "Field.hxx"

using namespace std;


int main( void )
{
  double a=-5.0;
  double b=5.0;
  int nx=1000;
  int ntmax=3;
  double dx = (b-a)/nx;
  double pi=3.1415927;
  // Transport velocity
  double cfl=0.5;
  double u=3.;
  double dt=cfl*dx/u;

  Mesh myMesh(a,b,nx);

  Field conc("Concentration",CELLS,myMesh,1);

  // Initial conditions
  double sigma=sqrt(0.2);
  for (int i=0 ; i<myMesh.getNumberOfCells() ; i++)
  {
   double x=myMesh.getCell(i).x();
   conc(i) = 0.5/(sigma*sqrt(2*pi))*exp(-0.5*pow((x/sigma),2));
  }

  double time=0.;
  double tmax=3.0;
  int iter=0;

  cout << "MED post-treatment of the solution at T=" << time << "…" << endl;
  string fileOutPut="EqTr1D";
  conc.setTime(time,iter);
  conc.writeMED(fileOutPut);
  conc.writeVTK(fileOutPut);
  conc.writeCSV(fileOutPut);
  int outputFreq=10;

  // Time loop
  while (iter<ntmax && time <= tmax )
  {
   cout << "-- Iter: " << iter << ", Time: " << time << ", dt: " << dt << endl;
   conc(0) = conc(0) -u*dt/dx*(conc(0)-conc(myMesh.getNumberOfCells()-1));
   for (int j=1 ; j<myMesh.getNumberOfCells() ; j++)
   {
    conc(j) = conc(j) -u*dt/dx*(conc(j)-conc(j-1));
   }
   time+=dt;
   iter+=1;
   if (iter%outputFreq==0)
   {
       conc.setTime(time,iter);
       conc.writeMED(fileOutPut,false);
       conc.writeVTK(fileOutPut,false);
       conc.writeCSV(fileOutPut);
   }
  }
  cout << "CDMATH calculation done." << endl;
  return 0;
}


