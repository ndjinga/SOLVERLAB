/*
 * meshtests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMATH
 */

#include "MeshTests.hxx"
#include "MEDLoader.hxx"
#include "Face.hxx"
#include "Cell.hxx"
#include "Node.hxx"

#include <string>
#include <cmath>

using namespace MEDCoupling;
using namespace std;

void
MeshTests::testNormals(Mesh mesh)
{
	double eps=1.E-10;
    int dim     = mesh.getSpaceDimension();
    int nbCells = mesh.getNumberOfCells();
    int nbFaces;//local number of faces around a cell
    Face Fk;
    Cell Cj;
    int indexFace;
    
    vector<double> test(dim,0);
    double norm;
    
    for(int j=0; j<nbCells; j++)
    {
        Cj = mesh.getCell(j);
        nbFaces = Cj.getNumberOfFaces();

        for(int i=0; i<dim; i++)
            test[i]=0;

        for(int k=0; k<nbFaces; k++)
        {
            indexFace = Cj.getFacesId()[k];
            Fk = mesh.getFace(indexFace);
            
            norm=0;
            for(int i=0; i<dim; i++)
            {
                test[i] += Cj.getNormalVector(k, i)*Fk.getMeasure();
                norm += Cj.getNormalVector(k, i) * Cj.getNormalVector(k, i);
            }      
            CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., norm, eps );
        }
        
        for(int i=0; i<dim; i++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL( 0., test[i], eps  );
    }
}
//----------------------------------------------------------------------
void
MeshTests::testClassMesh( void )
//----------------------------------------------------------------------
{
	double eps=1.E-10;

	// Testing Mesh(xmin, xMax, nx)
	Mesh M1(0.0,4.0,4);
	CPPUNIT_ASSERT_EQUAL( 1, M1.getSpaceDimension() );
	CPPUNIT_ASSERT_EQUAL( 5, M1.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 4, M1.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( 5, M1.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 4, M1.getNumberOfEdges() );
	CPPUNIT_ASSERT_EQUAL( 0., M1.getFace(0).x() );
	CPPUNIT_ASSERT_EQUAL( 0., M1.getNode(0).x() );
	CPPUNIT_ASSERT_EQUAL( 1., M1.getFace(1).x() );
	CPPUNIT_ASSERT_EQUAL( 1., M1.getNode(1).x() );
	CPPUNIT_ASSERT_EQUAL( 2., M1.getFace(2).x() );
	CPPUNIT_ASSERT_EQUAL( 2., M1.getNode(2).x() );
	CPPUNIT_ASSERT_EQUAL( 3., M1.getFace(3).x() );
	CPPUNIT_ASSERT_EQUAL( 3., M1.getNode(3).x() );
	CPPUNIT_ASSERT_EQUAL( 4., M1.getFace(4).x() );
	CPPUNIT_ASSERT_EQUAL( 4., M1.getNode(4).x() );
	double x11=M1.getCells()[1].x();
	double y11=M1.getCells()[1].y();
	CPPUNIT_ASSERT_EQUAL( x11, 1.5 );
	CPPUNIT_ASSERT_EQUAL( y11, 0.0 );
	M1.setGroupAtFaceByCoords(0.,0.,0.,1.E-14,"LeftEdge") ;
	M1.setGroupAtFaceByCoords(4.,0.,0.,1.E-14,"RightEdge") ;
	CPPUNIT_ASSERT(M1.getFace(0).isBorder()==true);
	CPPUNIT_ASSERT(M1.getFace(1).isBorder()==false);
	CPPUNIT_ASSERT(M1.getFace(2).isBorder()==false);
	CPPUNIT_ASSERT(M1.getFace(3).isBorder()==false);
	CPPUNIT_ASSERT(M1.getFace(4).isBorder()==true);
    CPPUNIT_ASSERT(M1.getNameOfFaceGroups().size() == 3);//There is a default group named "Boundary" that is created by the mesh class
	CPPUNIT_ASSERT(M1.getNameOfFaceGroups()[0].compare("RightEdge")==0);
	CPPUNIT_ASSERT(M1.getNameOfFaceGroups()[1].compare("LeftEdge")==0);
	CPPUNIT_ASSERT(M1.getNameOfFaceGroups()[2].compare("Boundary")==0);
    double dx1=M1.minRatioVolSurf();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dx1, 1., eps  );

    cout<<endl<<"Test mesh M1 normals"<<endl;
    testNormals(M1);

	// Testing Mesh(xmin, xMax, nx, ymin, yMax, ny)
	double xmin=0.0;
	double xmax=4.0;
	double ymin=0.0;
	double ymax=4.0;
	Mesh M2(xmin,xmax,4,ymin,ymax,4);
	CPPUNIT_ASSERT_EQUAL( 4, M2.getNx() );
	CPPUNIT_ASSERT_EQUAL( 4, M2.getNy() );
	CPPUNIT_ASSERT_EQUAL( 2, M2.getSpaceDimension() );
	CPPUNIT_ASSERT_EQUAL( 25, M2.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 16, M2.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( 40, M2.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 40, M2.getNumberOfEdges() );
	CPPUNIT_ASSERT(M2.isQuadrangular());
	int nbCellsM2 = M2.getNumberOfCells();
	double areaM2=0;
	for(int i=0; i<nbCellsM2; i++)
		areaM2+=M2.getCell(i).getMeasure();
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 16., areaM2, eps );
	double x1=M2.getCells()[4].x();
	double y1=M2.getCells()[4].y();
	CPPUNIT_ASSERT_EQUAL( x1, 0.5 );
	CPPUNIT_ASSERT_EQUAL( y1, 1.5 );

	double x2=M2.getNodes()[24].x();
	double y2=M2.getNodes()[24].y();
	CPPUNIT_ASSERT_EQUAL( x2, 4. );
	CPPUNIT_ASSERT_EQUAL( y2, 4. );

	M2.setGroupAtPlan(xmax,0,eps,"RightEdge");
	M2.setGroupAtPlan(xmin,0,eps,"LeftEdge");
	M2.setGroupAtPlan(ymin,1,eps,"BottomEdge");
	M2.setGroupAtPlan(ymax,1,eps,"TopEdge");
	CPPUNIT_ASSERT_EQUAL( 5, int(M2.getNameOfFaceGroups().size()) );//There is a default group named "Boundary" that is created by the mesh class
    CPPUNIT_ASSERT(M2.getNameOfFaceGroups().size() == 5);
	CPPUNIT_ASSERT(M2.getNameOfFaceGroups()[1].compare("BottomEdge")==0);
	int nbFaces=M2.getNumberOfFaces();
    M2.setPeriodicFaces();
	std::map<int,int> indexFaces=M2.getIndexFacePeriodic();
	for (int i=0;i<nbFaces;i++)
	{
		double x=M2.getFaces()[i].x();
		double y=M2.getFaces()[i].y();
		if (y==0. && x==0.5)
		{
			int indexFace=M2.getIndexFacePeriodic(i);
			double xi=M2.getFace(indexFace).x();
			double yi=M2.getFace(indexFace).y();
			CPPUNIT_ASSERT_EQUAL( xi, x );
			CPPUNIT_ASSERT_EQUAL( yi, ymax );
			CPPUNIT_ASSERT_EQUAL( true, M2.getFace(indexFace).isBorder() );
			CPPUNIT_ASSERT_EQUAL( indexFace, indexFaces[i] );
		}
	}
    double dx2=M2.minRatioVolSurf();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dx2, 1./4, eps  );

    cout<<"Test mesh M2 normals"<<endl;
    testNormals(M2);

    // Testing 2D simplexization (regular triangle mesh)
    int splittingPolicy =0;
	Mesh M2Triangle(xmin,xmax,4,ymin,ymax,4,splittingPolicy);
	CPPUNIT_ASSERT_EQUAL( 4, M2Triangle.getNx() );
	CPPUNIT_ASSERT_EQUAL( 4, M2Triangle.getNy() );
	CPPUNIT_ASSERT_EQUAL( 2, M2Triangle.getSpaceDimension() );
	CPPUNIT_ASSERT_EQUAL( 25, M2Triangle.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 32, M2Triangle.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( 40+16, M2Triangle.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 40+16, M2Triangle.getNumberOfEdges() );
	CPPUNIT_ASSERT(M2Triangle.isTriangular());
	int nbCellsM2Triangle = M2Triangle.getNumberOfCells();
	double areaM2Triangle=0;
	for(int i=0; i<nbCellsM2Triangle; i++)
		areaM2Triangle+=M2Triangle.getCell(i).getMeasure();
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 16., areaM2Triangle, eps );

	M2Triangle.setGroupAtPlan(xmax,0,eps,"RightEdge");
	M2Triangle.setGroupAtPlan(xmin,0,eps,"LeftEdge");
	M2Triangle.setGroupAtPlan(ymin,1,eps,"BottomEdge");
	M2Triangle.setGroupAtPlan(ymax,1,eps,"TopEdge");
	CPPUNIT_ASSERT_EQUAL( 5, int(M2Triangle.getNameOfFaceGroups().size()) );//There is a default group named "Boundary" that is created by the mesh class
	CPPUNIT_ASSERT(M2Triangle.getNameOfFaceGroups()[1].compare("BottomEdge")==0);
	std::map<int,int> indexFacesTriangle=M2Triangle.getIndexFacePeriodic();

    double dx2Triangle=M2Triangle.minRatioVolSurf();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dx2Triangle, 0.5/3.414, 0.01  );
    
    cout<<"Test mesh M2Triangle normals"<<endl;
    testNormals(M2Triangle);

	// Testing Mesh(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz) (hexaèdres)
	xmin=0.0;
	xmax=1.0;
	ymin=0.0;
	ymax=1.0;
	double zmin=0.0;
	double zmax=1.0;
    Mesh M3(xmin,xmax,4,ymin,ymax,4,zmin,zmax,4);
    CPPUNIT_ASSERT_EQUAL( 3, M3.getSpaceDimension() );
	CPPUNIT_ASSERT_EQUAL( 5*5*5, M3.getNumberOfNodes() );
	CPPUNIT_ASSERT_EQUAL( 4*4*4, M3.getNumberOfCells() );
	CPPUNIT_ASSERT_EQUAL( 5*4*4*3, M3.getNumberOfFaces() );
	CPPUNIT_ASSERT_EQUAL( 5*5*4*3, M3.getNumberOfEdges() );
    CPPUNIT_ASSERT(M3.isHexahedral());
    int nbCellsM3 = M3.getNumberOfCells();
    double volM3=0;
    for(int i=0; i<nbCellsM3; i++)
        volM3+=M3.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., volM3, eps );

	M3.setGroupAtPlan(xmax,0,eps,"RightEdge");
	M3.setGroupAtPlan(xmin,0,eps,"LeftEdge");
	M3.setGroupAtPlan(ymin,1,eps,"BottomEdge");
	M3.setGroupAtPlan(ymax,1,eps,"TopEdge");
	M3.setGroupAtPlan(zmin,2,eps,"DownEdge");
	M3.setGroupAtPlan(zmax,2,eps,"UpEdge");
	CPPUNIT_ASSERT_EQUAL( 7, int(M3.getNameOfFaceGroups().size()) );//There is a default group named "Boundary" that is created by the mesh class
	CPPUNIT_ASSERT(M3.getNameOfFaceGroups()[1].compare("DownEdge")==0);
	nbFaces=M3.getNumberOfFaces();
    M3.setPeriodicFaces();
	indexFaces=M3.getIndexFacePeriodic();
	for (int i=0;i<nbFaces;i++)
	{
		double x=M3.getFaces()[i].x();
		double y=M3.getFaces()[i].y();
		double z=M3.getFaces()[i].z();
		if (z==0. && x==1./8 && y==1./8)
		{
			int indexFace=M3.getIndexFacePeriodic(i);
			double xi=M3.getFace(indexFace).x();
			double yi=M3.getFace(indexFace).y();
			double zi=M3.getFace(indexFace).z();
			CPPUNIT_ASSERT_EQUAL( xi, x );
			CPPUNIT_ASSERT_EQUAL( yi, y );
			CPPUNIT_ASSERT_EQUAL( zi, zmax );
			CPPUNIT_ASSERT_EQUAL( true, M3.getFace(indexFace).isBorder() );
			CPPUNIT_ASSERT_EQUAL( indexFace, indexFaces[i] );
		}

		if (z==0.5 && y==0. && x==1.)
			CPPUNIT_ASSERT_EQUAL( -1, M3.getIndexFacePeriodic(i) );
	}

    double dx3=M3.minRatioVolSurf();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dx3, 0.25/6, eps  );

    cout<<"Test mesh M3 normals"<<endl;
    testNormals(M3);

    // Testing copies
    Mesh Mcopy1(M1);
    CPPUNIT_ASSERT_EQUAL( 1, Mcopy1.getSpaceDimension() );
    CPPUNIT_ASSERT_EQUAL( 5, Mcopy1.getNumberOfNodes() );
    CPPUNIT_ASSERT_EQUAL( 4, Mcopy1.getNumberOfCells() );
    CPPUNIT_ASSERT_EQUAL( 5, Mcopy1.getNumberOfFaces() );
    CPPUNIT_ASSERT_EQUAL( 4, Mcopy1.getNumberOfEdges() );

    Mcopy1=M2;
    CPPUNIT_ASSERT_EQUAL( 2, Mcopy1.getSpaceDimension() );
    CPPUNIT_ASSERT_EQUAL( 25, Mcopy1.getNumberOfNodes() );
    CPPUNIT_ASSERT_EQUAL( 16, Mcopy1.getNumberOfCells() );
    CPPUNIT_ASSERT_EQUAL( 40, Mcopy1.getNumberOfFaces() );
    CPPUNIT_ASSERT_EQUAL( 40, Mcopy1.getNumberOfEdges() );

    Mesh Mcopy2;
    Mcopy2=Mcopy1;
    CPPUNIT_ASSERT_EQUAL( 2, Mcopy2.getSpaceDimension() );
    CPPUNIT_ASSERT_EQUAL( 25, Mcopy2.getNumberOfNodes() );
    CPPUNIT_ASSERT_EQUAL( 16, Mcopy2.getNumberOfCells() );
    CPPUNIT_ASSERT_EQUAL( 40, Mcopy2.getNumberOfFaces() );
    CPPUNIT_ASSERT_EQUAL( 40, Mcopy2.getNumberOfEdges() );


    // Connection with MED
    string fileNameVTK="TestMesh";
    string fileNameMED="TestMesh";

    M2.writeMED(fileNameMED);
    Mesh M22(fileNameMED + ".med");
    CPPUNIT_ASSERT_EQUAL( 2, M22.getSpaceDimension() );
    CPPUNIT_ASSERT_EQUAL( 25, M22.getNumberOfNodes() );
    CPPUNIT_ASSERT_EQUAL( 16, M22.getNumberOfCells() );
    CPPUNIT_ASSERT_EQUAL( 40, M22.getNumberOfFaces() );
    CPPUNIT_ASSERT_EQUAL( 40, M22.getNumberOfEdges() );

    cout<<"Test mesh M22 normals "<<endl;
    testNormals(M22);

    // Testing 3D simplexization (regular tetrahedra mesh)
    splittingPolicy = 0;
    Mesh M3Tetra(xmin,xmax,4,ymin,ymax,4,zmin,zmax,4,splittingPolicy);
    CPPUNIT_ASSERT_EQUAL( 3, M3Tetra.getSpaceDimension() );
    CPPUNIT_ASSERT(M3Tetra.isTetrahedral());
    int nbCellsM3Tetra = M3Tetra.getNumberOfCells();
    double volM3Tetra=0;
    for(int i=0; i<nbCellsM3Tetra; i++)
        volM3Tetra+=M3Tetra.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., volM3Tetra, eps );

	M3Tetra.setGroupAtPlan(xmax,0,eps,"RightEdge");
	M3Tetra.setGroupAtPlan(xmin,0,eps,"LeftEdge");
	M3Tetra.setGroupAtPlan(ymin,1,eps,"BottomEdge");
	M3Tetra.setGroupAtPlan(ymax,1,eps,"TopEdge");
	M3Tetra.setGroupAtPlan(zmin,2,eps,"DownEdge");
	M3Tetra.setGroupAtPlan(zmax,2,eps,"UpEdge");
	CPPUNIT_ASSERT_EQUAL( 7, int(M3Tetra.getNameOfFaceGroups().size()) );//There is a default group named "Boundary" that is created by the mesh class
	CPPUNIT_ASSERT(M3Tetra.getNameOfFaceGroups()[1].compare("DownEdge")==0);
	indexFaces=M3Tetra.getIndexFacePeriodic();

    cout<<"Test mesh M3Tetra normals"<<endl;
    testNormals(M3Tetra);

    //Testing a 2D unstructured mesh (triangles)
    Mesh M23("meshSquare.med");
    CPPUNIT_ASSERT(M23.getNameOfFaceGroups().size() == 5);//There is a default group named "Boundary" that is created by the mesh class;
    CPPUNIT_ASSERT(M23.getNameOfFaceGroups()[3].compare("Bottom")==0);
    CPPUNIT_ASSERT(M23.getNameOfFaceGroups()[2].compare("Left")==0);
    CPPUNIT_ASSERT(M23.getNameOfFaceGroups()[1].compare("Right")==0);
    CPPUNIT_ASSERT(M23.getNameOfFaceGroups()[0].compare("Top")==0);
    CPPUNIT_ASSERT(M23.isTriangular());
    int nbCellsM23 = M23.getNumberOfCells();
    double areaM23=0;
    for(int i=0; i<nbCellsM23; i++)
        areaM23+=M23.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., areaM23 , eps);

    cout<<"Test mesh M23 normals"<<endl;
    testNormals(M23);

    Mcopy2.writeVTK(fileNameVTK);
    Mcopy2.writeMED(fileNameMED);
    Mesh M6(fileNameMED + ".med");
    CPPUNIT_ASSERT_EQUAL( 2, M6.getSpaceDimension() );
    CPPUNIT_ASSERT_EQUAL( 25, M6.getNumberOfNodes() );
    CPPUNIT_ASSERT_EQUAL( 16, M6.getNumberOfCells() );
    CPPUNIT_ASSERT_EQUAL( 40, M6.getNumberOfFaces() );
    CPPUNIT_ASSERT_EQUAL( 40, M6.getNumberOfEdges() );

    /*
    const MEDCouplingMesh* M1MEDMesh = M2.getMEDCouplingMesh();
     */

    //Test of a mesh with spaceDim=3 different from meshDim=2 (triangles)
    Mesh M4("meshSphere.med");
    CPPUNIT_ASSERT(M4.isTriangular());
    int nbCellsM4 = M4.getNumberOfCells();
    double areaM4=0;
    for(int i=0; i<nbCellsM4; i++)
        areaM4+=M4.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 4*3.14, areaM4, 1 );

    cout<<"Test mesh M4 normals"<<endl;
    testNormals(M4);

    //Testing a 3D unstructured mesh (tétraèdres)
    Mesh M5("meshCube.med");
    CPPUNIT_ASSERT(M5.isTetrahedral());
    int nbCellsM5 = M5.getNumberOfCells();
    double volM5=0;
    for(int i=0; i<nbCellsM5; i++)
        volM5+=M5.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., volM5, eps );
    
    cout<<"Test mesh M5 normals"<<endl;
    testNormals(M5);

    //Testing Mesh( std::vector<double> points, std::string meshName )
    int nbCellsM7 = 2*11;
    int nbNodes = nbCellsM7+1;
    vector<double> points (nbNodes);
    xmin=0;
    xmax=1;
    double dx_min = (xmax-xmin)*2/nbCellsM7/3;
    double dx_max = 2*dx_min;
    points[0]=0;
    for(int i=0; i<nbNodes-1; i++)
    {
        if(i%2==0)
            points[i+1] = (dx_min+dx_max)*(i/2) + dx_min;
        else
            points[i+1] = (dx_min+dx_max)*(i/2) + dx_min + dx_max;        
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(points[0],        xmin,eps);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(points[nbNodes-1],xmax,eps);

    Mesh M7(points, "Checkerboard mesh");

    double volM7=0;
    for(int i=0; i<nbCellsM7; i++)
        volM7+=M7.getCell(i).getMeasure();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., volM7, eps );
}
