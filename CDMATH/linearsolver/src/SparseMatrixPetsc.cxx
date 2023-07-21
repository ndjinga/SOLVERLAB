/*
 * SparseMatrixPetsc.cxx
 *
 *  Created on: 04/11/2017
 *      Author: mndjinga
 */

#include "SparseMatrixPetsc.hxx"
#include "CdmathException.hxx"
#include <petscdraw.h>

using namespace std;

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc()
//----------------------------------------------------------------------
{
	_numberOfColumns=0;
	_numberOfRows=0;
	_numberOfNonZeros=0;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int numberOfRows, int numberOfColumns)
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_isSparseMatrix=true;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqAIJ(MPI_COMM_SELF,_numberOfRows,_numberOfColumns,PETSC_DEFAULT,NULL,&_mat);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( Mat mat )
//----------------------------------------------------------------------
{
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	_isSparseMatrix=true;
	_mat=mat;
	//extract number of row and column
	MatGetSize(mat,&_numberOfRows,&_numberOfColumns);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int numberOfRows, int numberOfColumns, int nnz )
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_numberOfNonZeros=nnz;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqAIJ(MPI_COMM_SELF,_numberOfRows,_numberOfColumns,_numberOfNonZeros,NULL,&_mat);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int blockSize, int numberOfRows, int numberOfColumns, int nnz )
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_numberOfNonZeros=nnz;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqBAIJ(MPI_COMM_SELF,blockSize, _numberOfRows,_numberOfColumns,_numberOfNonZeros,NULL,&_mat);
}

SparseMatrixPetsc::SparseMatrixPetsc(std::string filename, bool hdf5BinaryMode)
{
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	_mat = NULL;
	readPETScMatrixFromFile( filename, hdf5BinaryMode);

	_isSparseMatrix=true;
	//extract number of row and column
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc(const SparseMatrixPetsc& matrix)
//----------------------------------------------------------------------
{
	_isSparseMatrix=matrix.isSparseMatrix();
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&_mat);
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
}

SparseMatrixPetsc
SparseMatrixPetsc::transpose() const
{
	Mat mattranspose;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatTranspose(_mat,MAT_INITIAL_MATRIX, &mattranspose);
	return SparseMatrixPetsc(mattranspose);
}

void
SparseMatrixPetsc::setValue( int i, int j, double value )
{
	MatSetValues(_mat,1, &i, 1, &j, &value, INSERT_VALUES);
}

void
SparseMatrixPetsc::addValue( int i, int j, double value )
{
	MatSetValues(_mat,1, &i, 1, &j, &value, ADD_VALUES);
}

void
SparseMatrixPetsc::setValue( int i, int j, Matrix M  )
{
    int I,J;
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
        {
            I=i+k;
            J=j+l;
            MatSetValues(_mat,1, &I, 1, &J, &M(k,l), INSERT_VALUES);
        }
}

void
SparseMatrixPetsc::addValue( int i, int j, Matrix M  )
{
    int I,J;
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
        {
            I=i+k;
            J=j+l;
            MatSetValues(_mat,1, &I, 1, &J, &M(k,l), ADD_VALUES);
        }
}

void
SparseMatrixPetsc::setValuesBlocked( int i, int j, Matrix M  )
{
    int blockSize;
    MatGetBlockSize(_mat,&blockSize);
    if(blockSize!=M.getNumberOfRows() || blockSize!=M.getNumberOfColumns())
        throw CdmathException("SparseMatrixPetsc::setValuesBlocked : matrix size is different from sparse matrix block structure");
    double petscValues[blockSize*blockSize];
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
            petscValues[k*blockSize+l]=M(k,l);
    MatSetValuesBlocked(_mat,1, &i, 1, &j, petscValues, INSERT_VALUES);
}

void
SparseMatrixPetsc::addValuesBlocked( int i, int j, Matrix M  )
{
    int blockSize;
    MatGetBlockSize(_mat,&blockSize);
    if(blockSize!=M.getNumberOfRows() || blockSize!=M.getNumberOfColumns())
        throw CdmathException("SparseMatrixPetsc::addValuesBlocked : matrix size is different from sparse matrix block structure");
    double petscValues[blockSize*blockSize];
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
            petscValues[k*blockSize+l]=M(k,l);
    MatSetValuesBlocked(_mat,1, &i, 1, &j, petscValues, ADD_VALUES);
}

//----------------------------------------------------------------------
double
SparseMatrixPetsc::operator()( int i, int j ) const
//----------------------------------------------------------------------
{
	double res;
	int idxm=i,idxn=j;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatGetValues(_mat,1,&idxm,1, &idxn,&res);
	return res;
}

//----------------------------------------------------------------------
SparseMatrixPetsc::~SparseMatrixPetsc()
//----------------------------------------------------------------------
{
	if(&_mat != NULL)
		MatDestroy(&_mat);
	//PetscFinalize();
}

Mat
SparseMatrixPetsc::getPetscMatrix() const
{
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	return (_mat);
}

bool
SparseMatrixPetsc::containsPetscMatrix() const
{
	return true;
}
//----------------------------------------------------------------------
const SparseMatrixPetsc&
SparseMatrixPetsc::operator= ( const SparseMatrixPetsc& matrix )
//----------------------------------------------------------------------
{
	_isSparseMatrix=matrix.isSparseMatrix();
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&_mat);
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
	return (*this);
}

SparseMatrixPetsc
operator+ (const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="SparseMatrixPetsc::operator()+(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): number of rows or columns of the matrices is different!";
		throw CdmathException(msg);
	}

	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatDuplicate(mat1, MAT_COPY_VALUES,&mat);
	MatAXPY(mat,1,mat2,DIFFERENT_NONZERO_PATTERN);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator- (const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="SparseMatrixPetsc::operator()-(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): number of rows or columns of the matrices is different!";
		throw CdmathException(msg);
	}
	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatDuplicate(mat1, MAT_COPY_VALUES,&mat);
	MatAXPY(mat,-1,mat2,DIFFERENT_NONZERO_PATTERN);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator* (double value , const SparseMatrixPetsc& matrix )
{
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator* (const SparseMatrixPetsc& matrix, double value )
{
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator/ (const SparseMatrixPetsc& matrix, double value)
{
	if(value==0.)
	{
		string msg="SparseMatrixPetsc SparseMatrixPetsc::operator()/(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): division by zero";
		throw CdmathException(msg);
	}
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, 1/value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator*(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatMatMult(mat1, mat2, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&mat);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator*= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	Mat mat2;
	MatMatMult(_mat, mat1, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&mat2);

	MatDestroy(&_mat);
	_mat=mat2;
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
	return (*this);
}

Vector
SparseMatrixPetsc::operator* (const Vector& vec) const
{
	int numberOfRows=vec.getNumberOfRows();
	Vec X=vectorToVec(vec);
	Vec Y;
	VecDuplicate (X,&Y);
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatMult(_mat,X,Y);
 
	//Clean memory
	VecDestroy(&X);

    Vector result=vecToVector(Y);

	//Clean memory
	VecDestroy(&Y);

	return result;
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator+= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatAXPY(_mat,1,mat1,DIFFERENT_NONZERO_PATTERN);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;

	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator-= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatAXPY(_mat,-1,mat1,DIFFERENT_NONZERO_PATTERN);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;

	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator*= (double value)
{
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatScale(_mat, value);
	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator/= (double value)
{
	if(value==0.)
	{
		string msg="SparseMatrixPetsc SparseMatrixPetsc::operator()/=(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): division by zero";
		throw CdmathException(msg);
	}
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatScale(_mat, 1/value);
	return (*this);
}

void
SparseMatrixPetsc::viewNonZeroStructure(double pause_lenght, std::string matrixName) const 
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	PetscViewer viewer = PETSC_VIEWER_DRAW_WORLD;
	std::string filename="MatrixNonZeroPlot_"+matrixName+".ppm";//To save non zero structure to a picture file
	
   PetscDraw draw;
   PetscViewerDrawGetDraw(viewer, 0, &draw);
   PetscDrawSetPause(draw, pause_lenght);
   PetscDrawSetSave( draw, filename.c_str());

   MatView(_mat, viewer);
}

void
SparseMatrixPetsc::printCoefficients(std::string matrixName) const 
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;
	std::string filename="MatrixCoefficients_"+matrixName+".txt";//To save matrix coefficients to a text file
	
	MatView(_mat, viewer);//Print matrix coefficients on screen
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);//Prepare dumping to file

    MatView(_mat, viewer);
	PetscViewerDestroy(&viewer);
}

void 
SparseMatrixPetsc::saveToFile(string filename, bool binaryMode) const
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	PetscViewer fileViewer;

	if( binaryMode)
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &fileViewer);
	else
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &fileViewer);
     
	MatView(_mat,fileViewer);
	PetscViewerDestroy(&fileViewer);
}

void 
SparseMatrixPetsc::readPETScMatrixFromFile(std::string filename, bool hdf5BinaryMode)
{
	//Create the viewer to read the file
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
	if(hdf5BinaryMode)
		PetscViewerSetType(viewer,PETSCVIEWERHDF5);
	else
		PetscViewerSetType(viewer,PETSCVIEWERBINARY);
	PetscViewerSetFromOptions(viewer);
	PetscViewerFileSetMode(viewer, FILE_MODE_READ);
	PetscViewerFileSetName(viewer,filename.c_str());
	
	//Empty the matrix (delete current content)
	if( _mat )
		MatDestroy(&_mat);
	
	//Create and load the matrix
	MatCreate(PETSC_COMM_WORLD, &_mat);
	MatSetFromOptions(_mat);
	MatLoad( _mat, viewer);	

	PetscViewerDestroy(&viewer);
}

double
SparseMatrixPetsc::getMatrixCoeff(int i, int j) const 
{
	double res;
	int idxm=i,idxn=j;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatGetValues(_mat,1,&idxm,1, &idxn,&res);
	return res;
}

std::vector< double > 
SparseMatrixPetsc::getArray()
{
	int size=_numberOfRows*_numberOfColumns;
	
	vector< double >  result(size);	
	double* values = result.data();
	
	int * idxm = new int[_numberOfRows];
	int * idxn = new int[_numberOfColumns];
    for (int i=0;i<_numberOfRows;i++) 
		idxm[i]=i;
    for (int i=0;i<_numberOfColumns;i++) 
		idxn[i]=i;
	
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatGetValues(_mat,_numberOfRows, idxm,_numberOfColumns, idxn,values);

	return result;
}

void 
SparseMatrixPetsc::diagonalShift(double lambda)
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

    MatShift(_mat, lambda);
}

void 
SparseMatrixPetsc::zeroEntries()
{
    MatZeroEntries(_mat);
}

Vector
SparseMatrixPetsc::vecToVector(const Vec& vec) const
{
	PetscInt numberOfRows;
	VecGetSize(vec,&numberOfRows);
    double * petscValues;
    VecGetArray(vec,&petscValues);
    
    DoubleTab values (numberOfRows,petscValues);
	Vector result(numberOfRows);
    result.setValues(values);

	return result;
}
Vec
SparseMatrixPetsc::vectorToVec(const Vector& myVector) const
{
	int numberOfRows=myVector.getNumberOfRows();
	const double* values = myVector.getValues().getValues();
	
	Vec result;
	VecCreate(PETSC_COMM_WORLD,&result);
	VecSetSizes(result,PETSC_DECIDE,numberOfRows);
	VecSetBlockSize(result,numberOfRows);
	VecSetFromOptions(result);
	int idx=0;//Index where to add the block of values
	VecSetValuesBlocked(result,1,&idx,values,INSERT_VALUES);

	VecAssemblyBegin(result);
	VecAssemblyEnd(result);

	return result;
}

bool SparseMatrixPetsc::isSymmetric(double tol) const
{
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	//Check that the matrix is symmetric
	PetscBool isSymetric;
	MatIsSymmetric(_mat,tol,&isSymetric);

	return isSymetric;
}

int 
SparseMatrixPetsc::computeSpectrum(int nev, double ** valPr, double ** valPi, double ***vecPr, double ***vecPi, EPSWhich which, double tol, EPSType type, bool viewEigenvaluesInXWindows, double pause_lenght, std::string matrixName) const
{
  EPS            eps;         /* eigenproblem solver context */
  PetscReal      error;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       i,maxit,its, nconv;

  SlepcInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);

  MatAssemblyBegin(_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_mat,MAT_FINAL_ASSEMBLY);

  MatCreateVecs(_mat,&xr,NULL);
  MatCreateVecs(_mat,&xi,NULL);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  EPSCreate(PETSC_COMM_WORLD,&eps);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  EPSSetOperators(eps,_mat,NULL);
  //if(isSymmetric(tol))
	//EPSSetProblemType(eps,EPS_HEP);
  //else
	EPSSetProblemType(eps,EPS_NHEP);
  EPSSetWhichEigenpairs(eps,which);
  EPSSetDimensions(eps,nev,PETSC_DEFAULT,PETSC_DEFAULT);
  EPSSetTolerances(eps,tol,PETSC_DEFAULT);
  EPSSetType(eps,type);
  /*
     Set solver parameters at runtime
  */
  EPSSetFromOptions(eps);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  EPSSolve(eps);
  
  if(viewEigenvaluesInXWindows)
  {
	PetscViewer viewer = PETSC_VIEWER_DRAW_WORLD;
	PetscDraw draw;
	PetscViewerDrawGetDraw(viewer, 0, &draw);
	PetscDrawSetPause(draw, pause_lenght); // time duration of the display. if pause_lenght = -1 then wait for user to press a key
	std::string filename="MatrixEigenvaluePlot_"+matrixName+".ppm";
	PetscDrawSetSave( draw, filename.c_str());
	EPSValuesView(eps, viewer);

  }
	/*
	 Optional: Get some information from the solver and display it on the screen
	*/
	EPSGetIterationNumber(eps,&its);
	PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
	EPSGetType(eps,&type);
	PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
	EPSGetDimensions(eps,&nev,NULL,NULL);
	PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);
	EPSGetTolerances(eps,&tol,&maxit);
	PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
					Display solution and clean up
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	 Get number of converged approximate eigenpairs
	*/
	EPSGetConverged(eps,&nconv);
	PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);

  *valPr=new double[nconv];
  *valPi=new double[nconv];
  *vecPr=new double * [nconv];
  *vecPi=new double * [nconv];
  double * myvecPr, * myvecPi;
    
  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    if(viewEigenvaluesInXWindows)
      PetscPrintf(PETSC_COMM_WORLD,
           "           k          ||Ax-kx||/||kx||\n"
           "   ----------------- ------------------\n");

    for (int i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
      if(viewEigenvaluesInXWindows)
      {
	      /*
	         Compute the relative error associated to each eigenpair
	      */
	      if(fabs(kr)>tol || fabs(ki)>tol)
	      {
		      EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
		
		      if (ki!=0.0)
		        PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)kr,(double)ki,(double)error);
		      else
		        PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)kr,(double)error);
		   }
		   else
	      {
		      if (ki!=0.0)
		        PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12s\n",(double)kr,(double)ki,"Null eigenvalue");
		      else
		        PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12s\n",(double)kr,"Null eigenvalue");
		   }
      }	   
	  /* Copy real part of the spectrum */
      *(*valPr + i)=kr;
      VecGetArray(xr,&myvecPr);
      *(*vecPr+  i)=new double [_numberOfRows];
      memcpy(*(*vecPr+  i),myvecPr,_numberOfRows*sizeof(double)) ;
	  /* Copy imaginary part of the spectrum */
      *(*valPi + i)=ki;
      VecGetArray(xi,&myvecPi);
      *(*vecPi+  i)=new double [_numberOfRows];
      memcpy(*(*vecPi+  i),myvecPi,_numberOfRows*sizeof(double)) ;
    }
    if(viewEigenvaluesInXWindows)
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    /*
     Free work space
    */
    EPSDestroy(&eps);
    VecDestroy(&xr);
    VecDestroy(&xi);
    SlepcFinalize();

    return nconv;
  }
  else
  {
	/*
     Free work space
    */
    EPSDestroy(&eps);
    VecDestroy(&xr);
    VecDestroy(&xi);
    SlepcFinalize();

    throw CdmathException("SparseMatrixPetsc::computeSpectrum : No eigenvector found");	
  }	
}

int 
SparseMatrixPetsc::computeSVD(int nsv, double ** valS, double ***vecS, SVDWhich which, double tol, SVDType type, bool viewSingularValuesInXWindows, double pause_lenght, std::string matrixName) const
{
  SVD            svd;         /* Singular value decomposition solver context */
  PetscReal      error;
  PetscScalar    sigma;
  Vec            u,v;
  PetscInt       i,maxit,its, nconv;

  SlepcInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);

  MatAssemblyBegin(_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  _mat,MAT_FINAL_ASSEMBLY);

  MatCreateVecs(_mat,&u,NULL);
  MatCreateVecs(_mat,NULL,&v);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create singular value decomposition context
  */
  SVDCreate(PETSC_COMM_WORLD,&svd);

  /*
     Set operators. In this case, it is a standard singular value problem
  */
#if (SLEPC_VERSION_MAJOR==3) && (SLEPC_VERSION_MINOR > 14)
	SVDSetOperators(svd,_mat,NULL);
#else
	SVDSetOperator(svd,_mat);
#endif
  
  SVDSetWhichSingularTriplets(svd,which);
  SVDSetDimensions(svd,nsv,PETSC_DEFAULT,PETSC_DEFAULT);
  SVDSetType(svd, type);
  SVDSetTolerances(svd,tol,PETSC_DEFAULT);

  /*
     Set solver parameters at runtime
  */
  SVDSetFromOptions(svd);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the SVD problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  SVDSolve(svd);

  if(viewSingularValuesInXWindows)
  {
    PetscViewer viewer = PETSC_VIEWER_DRAW_WORLD;
    PetscDraw draw;
    PetscViewerDrawGetDraw(viewer, 0, &draw);
    PetscDrawSetPause(draw, pause_lenght); // time duration of the display. If pause_lenght = -1 then wait for user to press a key
    std::string filename="MatrixSingularValuePlot_"+matrixName+".ppm";
    PetscDrawSetSave( draw, filename.c_str());
    SVDValuesView(svd, viewer);
  }

  /*
     Optional: Get some information from the solver and display it
  */
  SVDGetIterationNumber(svd,&its);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
  SVDGetType(svd,&type);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
  SVDGetDimensions(svd,&nsv,NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested singular values: %D\n",nsv);
  SVDGetTolerances(svd,&tol,&maxit);
  PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate singular values
  */
  SVDGetConverged(svd,&nconv);
  PetscPrintf(PETSC_COMM_WORLD," Number of converged singular values: %D\n\n",nconv);

  *valS=new double[nconv];
  *vecS=new double * [2*nconv];
  double * myVecS;
  
  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    if(viewSingularValuesInXWindows)
      PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");

    for (int i=0;i<nconv;i++) {
      /*
        Get converged singular values: i-th eigenvalue is stored in valS
      */
      SVDGetSingularTriplet(svd,i,*valS+i, u, v);
      if(viewSingularValuesInXWindows)
      {
	      if(fabs(*(*valS+i))>tol)
	      {
		      /*
		         Compute the relative error associated to each singular value
		      */
		      SVDComputeError( svd, i, SVD_ERROR_RELATIVE, &error );
		
		      PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)*(*valS+i),(double)error);
			}
	       else
	          PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12s\n",(double)*(*valS+i),"Null singular value");
	  }
      VecGetArray(u,&myVecS);
      *(*vecS + i)=new double [_numberOfRows];
      memcpy(*(*vecS + i),myVecS,_numberOfRows*sizeof(double)) ;
      VecGetArray(v,&myVecS);
      *(*vecS + i + nconv)=new double [_numberOfColumns];
      memcpy(*(*vecS + i + nconv),myVecS,_numberOfColumns*sizeof(double)) ;
    }
    if(viewSingularValuesInXWindows)
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    /*
     Free work space
    */
    SVDDestroy(&svd);
    VecDestroy(&u);
    VecDestroy(&v);
    SlepcFinalize();

    return nconv;
  }
  else
  {
	/*
     Free work space
    */
    SVDDestroy(&svd);
    VecDestroy(&u);
    VecDestroy(&v);
    SlepcFinalize();

    throw CdmathException("SparseMatrixPetsc::computeSVD : singular value decomposition failed");	
  }	
}

std::vector< double > 
SparseMatrixPetsc::getEigenvalues(int nev, EPSWhich which, double tol, EPSType type, bool viewEigenvaluesInXWindows, double pause_lenght, std::string matrixName) const
{
	int nconv;
	double * valPr, *  valPi;
	double **vecPr, ** vecPi;

	nconv=computeSpectrum(nev, &valPr, &valPi, &vecPr, &vecPi, which, tol,type, viewEigenvaluesInXWindows, pause_lenght, matrixName);
	
    std::vector< double > result(nconv);
	
    for (int i=0;i<nconv;i++) 
        result[i]=valPr[i];

	delete[] valPr;
	delete[] valPi;
    for (int i=0;i<nconv;i++) 
    {
		delete[] vecPr[i];
		delete[] vecPi[i];
	}
	delete[] vecPr;
	delete[] vecPi;	
	
    return result;
}

std::vector< std::vector< double > >
SparseMatrixPetsc::plotEigenvalues(std::string matrixName, int nev, double pause_lenght, double tol, EPSWhich which, EPSType type) const
{
	double * valPr, *  valPi;
	double **vecPr, ** vecPi;
	std::vector< std::vector< double > > result(2);
	
	int nconv;

	if(nev <=0)
		nconv=computeSpectrum(_numberOfRows, &valPr, &valPi, &vecPr, &vecPi, which, tol,type, true, pause_lenght, matrixName);	
	else
		nconv=computeSpectrum(          nev, &valPr, &valPi, &vecPr, &vecPi, which, tol,type, true, pause_lenght, matrixName);
	
    std::vector< double > result_r(nconv);//real parts of the eigenvalues
    std::vector< double > result_i(nconv);//imaginary parts of the eigenvalues

    for (int i=0;i<nconv;i++) 
    {
        result_r[i]=valPr[i];
        result_i[i]=valPi[i];
	}
	result[0]=result_r;
	result[1]=result_i;
	
	delete[] valPr;
	delete[] valPi;
    for (int i=0;i<nconv;i++) 
    {
		delete[] vecPr[i];
		delete[] vecPi[i];
	}
	delete[] vecPr;
	delete[] vecPi;	
	
	return result;
}

std::vector< Vector > 
SparseMatrixPetsc::getEigenvectors(int nev, EPSWhich which, double tol, EPSType type) const
{
	int nconv;
	double * valPr, *  valPi;
	double **vecPr, ** vecPi;

	nconv=computeSpectrum(nev, &valPr, &valPi, &vecPr, &vecPi, which, tol,type);
	
    std::vector< Vector > result(nconv);

    for (int i=0;i<nconv;i++) 
    {
		DoubleTab values (_numberOfRows,vecPr[i]);
        Vector myvecPr(_numberOfRows);
        myvecPr.setValues(values);
        result[i]=myvecPr;
	}

	delete[] valPr;
	delete[] valPi;
    for (int i=0;i<nconv;i++) 
    {
		delete[] vecPr[i];
		delete[] vecPi[i];
	}
	delete[] vecPr;
	delete[] vecPi;	
	
    return result;
}

MEDCoupling::DataArrayDouble *
SparseMatrixPetsc::getEigenvectorsDataArrayDouble(int nev, EPSWhich which, double tol, EPSType type) const
{
	int nconv;
	double * valPr, *  valPi;
	double **vecPr, ** vecPi;

	nconv=computeSpectrum(nev, &valPr, &valPi, &vecPr, &vecPi, which, tol);
	
#ifdef MEDCoupling_VERSION_VERSION_GREATER_9_4
	std::vector< long unsigned int > compoId(1);
#else
	std::vector< int > compoId(1);
#endif
	MEDCoupling::DataArrayDouble *arrays=MEDCoupling::DataArrayDouble::New();
	MEDCoupling::DataArrayDouble *array =MEDCoupling::DataArrayDouble::New();
	arrays->alloc(_numberOfRows,nconv);	
	
    for (int i=0;i<nconv;i++) 
    {
		array->useArray(vecPr[i],false, MEDCoupling::DeallocType::CPP_DEALLOC, _numberOfRows,1);//array vecPr[i] will not be deallocated
		compoId[0]=i;
		arrays->setSelectedComponents(array,compoId);
		arrays->setInfoOnComponent(i,std::to_string(valPr[i]));
	}
	delete[] valPr;
	delete[] valPi;
    for (int i=0;i<nconv;i++) 
    {
		delete[] vecPr[i];
		delete[] vecPi[i];
	}
	delete[] vecPr;
	delete[] vecPi;	
	
    array->decrRef();
    return arrays;
}

double 
SparseMatrixPetsc::getConditionNumber(bool isSingular, double tol) const
{
	if(isSymmetric(tol))//Eigenvalues computation is more robust that singular values computation
	{
		int nev, nconv;
		double lambda_max, lambda_min;
		std::vector< double > my_ev;
		
		/*** Lowest eigenvalue, first check if matrix is singular ***/
	    if(isSingular)
	        nev=2;
	    else
			nev=1 ;
	
		my_ev=getEigenvalues( nev, EPS_SMALLEST_MAGNITUDE, tol);
		nconv=my_ev.size();
		if(nconv<nev)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the smallest eigenvalue");
	    lambda_min=my_ev[nev-1];
	    
	    /*** Largest eigenvalue ***/
	    nev=1;
		my_ev=getEigenvalues( nev, EPS_LARGEST_MAGNITUDE, tol);
		nconv=my_ev.size();
		if(nconv<nev)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the largest eigenvalue");
	    lambda_max=my_ev[nev-1];
	    
	    return fabs(lambda_max)/fabs(lambda_min);		
	}
	else//Singular values computation is more robust than eigenvalues computation
	{
		int nsv, nconv;
		double * valS;
		double ** vecS;
		double sigma_max, sigma_min;
		
		/*** Lowest singular value, first check if matrix is singular ***/
	    if(isSingular)
	        nsv=2;
	    else
			nsv=1 ;
	
		nconv=computeSVD(nsv, &valS, &vecS, SVD_SMALLEST, tol, SVDCROSS);
		if(nconv<nsv)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the smallest singular value");
	    sigma_min=valS[nsv-1];

	    delete[] valS;
	    for (int i=0;i<2*nconv;i++) 
			delete[] vecS[i];
	    delete[] vecS;
	    
	    /*** Largest singular value ***/
	    nsv=1;
		nconv=computeSVD(nsv, &valS, &vecS, SVD_LARGEST, tol, SVDCROSS);
		if(nconv<nsv)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the largest singular value");
	    sigma_max=valS[nsv-1];

	    delete[] valS;
	    for (int i=0;i<2*nconv;i++) 
			delete[] vecS[i];
	    delete[] vecS;
	    
	    return sigma_max/sigma_min;
	}
}

std::vector< double > 
SparseMatrixPetsc::getSingularValues(int nsv, SVDWhich which, double tol, SVDType type, bool viewSingularValuesInXWindows, double pause_lenght, std::string matrixName) const
{
	int nconv;
	double * valS;
	double **vecS;
	
	nconv=computeSVD(nsv, &valS, &vecS, which, tol, type,  viewSingularValuesInXWindows, pause_lenght, matrixName);
	
    std::vector< double > result(nconv);
	
    for (int i=0;i<nconv;i++) 
        result[i]=valS[i];

	delete[] valS;
	    for (int i=0;i<2*nconv;i++) 
			delete[] vecS[i];
	delete[] vecS;
	
    return result;
}

std::vector< Vector > 
SparseMatrixPetsc::getSingularVectors(int nsv, SVDWhich which, double tol, SVDType type) const
{
	int nconv;
	double * valS;
	double **vecS;
	
	nconv=computeSVD(nsv, &valS, &vecS, which, tol, type);
	
    std::vector< Vector > result(2*nconv);
	
    for (int i=0;i<nconv;i++) 
    {
		DoubleTab values_left (_numberOfRows,vecS[i]);
        Vector myVecS_left(_numberOfRows);
        myVecS_left.setValues(values_left);
        result[i]=myVecS_left;
		DoubleTab values_right (_numberOfColumns,vecS[i+nconv]);
        Vector myVecS_right(_numberOfColumns);
        myVecS_right.setValues(values_right);
        result[i+nconv]=myVecS_right;
	}

	delete[] valS;
    for (int i=0;i<2*nconv;i++) 
		delete[] vecS[i];
	delete[] vecS;	

    return result;
}



void 
SparseMatrixPetsc::leftDiagonalScale(Vector v)
{
	Vec vec=vectorToVec(v);

	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatDiagonalScale(_mat, vec, NULL);
}
void 
SparseMatrixPetsc::rightDiagonalScale(Vector v)
{
	Vec vec=vectorToVec(v);

	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatDiagonalScale(_mat, NULL, vec);
}
