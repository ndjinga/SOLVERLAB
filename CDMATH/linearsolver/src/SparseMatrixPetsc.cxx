/*
 * SparseMatrixPetsc.cxx
 *
 *  Created on: 04/11/2017
 *      Author: mndjinga
 */

#include "SparseMatrixPetsc.hxx"
#include "CdmathException.hxx"

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
SparseMatrixPetsc::viewMatrix() const 
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatView(_mat,PETSC_VIEWER_STDOUT_SELF);
}

void 
SparseMatrixPetsc::saveMatrix(string filename, bool binaryMode) const
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	PetscViewer fileViewer;
	PetscViewerCreate(PETSC_COMM_WORLD,&fileViewer);
    PetscViewerFileSetMode(fileViewer,FILE_MODE_WRITE);
    PetscViewerFileSetName(fileViewer,filename.c_str());

	if( binaryMode)
	{
		PetscViewerSetType(fileViewer, PETSCVIEWERBINARY);		
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &fileViewer);
	}
	else
	{
		PetscViewerSetType(fileViewer, PETSCVIEWERASCII);		
		PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &fileViewer);
	}
     
	MatView(_mat,fileViewer);
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
SparseMatrixPetsc::computeSpectrum(int nev, double ** valP, double ***vecP, EPSWhich which, double tol, EPSType type) const
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
  /*
     Optional: Get some information from the solver and display it
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

  *valP=new double[nconv];
  *vecP=new double * [nconv];
  double * myVecp;
    
  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");

    for (int i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
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
	   
      *(*valP + i)=kr;
      VecGetArray(xr,&myVecp);
      *(*vecP+  i)=new double [_numberOfRows];
      memcpy(*(*vecP+  i),myVecp,_numberOfRows*sizeof(double)) ;
    }
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
SparseMatrixPetsc::computeSVD(int nsv, double ** valS, double ***vecS, SVDWhich which, double tol, SVDType type) const
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
    PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");

    for (int i=0;i<nconv;i++) {
      /*
        Get converged singular values: i-th eigenvalue is stored in valS
      */
      SVDGetSingularTriplet(svd,i,*valS+i, u, v);
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
	
      VecGetArray(u,&myVecS);
      *(*vecS + i)=new double [_numberOfRows];
      memcpy(*(*vecS + i),myVecS,_numberOfRows*sizeof(double)) ;
      VecGetArray(v,&myVecS);
      *(*vecS + i + nconv)=new double [_numberOfColumns];
      memcpy(*(*vecS + i + nconv),myVecS,_numberOfColumns*sizeof(double)) ;
    }
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
SparseMatrixPetsc::getEigenvalues(int nev, EPSWhich which, double tol, EPSType type) const
{
	int nconv;
	double * valP;
	double **vecP;

	nconv=computeSpectrum(nev, &valP, &vecP, which, tol,type);
	
    std::vector< double > result(nconv);
	
    for (int i=0;i<nconv;i++) 
        result[i]=valP[i];

	delete[] valP;
    for (int i=0;i<nconv;i++) 
		delete[] vecP[i];
	delete[] vecP;	
	
    return result;
}

std::vector< Vector > 
SparseMatrixPetsc::getEigenvectors(int nev, EPSWhich which, double tol, EPSType type) const
{
	int nconv;
	double * valP;
	double **vecP;

	nconv=computeSpectrum(nev, &valP, &vecP, which, tol,type);
	
    std::vector< Vector > result(nconv);

    for (int i=0;i<nconv;i++) 
    {
		DoubleTab values (_numberOfRows,vecP[i]);
        Vector myVecP(_numberOfRows);
        myVecP.setValues(values);
        result[i]=myVecP;
	}

	delete[] valP;
    for (int i=0;i<nconv;i++) 
		delete[] vecP[i];
	delete[] vecP;	
	
    return result;
}

MEDCoupling::DataArrayDouble *
SparseMatrixPetsc::getEigenvectorsDataArrayDouble(int nev, EPSWhich which, double tol, EPSType type) const
{
	int nconv;
	double * valP;
	double **vecP;

	nconv=computeSpectrum(nev, &valP, &vecP, which, tol);
	
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
		array->useArray(vecP[i],true, MEDCoupling::DeallocType::CPP_DEALLOC, _numberOfRows,1);
		compoId[0]=i;
		arrays->setSelectedComponents(array,compoId);
		arrays->setInfoOnComponent(i,std::to_string(valP[i]));
	}
	delete[] valP;
	delete[] vecP;	
	
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
	
		nconv=computeSVD(nsv, &valS, &vecS, SVD_SMALLEST, tol,SVDCYCLIC);
		if(nconv<nsv)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the smallest singular value");
	    sigma_min=valS[nsv-1];
	    delete[] valS, vecS;
	    
	    /*** Largest singular value ***/
	    nsv=1;
		nconv=computeSVD(nsv, &valS, &vecS, SVD_LARGEST, tol,SVDCYCLIC);
		if(nconv<nsv)
			throw CdmathException("SparseMatrixPetsc::getConditionNumber could not find the largest singular value");
	    sigma_max=valS[nsv-1];
	    delete[] valS, vecS;
	    
	    return sigma_max/sigma_min;
	}
}

std::vector< double > 
SparseMatrixPetsc::getSingularValues(int nsv, SVDWhich which, double tol, SVDType type) const
{
	int nconv;
	double * valS;
	double **vecS;
	
	nconv=computeSVD(nsv, &valS, &vecS, which, tol, type);
	
    std::vector< double > result(nconv);
	
    for (int i=0;i<nconv;i++) 
        result[i]=valS[i];

	delete[] valS, vecS;
	
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
    for (int i=0;i<nconv;i++) 
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
