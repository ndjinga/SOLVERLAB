/**
 * LinearSolver.cxx
 *
 *  Created on: 15 avr. 2013
 *      Author: mekkas
 */

#include <cmath>
#include <string>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include "CdmathException.hxx"
#include "LinearSolver.hxx"
#include "SparseMatrixPetsc.hxx"

using namespace std;


LinearSolver::LinearSolver ( void )
{
	_tol=1.E-15;
	_numberMaxOfIter=0;
	_residu=1.E30;
	_convergence=false;
	_numberOfIter=0;
	_isSingular=false;
	_nameOfPc="";
	_nameOfMethod="";
	_mat=NULL;
	_smb=NULL;
	_solution=NULL;
	_ksp=NULL;
	_prec=NULL;
	_isSparseMatrix=false;
	_computeConditionNumber=false;
    _conditionNumber=1e100;
}

void
LinearSolver::kspDuplicate(const KSP source, const Mat mat, KSP* destination) const
{
	KSPCreate(PETSC_COMM_WORLD,&(*destination));
#ifdef PETSC_VERSION_GREATER_3_5
	KSPSetOperators(*destination,mat,mat);
#else
	KSPSetOperators(*destination,mat,mat,SAME_NONZERO_PATTERN);
#endif
	KSPType type;
	KSPGetType(source,&type);
	KSPSetType(*destination,type);
	/*
    PetscReal tol1,tol2,tol3;
    PetscInt maxIter;
    KSPGetTolerances(source,&tol1,&tol2,&tol3,&maxIter);
    KSPSetTolerances(*destination,tol1,tol2,tol3,maxIter);
	// */
}

void
LinearSolver::precDuplicate(const PC source, const KSP ksp, PC* destination)  const
{
	KSPGetPC(ksp,destination);
	PCType type;
	PCGetType(source,&type);
	PCSetType(*destination,type);
}

LinearSolver::~LinearSolver ( void )
{
	//if(&_mat != NULL)
	//	MatDestroy(&_mat);
	if(_smb != NULL)
		VecDestroy(&_smb);
	if(_solution != NULL)
		VecDestroy(&_solution);
	KSPDestroy(&_ksp);
	//PetscFinalize();
}

void
LinearSolver::setTolerance(double tol)
{
	_tol=tol;
	KSPSetTolerances(_ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,getNumberMaxOfIter());
}

void
LinearSolver::setNumberMaxOfIter(int numberMaxOfIter)
{
	_numberMaxOfIter=numberMaxOfIter;
	KSPSetTolerances(_ksp,getTolerance(),PETSC_DEFAULT,PETSC_DEFAULT,numberMaxOfIter);
}

void
LinearSolver::setComputeConditionNumber(bool display)
{
	_computeConditionNumber=display;
}

double 
LinearSolver::getConditionNumber() const
{
    if(_computeConditionNumber and _convergence and !_isSingular)
        return _conditionNumber;
    else if(!_computeConditionNumber)
        throw CdmathException("LinearSolver::getConditionNumber(): Condition number can not be evaluated without prior call to setComputeConditionNumber()");
    else if(_isSingular)//matrix is singular
        throw CdmathException("LinearSolver::getConditionNumber(): Condition number can not be evaluated by PETSC for singular matrices (division by zero)");
    else //_convergence = false
        throw CdmathException("LinearSolver::getConditionNumber(): Condition number can not be evaluated without prior successful resolution of a linear system");
}

LinearSolver::LinearSolver( const GenericMatrix& matrix,
		const Vector& secondMember,
		int numberMaxOfIter,
		double tol,
		string nameOfMethod,
		string nameOfPc )
{
	_tol = tol;
	_numberMaxOfIter = numberMaxOfIter;
	_residu = 1.E30;
	_convergence = false;
	_numberOfIter = 0;
	_isSingular = false;
	_isSparseMatrix = matrix.isSparseMatrix();
	_computeConditionNumber=false;
	_nameOfPc = nameOfPc;
	_nameOfMethod = nameOfMethod;
	_secondMember = secondMember;
	_mat = NULL;
	_smb = NULL;
	_solution = NULL;
	_prec = NULL;
	_ksp = NULL;

	//setTolerance(tol);
	//setNumberMaxOfIter(numberMaxOfIter);
	setPreconditioner(nameOfPc);
	setMethod(nameOfMethod);
	setLinearSolver(matrix, secondMember);
}

LinearSolver::LinearSolver( const GenericMatrix& matrix,
		const std::vector<double>& secondMember,
		int numberMaxOfIter,
		double tol,
		string nameOfMethod,
		string nameOfPc )
{
	_tol = tol;
	_numberMaxOfIter = numberMaxOfIter;
	_residu = 1.E30;
	_convergence = false;
	_numberOfIter = 0;
	_isSingular = false;
	_isSparseMatrix = matrix.isSparseMatrix();
	_computeConditionNumber=false;
	_nameOfPc = nameOfPc;
	_nameOfMethod = nameOfMethod;
	_mat = NULL;
	_smb = NULL;
	_solution = NULL;
	_prec = NULL;
	_ksp = NULL;

	//setTolerance(tol);
	//setNumberMaxOfIter(numberMaxOfIter);
	setPreconditioner(nameOfPc);
	setMethod(nameOfMethod);
    //converting vector<double> to Vector
    int size=secondMember.size();
    Vector secondMemberVector(size);
 	for ( int i=0; i<size; i++)
		secondMemberVector(i)=secondMember[i];
	_secondMember = secondMemberVector;

	setLinearSolver(matrix, secondMemberVector);
}

void
LinearSolver::setPreconditioner(string pc)
{
	if((pc.compare("ICC") != 0) && (pc.compare("ILU") != 0) && (pc.compare("LU") != 0) && (pc.compare("CHOLESKY") != 0) && (pc.compare("") != 0))
	{
		string msg="LinearSolver::LinearSolver: preconditioner "+pc+" does not exist.\n";
		throw CdmathException(msg);
	}
	if ((pc.compare("ICC")==0 || pc.compare("ILU")==0) && _isSparseMatrix==false)
	{
		string msg="LinearSolver::LinearSolver: preconditioner "+pc+" is not compatible with dense matrix.\n";
		throw CdmathException(msg);
	}

	if ((pc.compare("ICC")==0 || pc.compare("ILU")==0) && (_nameOfMethod.compare("LU")==0 || _nameOfMethod.compare("CHOLESKY")==0 ))
	{
		string msg="LinearSolver::LinearSolver: preconditioner "+pc+" is not compatible with "+_nameOfMethod+".\n";
		throw CdmathException(msg);
	}
	_nameOfPc=pc;
}


void
LinearSolver::setMethod(string nameOfMethod)
{
	_nameOfMethod = nameOfMethod;

	if ((_nameOfPc.compare("ICC")==0 || _nameOfPc.compare("ILU")==0) && (_nameOfMethod.compare("LU")==0 || _nameOfMethod.compare("CHOLESKY")==0) )
	{
		string msg="LinearSolver::LinearSolver: preconditioner "+_nameOfPc+" is not compatible with "+_nameOfMethod+".\n";
		throw CdmathException(msg);
	}
}


void
LinearSolver::setLinearSolver(const GenericMatrix& matrix, const Vector& secondMember)
{
	if ((_nameOfPc.compare("ICC")==0 || _nameOfPc.compare("ILU")==0) && _isSparseMatrix==false)
	{
		string msg="LinearSolver::LinearSolver: preconditioner "+_nameOfPc+" is not compatible with dense matrix.\n";
		throw CdmathException(msg);
	}

	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);//All constructors lead here so we initialize petsc here
	setMatrix(matrix);
	setSndMember(secondMember);
	VecDuplicate(_smb,&_solution);
}


bool
LinearSolver::isSparseMatrix( void ) const
{
	return (_isSparseMatrix);
}

bool
LinearSolver::isMatrixSingular( void ) const
{
	return _isSingular;
}

int
LinearSolver::getNumberOfIter( void ) const
{
	return (_numberOfIter);
}

bool
LinearSolver::getStatus( void ) const
{
	return (_convergence);
}

double
LinearSolver::getResidu( void ) const
{
	return (_residu);
}

double
LinearSolver::getTolerance(void) const
{
	return (_tol);
}

int
LinearSolver::getNumberMaxOfIter(void) const
{
	return (_numberMaxOfIter);
}

void
LinearSolver::setMatrix(const GenericMatrix& matrix)
{
	/* matrix to mat */
	int numberOfRows=matrix.getNumberOfRows();
	int numberOfColumns=matrix.getNumberOfColumns();

	if (matrix.isSparseMatrix())//Matrix is already a petsc structure
	{
			const SparseMatrixPetsc& Smat = dynamic_cast<const SparseMatrixPetsc&>(matrix);
			_mat=Smat.getPetscMatrix();
	}
	else//Matrix is in A. Mekkas dense structure
	{
		MatCreate(PETSC_COMM_WORLD, &_mat);
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, numberOfRows, numberOfColumns);
		MatSetType(_mat,MATSEQDENSE);

		PetscScalar *a;
		PetscMalloc(numberOfRows*numberOfColumns*sizeof(PetscScalar),&a);
		MatSeqDenseSetPreallocation(_mat,a);

		for (int i=0;i<numberOfRows;i++)
			for (int j=0;j<numberOfColumns;j++)
				a[i+j*numberOfRows]=matrix(i,j);
	}
	//Assemblage final
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	KSPCreate(PETSC_COMM_WORLD, &_ksp);
#ifdef PETSC_VERSION_GREATER_3_5
KSPSetOperators(_ksp,_mat,_mat);
#else
	KSPSetOperators(_ksp,_mat,_mat,SAME_NONZERO_PATTERN);
#endif

	KSPGetPC(_ksp,&_prec);
}

void
LinearSolver::setSndMember(const Vector& secondMember)
{
	_secondMember=secondMember;
	if(&_smb!=NULL)
		VecDestroy(&_smb);
	_smb=vectorToVec(secondMember);

}

void
LinearSolver::setMatrixIsSingular(bool sing)
{
	_isSingular=sing;
}

Vector
LinearSolver::getSndMember(void) const
{
	return (_secondMember);
}

string
LinearSolver::getNameOfMethod(void) const
{
	return (_nameOfMethod);
}

string
LinearSolver::getNameOfPc(void) const
{
	return (_nameOfPc);
}

LinearSolver::LinearSolver ( const LinearSolver& LS )
{
	_tol=LS.getTolerance();
	_nameOfMethod=LS.getNameOfMethod();
	_numberMaxOfIter=LS.getNumberMaxOfIter();
	_secondMember=LS.getSndMember();
	_residu=LS.getResidu();
	_convergence=LS.getStatus();
	_numberOfIter=LS.getNumberOfIter();
	_isSingular=LS.isMatrixSingular();
	_nameOfPc=LS.getNameOfPc();
	_mat=NULL;
	MatDuplicate(LS.getPetscMatrix(),MAT_COPY_VALUES,&_mat);
	_smb=NULL;
	VecDuplicate(LS.getPetscVector(),&_smb);    				;
	_solution=NULL;
	VecDuplicate(LS.getKSPSolution(),&_solution);
	_ksp=NULL;
	kspDuplicate(LS.getPetscKsp(),_mat,&_ksp);
	_prec=NULL;
	precDuplicate(LS.getPetscPc(),_ksp,&_prec);
	//    _ksp=LS.getPetscKsp();
	//    _prec=LS.getPetscPc();
	_isSparseMatrix=LS.isSparseMatrix();
}

KSP
LinearSolver::getPetscKsp() const
{
	return (_ksp);
}

Mat
LinearSolver::getPetscMatrix() const
{
	return (_mat);
}

Vec
LinearSolver::getPetscVector() const
{
	return (_smb);
}

void
LinearSolver::viewPetscMatrix() const 
{
	MatView(_mat,PETSC_VIEWER_STDOUT_SELF);
}
void
LinearSolver::viewPetscRHS() const 
{
	VecView(_smb,PETSC_VIEWER_STDOUT_SELF);
}
double
LinearSolver::getPetscMatValue(int i, int j) const 
{
	double res;
	int idxm=i,idxn=j;
	MatGetValues(_mat,1,&idxm,1, &idxn,&res);
	return res;
}
double
LinearSolver::getPetscRHSValue(int i) const 
{
	double res;
	int idxm=i;
	VecGetValues(_smb,1,&idxm,&res);
	return res;
}
PC
LinearSolver::getPetscPc() const
{
	return (_prec);
}

Vector
LinearSolver::solve( Vector X0 )
{
	if (_nameOfMethod.compare("GMRES")==0)
		KSPSetType(_ksp,KSPGMRES);
	else if (_nameOfMethod.compare("FGMRES")==0)
		KSPSetType(_ksp,KSPFGMRES);
	else if (_nameOfMethod.compare("CG")==0)
		KSPSetType(_ksp,KSPCG);
	else if (_nameOfMethod.compare("CGNE")==0)
		KSPSetType(_ksp,KSPCGNE);
	else if (_nameOfMethod.compare("BCGS")==0)
		KSPSetType(_ksp,KSPBCGS);
	else if (_nameOfMethod.compare("CR")==0)
		KSPSetType(_ksp,KSPCR);
	else if (_nameOfMethod.compare("CGS")==0)
		KSPSetType(_ksp,KSPCGS);
	else if (_nameOfMethod.compare("BICG")==0)
		KSPSetType(_ksp,KSPBICG);
	else if (_nameOfMethod.compare("GCR")==0)
		KSPSetType(_ksp,KSPGCR);
	else if (_nameOfMethod.compare("LSQR")==0)
		KSPSetType(_ksp,KSPLSQR);
	else if (_nameOfMethod.compare("CHOLESKY")==0)
	{
		KSPSetType(_ksp,KSPGMRES);
		PCSetType(_prec,PCCHOLESKY);
	}
	else if (_nameOfMethod.compare("LU")==0)
	{
		KSPSetType(_ksp,KSPGMRES);
		PCSetType(_prec,PCLU);
	}
	else
	{
		string msg="Vector LinearSolver::solve( void ) : The method "+_nameOfMethod+" is not yet implemented.\n";
		msg+="The methods implemented are : GMRES, BICG, CG, CGNE, CHOLESKY, LU, BCGS, FGMRES, LSQR, CR, CGS and GCR.\n";
		throw CdmathException(msg);
	}

	if (_nameOfPc.compare("ILU")==0)
		PCSetType(_prec,PCILU);
	else if (_nameOfPc.compare("LU")==0)
		PCSetType(_prec,PCLU);
	else if (_nameOfPc.compare("ICC")==0)
		PCSetType(_prec,PCICC);
	else if (_nameOfPc.compare("CHOLESKY")==0)
		PCSetType(_prec,PCCHOLESKY);
#if SOLVERLAB_WITH_MPI 
	else if (_nameOfPc.compare("euclid")==0 || _nameOfPc.compare("pilut")==0 || _nameOfPc.compare("parasails")==0 || _nameOfPc.compare("ams")==0 || _nameOfPc.compare("ads")==0)
	{
		PCSetType(_prec,PCHYPRE);
		PCHYPRESetType(_prec,_nameOfPc.c_str());
	}
	else if (_nameOfPc.compare("BOOMERAMG")==0)
	{
		PCSetType(_prec,PCHYPRE);
		PetscOptionsInsertString(NULL,"-pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7  -pc_hypre_boomeramg_agg_nl 4 -pc_hypre_boomeramg_agg_num_paths 5 -pc_hypre_boomeramg_max_levels 25 -pc_hypre_boomeramg_coarsen_type HMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_P_max 2 -pc_hypre_boomeramg_truncfactor 0.3");
	}
#endif
	else if (_nameOfPc.compare("")==0)
		PCSetType(_prec,PCNONE);
	else
	{
		string msg="Vector LinearSolver::solve( void ) : The preconditioner "+_nameOfPc+" is not yet available.\n";
		msg+="The preconditioners available are : void, ICC, ILU, CHOLESKY, LU";
#if SOLVERLAB_WITH_MPI 
		msg+=", euclid, pilut, parasails, ams, ads, BOOMERAMG";
#endif
		msg+=".\n";
		throw CdmathException(msg);
	}

	KSPSetTolerances(_ksp,_tol*_tol*_tol*_tol,_tol,PETSC_DEFAULT,_numberMaxOfIter);

	PetscInt its;
	PetscReal atol;
	PetscInt maxits;

	if( X0.size()==0)
		KSPSetInitialGuessNonzero( _ksp,PETSC_FALSE);
	else
	{
		KSPSetInitialGuessNonzero( _ksp,PETSC_TRUE);
		_solution=vectorToVec(X0);
	}
	
	if (_isSingular)//Matrix should be symmetric semi-definite with constant vectors in its kernel
	{
		//Check that the matrix is symmetric
		PetscBool isSymetric;
		MatIsSymmetric(_mat,_tol,&isSymetric);
		if(!isSymetric)
			{
				cout<<"Singular matrix is not symmetric, tolerance= "<< _tol<<endl;
				throw CdmathException("Singular matrix should be symmetric with kernel composed of constant vectors");
			}
		else
		{
			MatSetOption(_mat,MAT_HERMITIAN, PETSC_TRUE);
			MatNullSpace nullsp;
			MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);//Declaration of a kernel containing exclusively constant vectors
			MatSetTransposeNullSpace(_mat, nullsp);//Transpose matrix has the same kernel since the matrix is symmetric
			MatSetNullSpace(_mat, nullsp);
			MatNullSpaceDestroy(&nullsp);
		}
	}

	if(_computeConditionNumber)
		KSPSetComputeSingularValues(_ksp,PETSC_TRUE);

	KSPSolve(_ksp,_smb,_solution);

	KSPGetResidualNorm(_ksp,&atol);
	_residu=(double)atol;

	KSPGetIterationNumber(_ksp,&its);
	_numberOfIter=(int)its;

	KSPConvergedReason reason;
	KSPGetConvergedReason(_ksp,&reason);

	if (reason>0)
		_convergence=true;
	else{
		_convergence=false;
		cout<<endl<<"!!!!!!!!!!!!!!!!!!  Linear system algorithm did not converge  !!!!!!!!!!!!!!" <<endl;
		if( reason == -3)
		    cout<<"Maximum number of iterations "<<_numberMaxOfIter<<" reached"<<endl;
		else if( reason == -11)
		    cout<<"!!!!!!! Construction of preconditioner failed !!!!!!"<<endl;
		else if( reason == -5)
		    cout<<"!!!!!!! Generic breakdown of the linear solver (Could be due to a singular matrix or preconditioner)!!!!!!"<<endl;
		else
		{
		    cout<<"PETSc divergence reason  "<< reason <<endl;
			cout<<"Final iteration= "<<_numberOfIter<<". Maximum allowed was " << _numberMaxOfIter<<endl;
		}
        cout<<"Solver used "<<  _nameOfMethod<<", with preconditioner "<<_nameOfPc<<endl;
		cout<<"Final residual "<< _residu<< ". Objective was "<< _tol<<endl;
		string msg="Linear system algorithm did not converge";
		throw CdmathException(msg);
	}

	if(_computeConditionNumber and !_isSingular)
	{
		PetscReal sv_max, sv_min;
		KSPComputeExtremeSingularValues(_ksp, &sv_max, &sv_min);
		cout<<" Maximal singular value = " << sv_max <<", Minimal singular value = " << sv_min <<", Condition number = " << sv_max/sv_min <<endl;
        _conditionNumber=sv_max/sv_min;
	}

	Vector X1=vecToVector(_solution);

	return X1;
}

Vec
LinearSolver::vectorToVec(const Vector& myVector) const
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

Vector
LinearSolver::vecToVector(const Vec& vec) const
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

Vector
LinearSolver::getResidual( Vector X ) const
{//This is adapted from PETSc source code of KSPInitialResidual
 Mat  Amat, Pmat;
 Vec  vb, vt1, vt2, vres;//vb=rhs, vres==resultat, vt1,vt2=temporary storage
 PCSide side;

 KSPGetPCSide(_ksp,&side);
 PCGetOperators(_prec,&Amat,&Pmat);
 KSPGetRhs( _ksp,&vb);

    if ( X.size()!=0 ) {
		//check size compatibility
		int numberOfRows;
		MatGetSize(_mat,&numberOfRows,NULL);
		if( X.size()!= numberOfRows )
			throw CdmathException("LinearSolver::getResidual input vector has incorrect size");

      Vec vsoln = vectorToVec(X);
      MatMult(Amat,vsoln,vt1);
      VecCopy(vb,vt2);
      VecAXPY(vt2,-1.0,vt1);
      if ( side == PC_RIGHT) {
        PCDiagonalScaleLeft(_prec,vt2,vres);
      } else if ( side == PC_LEFT) {
        PCApply(_prec,vt2,vres);
        PCDiagonalScaleLeft(_prec,vres,vres);
      } else if ( side == PC_SYMMETRIC) {
        PCApplySymmetricLeft(_prec,vt2,vres);
      } else throw CdmathException( "LinearSolver::getResidual : Invalid preconditioning side ");
    } else {
		//check initial guess zero
		PetscBool has_non_zero_guess;
		KSPGetInitialGuessNonzero( _ksp, &has_non_zero_guess);
		if(has_non_zero_guess)
			throw CdmathException( "LinearSolver::getResidual : initial guess was not zero ");

      VecCopy(vb,vt2);//A quoi sert cette ligne ???
      if ( side == PC_RIGHT) {
        PCDiagonalScaleLeft(_prec,vb,vres);
      } else if ( side == PC_LEFT) {
        PCApply(_prec,vb,vres);
        PCDiagonalScaleLeft(_prec,vres,vres);
      } else if ( side == PC_SYMMETRIC) {
        PCApplySymmetricLeft(_prec, vb, vres);
      } 
      else 
		throw CdmathException( "LinearSolver::getResidual : Invalid preconditioning side ");
    }	
    
    return vecToVector(vres);
}

Vector
LinearSolver::getInitialResidual( ) const
{
	PetscBool has_non_zero_guess;
	KSPGetInitialGuessNonzero( _ksp, &has_non_zero_guess);
	if(has_non_zero_guess)
		throw CdmathException( "LinearSolver::getInitialResidual : initial guess was not zero ");
	
 Mat  Amat, Pmat;
 Vec  vb, vres;//vb=rhs, vres==resultat
 PCSide side;

 KSPGetPCSide(_ksp,&side);
 KSPGetRhs( _ksp,&vb);
 PCGetOperators(_prec,&Amat,&Pmat);

  if ( side == PC_RIGHT)
        PCDiagonalScaleLeft(_prec,vb,vres);
  else if ( side == PC_LEFT){
        PCApply(_prec,vb,vres);//Pb
        PCDiagonalScaleLeft(_prec,vres,vres);}
  else if ( side == PC_SYMMETRIC)
        PCApplySymmetricLeft(_prec, vb, vres);
  else 
        throw CdmathException( "LinearSolver::getInitialResidual : Invalid preconditioning side ");
    
    return vecToVector(vres);
}

Vector
LinearSolver::getFinalResidual( ) const
{
	if( !_convergence )
        throw CdmathException( "LinearSolver::getFinalResidual : solve a linear system first ");

	Mat  Amat, Pmat;
	Vec  vb, vsoln, vt1, vres;//vsoln=solution, vb=rhs, vres==resultat, vt12=temporary storage
	PCSide side;
	
	KSPGetPCSide(_ksp,&side);
	KSPGetSolution( _ksp,&vsoln);
	KSPGetRhs( _ksp,&vb);
	PCGetOperators(_prec,&Amat,&Pmat);
	
	MatMult(Amat,vsoln,vt1);//AX
	VecAXPY(vb,-1.0,vt1);//b-AX
	if ( side == PC_RIGHT) {
		PCDiagonalScaleLeft(_prec,vb,vres);
	} 
	else if ( side == PC_LEFT) {
		PCApply(_prec,vb,vres);//P(b-AX)
		PCDiagonalScaleLeft(_prec,vres,vres);
	} 
	else if ( side == PC_SYMMETRIC) {
		PCApplySymmetricLeft(_prec,vb,vres);
	} else 
		throw CdmathException( "LinearSolver::getFinalResidual : Invalid preconditioning side ");
    
    return vecToVector(vres);
}

//----------------------------------------------------------------------
const LinearSolver&
LinearSolver::operator= ( const LinearSolver& linearSolver )
//----------------------------------------------------------------------
{
	_tol=linearSolver.getTolerance();
	_numberMaxOfIter=linearSolver.getNumberMaxOfIter();
	_nameOfMethod=linearSolver.getNameOfMethod();
	_residu=linearSolver.getResidu();
	_convergence=linearSolver.getStatus();
	_numberOfIter=linearSolver.getNumberOfIter();
	_isSingular=linearSolver.isMatrixSingular();
	_secondMember=linearSolver.getSndMember();
	_isSparseMatrix=linearSolver.isSparseMatrix();
	_nameOfPc="";
	setPreconditioner(linearSolver.getNameOfPc());
	_mat=NULL;
	MatDuplicate(linearSolver.getPetscMatrix(),MAT_COPY_VALUES,&_mat);
	_smb=NULL;
	VecDuplicate(linearSolver.getPetscVector(),&_smb);    				;
	_solution=NULL;
	VecDuplicate(linearSolver.getKSPSolution(),&_solution);    				;
	_ksp=NULL;
	kspDuplicate(linearSolver.getPetscKsp(),_mat,&_ksp);
	_prec=NULL;
	precDuplicate(linearSolver.getPetscPc(),_ksp,&_prec);
	return (*this);
}
