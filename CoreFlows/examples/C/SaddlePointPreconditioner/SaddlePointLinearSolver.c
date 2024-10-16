static char help[] = "Read a PETSc matrix from a file -f0 <input file>\n Parameters : \n -f0 : matrix fileName \n -nU :number of velocity lines \n -nP : number of pressure lines \n -mat_type : PETSc matrix type \n";

/*************************************************************************************************/
/* Sequential implementation of a new preconditioner for the linear system A_{input} X_{output} = b_{input} */
/*                                                                                               */
/* Input  : - Matrix A_{input}    (system matrix, loaded from a file)                            */
/*          - Vector b_{input}    (right hand side, made up for testing)                         */
/* Output : - Vector X_{output}   (unknown vector, to be determined                              */
/*                                                                                               */
/* Auxilliary variables : - A_hat (transformed matrix)                                           */
/*                        - X_hat (unknown of the transformed system)                            */
/*                        - b_hat (RHS of the transformed system)                                */
/*                        - Pmat  (preconditioning matrix)                                       */
/*                        - M top    left  submatrix of A_{input}                                */
/*                        - G top    right submatrix of A_{input}                                */
/*                        - D bottom left  submatrix of A_{input}                                */
/*                        - C bottom right submatrix of A_{input}                                */
/*                                                                                               */
/*                                 *M   G*                                                       */
/*                        A     = *       *                                                      */
/*                                 *D   C*                                                       */
/*                                                                                               */
/*                                 *C_hat  -D*                                                   */
/*                        A_hat = *           *                                                  */
/*                                 *G_hat   M*                                                   */
/*                                                                                               */
/*                                 *C_hat        0   *                                           */
/*                        Pmat  = *                   *                                          */
/*                                 *G_hat   2 diag(M)*                                           */
/*                                                                                               */
/*************************************************************************************************/

#include <petscis.h>
#include <petscksp.h>

int main( int argc, char **args ){
	PetscInitialize(&argc,&args, (char*)0,help);
	PetscMPIInt    size;        /* size of communicator */
	PetscMPIInt    rank;        /* processor rank */
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	PetscErrorCode ierr=0;

	PetscCheck( size == 1, PETSC_COMM_WORLD, ierr, "Incorrect number of procs nprocs = %d.\n !!! This is a sequential implementation !!! \n", size);

//##### Load the matrix A in the file given in the command line
	char file[1][PETSC_MAX_PATH_LEN], mat_type[256]; // File to load, matrix type
	PetscViewer viewer;
	Mat A_input;
	PetscBool flg;

	PetscOptionsGetString(NULL,NULL,"-f0",file[0],PETSC_MAX_PATH_LEN,&flg);
	PetscStrcpy(mat_type,MATAIJ);// Default value for PETSc Matrix type
	PetscOptionsGetString(NULL,NULL,"-mat_type",mat_type,sizeof(mat_type),NULL);

	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix type %s from file %s on %d processor(s)...\n", mat_type, file[0], size);	

	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);	
	PetscViewerSetType(viewer,PETSCVIEWERBINARY);
	PetscViewerFileSetMode(viewer,FILE_MODE_READ);
	PetscViewerFileSetName(viewer,file[0]);
	
	MatCreate(PETSC_COMM_WORLD, &A_input);
	MatSetType(A_input,mat_type);
	MatLoad(A_input,viewer);
	PetscViewerDestroy(&viewer);

	PetscPrintf(PETSC_COMM_WORLD,"... matrix Loaded \n");	
	

//####	Decompose the matrix A_input into 4 blocks M, G, D, C
	Mat M, G, D, C;
	PetscInt nrows, ncolumns;//Total number of rows and columns of A_input
	PetscInt n_u, n_p, n;//Total number of velocity and pressure lines. n = n_u+ n_p
	IS is_U,is_P;

	PetscOptionsGetInt(NULL,NULL,"-nU",&n_u,NULL);
	PetscOptionsGetInt(NULL,NULL,"-nP",&n_p,NULL);
	n=n_u+n_p;
	MatGetSize( A_input, &nrows, &ncolumns);

	PetscCheck( nrows == ncolumns, PETSC_COMM_WORLD, ierr, "Matrix is not square !!!\n");
	PetscCheck( n == ncolumns, PETSC_COMM_WORLD, ierr, "Inconsistent data : the matrix has %d lines but only %d velocity lines and %d pressure lines declared\n", ncolumns, n_u,n_p);
	PetscPrintf(PETSC_COMM_WORLD,"The matrix has %d lines : %d velocity lines and %d pressure lines\n", n, n_u,n_p);
	PetscPrintf(PETSC_COMM_SELF,"Matrix size : %d x %d, n_u = %d, n_p = %d \n", nrows, ncolumns, n_u, n_p);
	
	PetscPrintf(PETSC_COMM_WORLD,"Extraction of the 4 blocks \n");
	ISCreateStride(PETSC_COMM_WORLD, n_u, 0, 1, &is_U);
	ISCreateStride(PETSC_COMM_WORLD, n_p, n_u, 1, &is_P);

	MatCreateSubMatrix(A_input,is_U, is_U,MAT_INITIAL_MATRIX,&M);
	MatCreateSubMatrix(A_input,is_U, is_P,MAT_INITIAL_MATRIX,&G);
	MatCreateSubMatrix(A_input,is_P, is_U,MAT_INITIAL_MATRIX,&D);
	MatCreateSubMatrix(A_input,is_P, is_P,MAT_INITIAL_MATRIX,&C);
	PetscPrintf(PETSC_COMM_WORLD,"... end of extraction\n");

//##### Definition of the right hand side to test the preconditioner
	Vec b_input, b_input_p, b_input_u, b_hat, X_hat, X_anal;
	Vec X_array[2];
	PetscScalar y[n_p];
	PetscInt i_p[n_p];

	PetscPrintf(PETSC_COMM_WORLD,"Creation of the RHS, exact and numerical solution vectors...\n");
	VecCreate(PETSC_COMM_WORLD,&b_input);
	VecSetSizes(b_input,PETSC_DECIDE,n_u+n_p);
	VecSetFromOptions(b_input);

	VecDuplicate(b_input,&X_anal);//X_anal will store the exact solution
	VecDuplicate(b_input,&X_hat);// X_hat will store the numerical solution of the transformed system
	VecDuplicate(b_input,&b_hat);// b_hat will store the right hand side of the transformed system
	
	VecSet(X_anal,0.0);
	for (int i = n_u;i<n;i++){
		y[i-n_u]=1.0/i;
		i_p[i-n_u]=i;
	}
	VecSetValues(X_anal,n_p,i_p,y,INSERT_VALUES);
	VecAssemblyBegin(X_anal);
	VecAssemblyEnd(X_anal);
	VecNormalize( X_anal, NULL);
	MatMult( A_input, X_anal, b_input);
	PetscPrintf(PETSC_COMM_WORLD,"... vectors created \n");	
	//MatDestroy(&A_input);//Early destruction since A_input is a sequential matrix stored on processed 0

	//Swap the pressure and velocity components + change the sign of the pressure components of b_input (this is due to the change in ordering of the variable in pierre-loic original script)
	VecGetSubVector( b_input, is_P, &b_input_p);
	VecGetSubVector( b_input, is_U, &b_input_u);
	//VecScale(b_input_p, -1);
	X_array[0] = b_input_p;
	X_array[1] = b_input_u;

	//VecCreateNest( PETSC_COMM_WORLD, 2, NULL, X_array, &b_hat);//This may generate an error message : "Nest vector argument 3 not setup "
	VecConcatenate(2, X_array, &b_hat, NULL);
	
//##### Application of the transformation A -> A_hat
	// Declaration
	Mat D_M_inv_G, Mat_array[4];
	Mat A_hat, Pmat, C_hat, G_hat;
	Mat diag_2M;//Will store 2*diagonal part of M (to approximate the Schur complement)
	Vec v;
	
	Mat_array[3]=M;//Bottom left block of A_hat

	//Extraction of the diagonal of M
	MatCreateVecs(M,NULL,&v);//v has the size of M
	MatGetDiagonal(M,v);

	//Creation of matrix 2*diag(M). Why not use MatCreateDiagonal ???
	MatCreateConstantDiagonal(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n_u, n_u, 2, &diag_2M);
	MatConvert(diag_2M,  MATAIJ, MAT_INPLACE_MATRIX, &diag_2M);
	MatDiagonalScale(diag_2M, v, NULL);//store 2*diagonal part of M
	VecReciprocal(v);//Must first check that all the coefficients are non zero
	
	// Creation of D_M_inv_G = D_M_inv*G
	MatDuplicate(G,MAT_COPY_VALUES,&D_M_inv_G);//D_M_inv_G contains G
	MatDiagonalScale( D_M_inv_G, v, NULL);//D_M_inv_G contains D_M_inv*G

	// Creation of C_hat
	MatMatMult(D,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C_hat);//C_hat contains D*D_M_inv*G
	MatAXPY(C_hat,1.0,C,SUBSET_NONZERO_PATTERN);//C_hat contains C + D*D_M_inv*G
	Mat_array[0]=C_hat;//Top left block of A_hat

	// Creation of G_hat
	MatMatMult(M,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&G_hat);//G_hat contains M*D_M_inv*G
	MatAYPX(G_hat,-1.0,G,UNKNOWN_NONZERO_PATTERN);//G_hat contains G - M*D_M_inv*G
	Mat_array[2]=G_hat;//Bottom left block of A_hat

	// Creation of -D
	MatScale(D,-1.0);
	Mat_array[1]=D;//Top right block of A_hat

	// Creation of A_hat = reordered A_input
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,Mat_array,&A_hat);

	// Creation of Pmat
	Mat_array[3]=diag_2M;
	Mat_array[1]=NULL;//Cancel top right block
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,Mat_array,&Pmat);


	// Finalisation of the preconditioner	
	IS is_U_hat,is_P_hat;
	
	ISCreateStride(PETSC_COMM_WORLD, n_u, n_p, 1, &is_U_hat);
	ISCreateStride(PETSC_COMM_WORLD, n_p, 0, 1, &is_P_hat);

//##### Calling KSP solver and monitor convergence
	KSP ksp, *subksp;
	PC pc, subpc0, subpc1;
	KSPType ksp_type = KSPFGMRES, ksp_type0, ksp_type1;
	PCType pc_type=PCFIELDSPLIT, pc_type0, pc_type1;
	int nsplit = 2;
	PCCompositeType pc_composite_type = PC_COMPOSITE_MULTIPLICATIVE;//or ADDITIVE ???

	double residu, abstol, rtol=1e-7, dtol;
	int iter, numberMaxOfIter;

	PetscPrintf(PETSC_COMM_WORLD,"Definition of the solver ...\n");
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetType(ksp, ksp_type);
	KSPSetOperators(ksp,A_hat,Pmat);
	KSPSetTolerances(ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT, PETSC_DEFAULT);
	KSPGetPC(ksp,&pc);
	PetscPrintf(PETSC_COMM_WORLD,"Setting the preconditioner %s...\n", pc_type);
	PCSetType(pc,pc_type);
	if( strcmp(pc_type , PCFIELDSPLIT)==0 ){
		PCFieldSplitSetType(pc, pc_composite_type);
		PCFieldSplitSetIS(pc, "0",is_P_hat);
		PCFieldSplitSetIS(pc, "1",is_U_hat);
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD,"Using PCILU\n");
		PCSetType(pc,PCILU);//This prec works fine in sequential
	}
	PCSetFromOptions(pc);
	PCSetUp(pc);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	PetscPrintf(PETSC_COMM_WORLD,"Solving the linear system...\n");
	KSPSolve(ksp,b_hat,X_hat);

	PCFieldSplitGetType(pc, &pc_composite_type);
	KSPGetType(ksp,&ksp_type);
		PCGetType(pc,&pc_type);
	PCFieldSplitSchurGetSubKSP( pc, &nsplit, &subksp);
	KSPGetType(subksp[0],&ksp_type0);
	KSPGetType(subksp[1],&ksp_type1);
	KSPGetPC(subksp[0], &subpc0);
	KSPGetPC(subksp[1], &subpc1);
	PCGetType( subpc0, &pc_type0);
	PCGetType( subpc1, &pc_type1);
	if(pc_composite_type==PC_COMPOSITE_MULTIPLICATIVE)
		PetscPrintf(PETSC_COMM_WORLD,"... linear system solved with ksp_type %s, pc_composite_type PC_COMPOSITE_MULTIPLICATIVE\n",ksp_type);
	else
		PetscPrintf(PETSC_COMM_WORLD,"... linear system solved with ksp_type %s, pc_composite_type %d (different from PC_COMPOSITE_MULTIPLICATIVE)\n",ksp_type,pc_composite_type);

	//Extract informations about the convergence
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);
	KSPGetIterationNumber(ksp,&iter);
	KSPGetResidualNorm( ksp, &residu);
	KSPGetTolerances( ksp, &rtol, &abstol, &dtol, &numberMaxOfIter);

	if (reason>0)
		PetscPrintf(PETSC_COMM_WORLD, "Linear system converged in %d iterations \n", iter);
	else
		PetscPrintf(PETSC_COMM_WORLD, "!!!!!!!!!!!!!!!!!! Linear system diverged  after %d iterations !!!!!!!!!!!!!!\n", iter);
		
	switch(reason){
		case 2:
		    PetscPrintf(PETSC_COMM_WORLD, "Residual 2-norm < rtol*||RHS||_2 with rtol = %e, final residual = %e\n", rtol, residu);
		    break;
		case 3:
		    PetscPrintf(PETSC_COMM_WORLD, "Residual 2-norm < atol with atol = %e, final residual = %e\n", abstol, residu);
		    break;
		case -4:
		    PetscPrintf(PETSC_COMM_WORLD, "!!!!!!! Residual 2-norm > dtol*||RHS||_2 with dtol = %e, final residual = %e !!!!!!! \n", dtol, residu);
		    break;
		case -3:
		    PetscPrintf(PETSC_COMM_WORLD, "!!!!!!! Maximum number of iterations %d reached !!!!!!! \n", numberMaxOfIter);
		    break;
		case -11:
		    PetscPrintf(PETSC_COMM_WORLD, "!!!!!!! Construction of preconditioner failed !!!!!! \n");
		    break;
		case -5:
		    PetscPrintf(PETSC_COMM_WORLD, "!!!!!!! Generic breakdown of the linear solver (Could be due to a singular matrix or preconditioner)!!!!!! \n");
		    break;
		default:
			if (reason>0)
			    PetscPrintf(PETSC_COMM_WORLD, "PETSc convergence reason %d \n", reason);
			else
			    PetscPrintf(PETSC_COMM_WORLD, "PETSc divergence reason %d \n" , reason);
		}

//##### Compute X from X_hat
	Vec X_hat_p;//Pressure components of the transformed unknown
	Vec X_hat_u;//Velocity components of the transformed unknown
	Vec X_p;//Pressure components of the main unknown
	Vec X_u;//Velocity components of the transformed unknown
	Vec X_output;
	
	VecGetSubVector( X_hat, is_P_hat, &X_hat_p);
	VecGetSubVector( X_hat, is_U_hat, &X_hat_u);

	VecDuplicate(X_hat_u,&X_u);
	VecDuplicate(X_hat_p,&X_p);
	VecCopy(X_hat_p,X_p);
	MatMult( G, X_hat_p, X_u);
	VecPointwiseMult(X_u,X_u,v);
	VecAYPX( X_u, -1, X_hat_u);

	X_array[0] = X_u;
	X_array[1] = X_p;

	//VecCreateNest( PETSC_COMM_WORLD, 2, NULL, X_array, &X_output);//This generate an error message : "Nest vector argument 3 not setup "
	VecConcatenate(2, X_array, &X_output, NULL);
	
//##### Compute the error and check it is small
	Vec X_anal_p, X_anal_u;//Pressure and velocity components of the analitic solution
	double error, error_p, error_u;
	
	VecGetSubVector( X_anal, is_P, &X_anal_p);
	VecGetSubVector( X_anal, is_U, &X_anal_u);

	VecAXPY(  X_p, -1, X_anal_p);
	VecNorm(  X_p, NORM_2, &error_p);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error p : ||X_anal_p - X_num_p|| = %e\n", error_p);
	VecAXPY(  X_u, -1, X_anal_u);
	VecNorm(  X_u, NORM_2, &error_u);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error u : ||X_anal_u - X_num_u|| = %e \n", error_u);

	VecAXPY(X_output, -1, X_anal);
	VecNorm( X_output, NORM_2, &error);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error : ||X_anal - X_num|| = %e, (remember ||X_anal||=1)\n", error);

	PetscCheck( error < 1.e-5, PETSC_COMM_WORLD, ierr, "Linear system did not return accurate solution. Error is too high\n");
	
//##### Cleaning of the memory
	MatDestroy(&A_input);
	MatDestroy(&A_hat);
	MatDestroy(&Pmat);
	MatDestroy(&M);
	MatDestroy(&G_hat);
	MatDestroy(&C_hat);
	MatDestroy(&D);	
	MatDestroy(&G);
	MatDestroy(&C);
	MatDestroy(&D_M_inv_G);
	MatDestroy(&diag_2M);

	VecDestroy(&b_input);
	VecDestroy(&b_hat);
	VecDestroy(&X_hat);
	VecDestroy(&X_anal);
	VecDestroy(&v);

	ISDestroy(&is_U);
	ISDestroy(&is_P);
	ISDestroy(&is_U_hat);
	ISDestroy(&is_P_hat);

	KSPDestroy(&ksp);
	PetscFree(subksp);
	
	PetscFinalize();
	return ierr;
}
