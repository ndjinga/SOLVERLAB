static char help[] = "Read a PETSc matrix from a file -f0 <input file>\n Parameters : \n -f0 : matrix fileName \n -nU :number of velocity lines \n -nP : number of pressure lines \n -mat_type : PETSc matrix type \n";

#include <petscis.h>
#include <petscksp.h>

int main( int argc, char **args ){
	PetscInitialize(&argc,&args, (char*)0,help);
	PetscMPIInt    size;        /* size of communicator */
	PetscMPIInt    rank;        /* processor rank */
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	PetscErrorCode ierr=0;

//##### Load the matrix A in the file given in the command line
	char file[1][PETSC_MAX_PATH_LEN], mat_type[256]; // File to load, matrix type
	PetscViewer viewer;
	Mat A;
	PetscBool flg;

	PetscOptionsGetString(NULL,NULL,"-f0",file[0],PETSC_MAX_PATH_LEN,&flg);
	PetscStrcpy(mat_type,MATAIJ);
	PetscOptionsGetString(NULL,NULL,"-mat_type",mat_type,sizeof(mat_type),NULL);

	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix type %s from file %s on %d processor(s)...\n", mat_type, file[0], size);	
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);	
	PetscViewerSetType(viewer,PETSCVIEWERBINARY);//Use PETSCVIEWERHDF5 for better parallel performance
	PetscViewerFileSetMode(viewer,FILE_MODE_READ);
	PetscViewerFileSetName(viewer,file[0]);
	
	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetType(A,mat_type);
	MatLoad(A,viewer);
	PetscViewerDestroy(&viewer);
	PetscPrintf(PETSC_COMM_WORLD,"... matrix Loaded \n");	
	

//####	Decompose the matrix A into 4 blocks M, G, D, C
	Mat M, G, D, C;
	PetscInt nrows, ncolumns;//Total number of rows and columns of A
	PetscInt irow_min, irow_max;//min and max indices of rows stored locally on this process
	PetscInt n_u, n_p, n;//Total number of velocity and pressure lines. n = n_u+ n_p
	IS is_U,is_P;
	PC pc;

	PetscOptionsGetInt(NULL,NULL,"-nU",&n_u,NULL);
	PetscOptionsGetInt(NULL,NULL,"-nP",&n_p,NULL);
	n=n_u+n_p;
	MatGetOwnershipRange( A, &irow_min, &irow_max);
	MatGetSize( A, &nrows, &ncolumns);
	int nb_pressure_lines = irow_max >= n_u ? irow_max - n_u : 0;
	int nb_velocity_lines = irow_min <= n_u ? n_u - irow_min : 0;
	PetscInt i_p[nb_pressure_lines],i_u[nb_velocity_lines];

	PetscCheck( nrows == ncolumns, PETSC_COMM_WORLD, ierr, "Matrix is not squared !!!\n");
	PetscCheck( n == ncolumns, PETSC_COMM_WORLD, ierr, "Inconsistent data : the matrix has %d lines but only %d velocity lines and %d pressure lines declared\n", ncolumns, n_u,n_p);
	PetscPrintf(PETSC_COMM_WORLD,"The matrix has %d lines : %d velocity lines and %d pressure lines\n", n, n_u,n_p);
	PetscPrintf(PETSC_COMM_SELF,"irow_min = %d, irow_max = %d \n", irow_min, irow_max);
	
	for (int i = n_u;i<irow_max;i++){
		i_p[i-n_u]=i;
	}
	for (int i=irow_min;i<n_u;i++){
		i_u[i]=i;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Extraction of the 4 blocks \n");
	ISCreateGeneral(PETSC_COMM_WORLD,n_u,(const PetscInt *) i_u,PETSC_OWN_POINTER,&is_U);
	ISCreateGeneral(PETSC_COMM_WORLD,n_p,(const PetscInt *) i_p,PETSC_OWN_POINTER,&is_P);
	MatCreateSubMatrix(A,is_U, is_U,MAT_INITIAL_MATRIX,&M);
	MatCreateSubMatrix(A,is_U, is_P,MAT_INITIAL_MATRIX,&G);
	MatCreateSubMatrix(A,is_P, is_U,MAT_INITIAL_MATRIX,&D);
	MatCreateSubMatrix(A,is_P, is_P,MAT_INITIAL_MATRIX,&C);
	MatDestroy(&A);//Early destruction since A is a sequential matrix stored on processed 0
	PetscPrintf(PETSC_COMM_WORLD,"... end of extraction\n");

//##### Application of the transformation A -> A_hat
	// Declaration
	Mat D_M_inv_G,array[4];
	Mat A_hat,Pmat, S_hat, C_hat,G_hat;
	Vec v;
	array[3]=M;

	//Extraction of the diagonal of M
	VecCreate(PETSC_COMM_WORLD,&v);
	VecSetSizes(v,n_u,PETSC_DECIDE);
	VecSetFromOptions(v);
	VecSetUp(v);
	MatGetDiagonal(M,v);
	VecReciprocal(v);
	
	// Creation of D_M_inv_G = D_M_inv*G
	MatDuplicate(G,MAT_COPY_VALUES,&D_M_inv_G);//D_M_inv_G contains G
	MatDiagonalScale( D_M_inv_G, v, NULL);//D_M_inv_G contains D_M_inv*G

	// Creation of C_hat
	MatMatMult(D,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C_hat);//C_hat contains D*D_M_inv*G
	MatAXPY(C_hat,1.0,C,SUBSET_NONZERO_PATTERN);//C_hat contains C + D*D_M_inv*G
	array[0]=C_hat;

	// Creation of G_hat
	MatMatMult(M,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&G_hat);//G_hat contains M*D_M_inv*G
	MatAYPX(G_hat,-1.0,G,UNKNOWN_NONZERO_PATTERN);//G_hat contains G - M*D_M_inv*G
	array[2]=G_hat;

	// Creation of -D
	MatScale(D,-1.0);
	array[1]=D;

	// Creation of A_hat = reordered A
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,array,&A_hat);
	MatConvert(A_hat,MATAIJ,MAT_INPLACE_MATRIX,&A_hat);

	// Creation of S_hat=2M: spectral equivalent to the Schur complement
	MatDuplicate(M,MAT_COPY_VALUES,&S_hat);
	MatScale(S_hat,2.0);

	// Finalisation of the preconditioner	
	PetscInt i_p_hat[n_p],i_u_hat[n_u];
	IS is_U_hat,is_P_hat;
	
	for (int i=0;i<n_p;i++){
		i_p_hat[i]=i;
	}
	for (int i=n_p;i<n;i++){
		i_u_hat[i-n_p]=i;
	}

	ISCreateGeneral(PETSC_COMM_WORLD,n_p,(const PetscInt *) i_p_hat,PETSC_OWN_POINTER,&is_P_hat);
	ISCreateGeneral(PETSC_COMM_WORLD,n_u,(const PetscInt *) i_u_hat,PETSC_OWN_POINTER,&is_U_hat);

	// Creation of Pmat
	array[1]=NULL;
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,array,&Pmat);
	MatConvert(Pmat,MATAIJ,MAT_INPLACE_MATRIX,&Pmat);


//##### Definition of the right hand side to test the preconditioner
	Vec b,u,x_anal;
	KSP ksp;
	PetscScalar y[n_p];
	double residu, abstol, rtol, norm_x_anal;
	int iter;

	PetscPrintf(PETSC_COMM_WORLD,"Creation of the RHS, exact and numerical solution vectors...\n");
	VecCreate(PETSC_COMM_WORLD,&b);
	VecSetSizes(b,n_u+n_p,PETSC_DECIDE);
	VecSetFromOptions(b);

	VecDuplicate(b,&x_anal);//x_anal will store the exact solution
	VecDuplicate(b,&u);//u will store the numerical solution
	VecSet(x_anal,0.0);

	for (int i = n_u;i<n;i++){
		y[i-n_u]=1.0;
	}
	VecSetValues(x_anal,n_p,i_p,y,INSERT_VALUES);
	MatMult( A_hat, x_anal, b);
	PetscPrintf(PETSC_COMM_WORLD,"... vectors created \n");	

	PetscPrintf(PETSC_COMM_WORLD,"Definition of the KSP solver to test the preconditioner...\n");
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetType(ksp,KSPFGMRES);
	KSPSetOperators(ksp,A_hat,Pmat);
	KSPGetPC(ksp,&pc);
	PetscPrintf(PETSC_COMM_WORLD,"Setting the preconditioner...\n");
	PCSetType(pc,PCFIELDSPLIT);
	PCFieldSplitSetIS(pc, "0",is_P_hat);
	PCFieldSplitSetIS(pc, "1",is_U_hat);
	PCFieldSplitSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
//	PCSetType(pc,PCILU);
	PCSetFromOptions(pc);
	PCSetUp(pc);
	KSPSetFromOptions(ksp);
	PetscPrintf(PETSC_COMM_WORLD,"Solving the linear system...\n");
	KSPSolve(ksp,b,u);
	//Extract informations about the convergence
	KSPGetIterationNumber(ksp,&iter);
	KSPGetResidualNorm( ksp, &residu);
	KSPGetTolerances( ksp, &rtol, &abstol, NULL, NULL);
	PetscPrintf(PETSC_COMM_WORLD,"... linear system solved in %d iterations, final residual %e, relative tolerance %e, absolute tolerance %e\n", iter, residu, rtol, abstol);

//	MatView(M,PETSC_VIEWER_STDOUT_WORLD);
//	VecView(u,PETSC_VIEWER_STDOUT_WORLD);
	
//##### Compute the error and check it is small
	double error = 0.;

	VecAXPY(u, -1, x_anal);
	VecNorm( u, NORM_2, &error);
	VecNorm( x_anal, NORM_2, &norm_x_anal);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error : ||x_anal - x_num|| = %e, ||x_anal - x_num||/||x_anal|| = %e\n", error, error/norm_x_anal);

	PetscCheck( error/norm_x_anal < 1.e-5, PETSC_COMM_WORLD, ierr, "Linear system did not return accurate solution. Error is too high\n");
	
	// Cleaning of the code
	MatDestroy(&M);
	MatDestroy(&S_hat);
	MatDestroy(&G_hat);
	MatDestroy(&C_hat);
	MatDestroy(&D);	
	MatDestroy(&G);
	MatDestroy(&C);
	VecDestroy(&b);
	VecDestroy(&u);
	KSPDestroy(&ksp);

	PetscFinalize();
	return ierr;
}
