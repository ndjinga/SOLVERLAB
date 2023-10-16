static char help[] = "Read a PETSc matrix from a file -f0 <input file>";

#include <petscis.h>
#include <petscksp.h>

int main( int argc, char **args ){
	//PetscErrorCode ierr;
	//Defintion of the main objects and intialization of the problem
	char file[1][PETSC_MAX_PATH_LEN]; // File to load
	PetscViewer viewer;
	Mat A, S_hat, M,G,D,C,C_hat,G_hat;
	PetscBool flg;
	Vec b,u;
	PetscInt n_u,n_p,n,iter;
	IS is_U,is_P;
	KSP ksp;
	PC pc;
	//int rank, size;
//	PetscObjectSetName((PetscObject)A, "A");
//	PetscObjectSetName((PetscObject)v, "v");
	PetscInitialize(&argc,&args, (char*)0,help);
	PetscOptionsGetString(NULL,NULL,"-f0",file[0],PETSC_MAX_PATH_LEN,&flg);
	PetscOptionsGetInt(NULL,NULL,"-nU",&n_u,NULL);
	PetscOptionsGetInt(NULL,NULL,"-nP",&n_p,NULL);
	n=n_u+n_p;
	PetscScalar y[n_p];
	PetscInt i_p[n_p],i_u[n_u];

	// Loading of the context of the matrix
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);	
	PetscViewerSetType(viewer,PETSCVIEWERBINARY);
	PetscViewerFileSetMode(viewer,FILE_MODE_READ);
	PetscViewerFileSetName(viewer,file[0]);
	
//	PetscPrintf(PETSC_COMM_WORLD,"Loading started...\n");	

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetType(A,MATSEQAIJ);
	MatLoad(A,viewer);
	PetscViewerDestroy(&viewer);
	

	// Definition of the right hand side
//	PetscPrintf(PETSC_COMM_WORLD,"Definition of the vectors...\n");
	VecCreate(PETSC_COMM_WORLD,&b);
	//PetscObjectSetName((PetscObject)b,"RHS");
	VecSetSizes(b,n,PETSC_DECIDE);
	VecSetFromOptions(b);
	VecSetUp(b);

	VecDuplicate(b,&u);
	//PetscObjectSetName((PetscObject)u,"Solution");
	VecSet(b,0.0);
	// Definition of the transformation
	for (int i = n_u;i<n;i++){
		i_p[i-n_u]=i;
		y[i-n_u]=1.0;
	}
	for (int i=0;i<n_u;i++){
		i_u[i]=i;
	}
	VecSetValues(b,n_p,i_p,y,INSERT_VALUES);

// 	Decomposition into blocks	
	ISCreateGeneral(PETSC_COMM_WORLD,n_u,(const PetscInt *) i_u,PETSC_OWN_POINTER,&is_U);
	ISCreateGeneral(PETSC_COMM_WORLD,n_p,(const PetscInt *) i_p,PETSC_OWN_POINTER,&is_P);
	MatCreateSubMatrix(A,is_U, is_U,MAT_INITIAL_MATRIX,&M);
	MatCreateSubMatrix(A,is_U, is_P,MAT_INITIAL_MATRIX,&G);
	MatCreateSubMatrix(A,is_P, is_U,MAT_INITIAL_MATRIX,&D);
	MatCreateSubMatrix(A,is_P, is_P,MAT_INITIAL_MATRIX,&C);
	MatDestroy(&A);

//	Application of the transformation
	// Declaration
	Mat D_M_inv,I_u,array[4],A_hat,Pmat;
	Vec v,v_I_u,v_D_inv;
	array[3]=M;

	// Creation of identity of n_u
	VecCreate(PETSC_COMM_WORLD,&v_I_u);
	VecSetSizes(v_I_u,n_u,PETSC_DECIDE);
	VecSetFromOptions(v_I_u);
	VecSetUp(v_I_u);
	VecSet(v_I_u,1.0);
	MatCreate(PETSC_COMM_WORLD,&I_u);
	MatSetSizes(I_u,PETSC_DECIDE,PETSC_DECIDE,n_u,n_u);
	MatSetFromOptions(I_u);
	MatSetUp(I_u);
	MatDiagonalSet(I_u,v_I_u,INSERT_VALUES);
	
	
	// Creation of D_M_inv
	MatDuplicate(I_u,MAT_SHARE_NONZERO_PATTERN,&D_M_inv);
	MatZeroEntries(D_M_inv);
	VecCreate(PETSC_COMM_WORLD,&v);
	VecSetSizes(v,n_u,PETSC_DECIDE);
	VecSetFromOptions(v);
	VecSetUp(v);
	VecDuplicate(v,&v_D_inv);
	MatGetDiagonal(M,v);
	VecPointwiseDivide(v_D_inv,v_I_u,v);
	MatDiagonalSet(D_M_inv,v_D_inv,INSERT_VALUES);

	// Creation of C_hat
	MatMatMatMult(D,D_M_inv,G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C_hat);
	MatAXPY(C_hat,1.0,C,SUBSET_NONZERO_PATTERN);
	array[0]=C_hat;

	// Creation of G_hat
	MatMatMatMult(M,D_M_inv,G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&G_hat);
	MatAYPX(G_hat,-1.0,G,UNKNOWN_NONZERO_PATTERN);
	array[2]=G_hat;

	// Creation of -D
	MatScale(D,-1.0);
	array[1]=D;

	// Creation of reordered A
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,array,&A_hat);
//	MatConvert(A_hat,MATAIJ,MAT_INPLACE_MATRIX,&A_hat);

	// Creation of S_hat=2M: spectral equivalent to the Schur complement
//	MatDuplicate(M,MAT_COPY_VALUES,&S_hat);
//	MatScale(S_hat,2.0);
	// Creation of S_hat=2diag(M)
	MatDuplicate(D_M_inv,MAT_SHARE_NONZERO_PATTERN,&S_hat);
	MatDiagonalSet(S_hat,v,INSERT_VALUES);
	MatScale(S_hat,2.0);



	// Check of the values of the right hand side
//	PetscPrintf(PETSC_COMM_WORLD,"Check of the value of the vector...");
//	VecView(b,PETSC_VIEWER_STDOUT_WORLD);
//	Preparation of the preconditioner	
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
//	MatConvert(Pmat,MATAIJ,MAT_INPLACE_MATRIX,&Pmat);


//	PetscPrintf(PETSC_COMM_WORLD,"Definition of the solver...\n");
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetType(ksp,KSPFGMRES);
	KSPSetOperators(ksp,A_hat,A_hat);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCFIELDSPLIT);
	PCFieldSplitSetIS(pc, "0",is_P_hat);
	PCFieldSplitSetIS(pc, "1",is_U_hat);
	PCFieldSplitSetSchurPre(pc,PC_FIELDSPLIT_SCHUR_PRE_USER,S_hat);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	KSPSolve(ksp,b,u);
	KSPGetIterationNumber(ksp,&iter);



	PetscPrintf(PETSC_COMM_WORLD,"Number of iterations : %d...\n",iter);


	
	// Cleaning of the code
	MatDestroy(&M);
	MatDestroy(&S_hat);
	MatDestroy(&G_hat);
	MatDestroy(&C_hat);
	MatDestroy(&D);	
	MatDestroy(&G);
	MatDestroy(&C);
	MatDestroy(&D_M_inv);
	MatDestroy(&I_u);
	MatDestroy(&A_hat);
	MatDestroy(&Pmat);
	VecDestroy(&b);
	VecDestroy(&u);
	VecDestroy(&v);
	VecDestroy(&v_I_u);
	VecDestroy(&v_D_inv);
	KSPDestroy(&ksp);

	PetscFinalize();
	return 0;
}
