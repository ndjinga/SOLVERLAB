static char help[] = "Read a PETSc matrix from a file -f0 <input file>\n Parameters : \n -f0 : matrix fileName \n -nU :number of velocity lines \n -nP : number of pressure lines \n -nprocs : number of processors \n -mat_type : PETSc matrix type";

#include <petscis.h>
#include <petscksp.h>

int main( int argc, char **args ){
	//PetscErrorCode ierr;
	//Defintion of the main objects and intialization of the problem
	char file[1][PETSC_MAX_PATH_LEN], mat_type[256]; // File to load, matrix type
	PetscViewer viewer;
	Mat A, S_hat, M,G,D,C,C_hat,G_hat;
	PetscBool flg;
	Vec b,u;
	PetscInt n_u,n_p,n,iter, nprocs=1;
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
	PetscOptionsGetInt(NULL,NULL,"-nprocs",&nprocs,NULL);
	PetscStrcpy(mat_type,MATAIJ);
	PetscOptionsGetString(NULL,NULL,"-mat_type",mat_type,sizeof(mat_type),NULL);
	n=n_u+n_p;
	PetscScalar y[n_p];
	PetscInt i_p[n_p],i_u[n_u];

	// Loading of the context of the matrix
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix type %s from file %s on %d processor(s)...\n", mat_type, file[0], nprocs);	
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);	
	PetscViewerSetType(viewer,PETSCVIEWERBINARY);//Use PETSCVIEWERHDF5 for better parallel performance
	PetscViewerFileSetMode(viewer,FILE_MODE_READ);
	PetscViewerFileSetName(viewer,file[0]);
	
	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetType(A,mat_type);
	MatLoad(A,viewer);
	PetscViewerDestroy(&viewer);
	PetscPrintf(PETSC_COMM_WORLD,"... matrix Loaded \n");	
	

	// Definition of the right hand side
	PetscPrintf(PETSC_COMM_WORLD,"Creation of the vectors...\n");
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
	PetscPrintf(PETSC_COMM_WORLD,"... vectors created \n");	

// 	Decomposition into blocks	
	PetscPrintf(PETSC_COMM_WORLD,"Extraction of the 4 blocks \n");
	ISCreateGeneral(PETSC_COMM_WORLD,n_u,(const PetscInt *) i_u,PETSC_OWN_POINTER,&is_U);
	ISCreateGeneral(PETSC_COMM_WORLD,n_p,(const PetscInt *) i_p,PETSC_OWN_POINTER,&is_P);
	MatCreateSubMatrix(A,is_U, is_U,MAT_INITIAL_MATRIX,&M);
	MatCreateSubMatrix(A,is_U, is_P,MAT_INITIAL_MATRIX,&G);
	MatCreateSubMatrix(A,is_P, is_U,MAT_INITIAL_MATRIX,&D);
	MatCreateSubMatrix(A,is_P, is_P,MAT_INITIAL_MATRIX,&C);
	PetscPrintf(PETSC_COMM_WORLD,"... end of extraction\n");

//	Application of the transformation
	// Declaration
	Mat D_M_inv_G,array[4],A_hat,Pmat;
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
	MatAYPX(G_hat,-1.0,G,UNKNOWN_NONZERO_PATTERN);//G_hat contains -G + M*D_M_inv*G
	array[2]=G_hat;

	// Creation of -D
	MatScale(D,-1.0);
	array[1]=D;

	// Creation of reordered A
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,array,&A_hat);
	MatConvert(A_hat,MATAIJ,MAT_INPLACE_MATRIX,&A_hat);

	// Creation of S_hat=2M: spectral equivalent to the Schur complement
	MatDuplicate(M,MAT_COPY_VALUES,&S_hat);
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
	MatConvert(Pmat,MATAIJ,MAT_INPLACE_MATRIX,&Pmat);


	PetscPrintf(PETSC_COMM_WORLD,"Definition of the KSP solver...\n");
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetType(ksp,KSPFGMRES);
	KSPSetOperators(ksp,A_hat,Pmat);
	KSPGetPC(ksp,&pc);
	PetscPrintf(PETSC_COMM_WORLD,"Setting thepreconditioner...\n");
	PCSetType(pc,PCFIELDSPLIT);
	PCFieldSplitSetIS(pc, "0",is_P_hat);
	PCFieldSplitSetIS(pc, "1",is_U_hat);
	PCFieldSplitSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
//	PCSetType(pc,PCILU);
	PCSetFromOptions(pc);
	PCSetUp(pc);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	PetscPrintf(PETSC_COMM_WORLD,"Solving the linear system...\n");
	KSPSolve(ksp,b,u);
	KSPGetIterationNumber(ksp,&iter);
	PetscPrintf(PETSC_COMM_WORLD,"... linear system solved\n");





//	PetscPrintf(PETSC_COMM_WORLD,"Print of the solution...\n");
//	PetscPrintf(PETSC_COMM_WORLD,"Print M...\n");
//	MatView(A_hat,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"Print G...\n");
//	MatView(Pmat,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"Print D...\n");
//	MatView(D,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"Print C...\n");
//	MatView(G_hat,PETSC_VIEWER_STDOUT_WORLD);
//	MatView(M,PETSC_VIEWER_STDOUT_WORLD);
	VecView(u,PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"Number of iterations : %d...\n",iter);


	
	// Cleaning of the code
	MatDestroy(&A);
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
	return 0;
}
