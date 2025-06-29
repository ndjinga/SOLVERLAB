static char help[] = "Read a PETSc matrix from a file Parameters : \n -f0 : matrix fileName (mandatory) \n -nU : global number of velocity lines (assumption of contiguous velocity line numbers) \n -nP : global number of pressure lines (assumption of contiguous pressure line numbers) \n -mat_type : PETSc matrix type (optional) \n";

/*************************************************************************************************/
/* Parallel implementation of a new preconditioner for the linear system A_{input} X_{output} = b_{input} */
/* The coupling between two variables u and p displays a saddle point structure  ( M G \\ D C )  */
/*                                                                                               */
/* Input  : - Matrix A_{input}    (system matrix, loaded from a file)                            */
/*          - Vector b_{input}    (right hand side, made up for testing)                         */
/* Output : - Vector X_{output}   (unknown vector, to be determined)                             */
/*                                                                                               */
/* Auxilliary variables : - A_hat (transformed matrix)                                           */
/*                        - X_hat (unknown of the transformed system)                            */
/*                        - b_hat (RHS of the transformed system)                                */
/*                        - Pmat  (preconditioning matrix)                                       */
/*                        - M submatrix u-u of A_{input}                                         */
/*                        - G submatrix u-p of A_{input}                                                  */
/*                        - D submatrix p-u of A_{input}                                                  */
/*                        - C submatrix p-p of A_{input}                                                  */
/*                        - R submatrix (neither p nor u lines) - (neither p nor u columns) of A_{input}  */
/*                                                                                               */
/*                                 * M G X *                                                     */
/*                        A     = *  D C X  *                                                    */
/*                                 * X X R *                                                     */
/*                                                                                               */
/*                                 * C_hat -D X *                                                */
/*                        A_hat = *  G_hat  M X  *                                               */
/*                                 * X      X R *                                                */
/*                                                                                               */
/*                                 * C_hat 0          0 *                                        */
/*                        Pmat  = *  G_hat 2 diag(M)  0  *                                       */
/*                                 * X     X          R *                                        */
/*                                                                                               */
/*************************************************************************************************/

#include <petscis.h>
#include <petscksp.h>

int main( int argc, char **args ){
  /*
    Every PETSc routine should begin with the PetscInitialize() routine.
    argc, argv - These command line arguments are taken to extract the options
                 supplied to PETSc and options supplied to MPI.
    help       - When PETSc executable is invoked with the option -help,
                 it prints the various options that can be applied at
                 runtime.  The user can use the "help" variable place
                 additional help messages in this printout.
  */
 	PetscInitialize(&argc,&args, (char*)0,help);
	PetscMPIInt    size;        /* size of communicator */
	PetscMPIInt    rank;        /* rank of processor */
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);/* Get rank of processor */
	MPI_Comm_size(PETSC_COMM_WORLD,&size);/* Get size of communicator */
	PetscErrorCode ierr=0;

//##### Load the matrix A contained in the file given in the command line
	char file[1][PETSC_MAX_PATH_LEN], mat_type[256]; // File to load, matrix type
	PetscViewer viewer;
	Mat A_input;
	PetscBool setFileName;//variable to store whether the file containing the matrix was found

	PetscOptionsGetString(NULL,NULL,"-f0",file[0],PETSC_MAX_PATH_LEN,&setFileName);//Get the file name from command line
	if( !setFileName )//Check file name was found
	{
          PetscErrorPrintf(PETSC_COMM_WORLD," Error : no file name provided. Use the keyword -f0 followed by the file name in the command line.\n");	
          return  PETSC_ERR_ARG_WRONG;
        }
	PetscStrcpy(mat_type,MATAIJ);// Default value for PETSc Matrix type
	PetscOptionsGetString(NULL,NULL,"-mat_type",mat_type,sizeof(mat_type),NULL);//Changes the value of mat_type if given in the command line

	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix type %s from file %s on %d processor(s)...\n", mat_type, file[0], size);	
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);	
	PetscViewerSetType(viewer,PETSCVIEWERBINARY);//Use PETSCVIEWERHDF5 for better parallel performance
	PetscViewerFileSetMode(viewer,FILE_MODE_READ);
	PetscViewerFileSetName(viewer,file[0]);
	
	MatCreate(PETSC_COMM_WORLD, &A_input);
	MatSetType(A_input,mat_type);
	MatLoad(A_input,viewer);
	PetscViewerDestroy(&viewer);
	PetscPrintf(PETSC_COMM_WORLD,"... matrix Loaded \n");	

//####	Decompose the matrix A_input into 9=3x3 blocks and extract the 5 blocks M, G, D, C and Remaining_Diagonal_Block
	Mat M, G, D, C, Remaining_Diagonal_Block;
	PetscInt nrows, ncolumns;//Total number of rows and columns of A_input
	PetscInt irow_min, irow_max;//min and max indices of rows stored locally on this process
	IS is_U,is_P;
	PetscInt n_u, n_p, n_neither_U_nor_P;//Total number of velocity (n_u), pressure (n_p) and remaining (n_neither_U_nor_P) lines. n=matrix size = n_u+ n_p+ n_neither_U_nor_P
	PetscBool setNbU, setNbP;

	MatGetOwnershipRange( A_input, &irow_min, &irow_max);
	MatGetSize( A_input, &nrows, &ncolumns);
	PetscCheck( nrows == ncolumns, PETSC_COMM_WORLD, ierr, "Matrix is not square !!!\n");
	PetscOptionsGetInt(NULL,NULL,"-nU",&n_u,setNbU);
	PetscOptionsGetInt(NULL,NULL,"-nP",&n_p,setNbP);

	if( setNbU && setNbP ) //build is_U and is_P as a stride
	{	
	n_neither_U_nor_P = nrows - n_u - n_p
	//Contiguous velocity lines followed by contiguous pressure lines : pressure indices must come after the velocity indices
	PetscInt min_pressure_lines = irow_min <= n_u ? n_u : irow_min;//max(irow_min, n_u)
	PetscInt max_velocity_lines = irow_max >= n_u ? n_u : irow_max;//min(irow_max, n_u)
	//velocity (resp. pressure) indices are assumed to be consecutive, and nu+np+ n_neither_U_nor_P = irow_max - irow_min
	PetscInt nb_pressure_lines = irow_max >= n_u ? min(irow_max, n_neither_U_nor_P) - min(min_pressure_lines, n_neither_U_nor_P) : 0;
	PetscInt nb_velocity_lines = irow_min <= n_u ? max_velocity_lines - irow_min : 0;
	ISCreateStride(PETSC_COMM_WORLD, nb_velocity_lines, max_velocity_lines - nb_velocity_lines, 1, &is_U);
	ISCreateStride(PETSC_COMM_WORLD, nb_pressure_lines, min_pressure_lines                    , 1, &is_P);
	PetscPrintf(PETSC_COMM_WORLD,"-nU and -nP set, so contiguous velocity line numbers followed by contiguous pressure line numbers);	
	PetscPrintf(PETSC_COMM_SELF,"Process %d local rows : irow_min = %d, irow_max = %d, min_pressure_lines = %d, max_velocity_lines = %d, local nb_pressure_lines = %d, local nb_velocity_lines = %d \n", rank, irow_min, irow_max, irow_min, irow_max, min_pressure_lines, max_velocity_lines, nb_pressure_lines, nb_velocity_lines);
	}
	else
	{
	n_intersect;//Number of indices belonging to both is_U and is_P
	IS is_intersect, is_union, is_neither_U_nor_P;
	
	ISIntersect( is_U, is_P, &is_intersect);
	ISExpand( is_U, is_P, &is_union);
	ISComplement(is_union, 0, n_rows, is_neither_U_nor_P);

	ISGetSize(is_U, &n_u);//Total number of velocity lines.
	ISGetSize(is_P, &n_p);//Total number of pressure lines.
	ISGetSize(is_neither_U_nor_P, &n_neither_U_nor_P);//Total number of remaining lines.
	ISGetSize(is_intersect, &n_intersect);//n_intersect should equal zero
	
	PetscCheck( n_intersect == 0, PETSC_COMM_WORLD, ierr, "is_U and is_P should not intersect (common row indices for pressure and velocity) !!!\n");	
	PetscCheck( n_u+n_p + n_neither_U_nor_P == ncolumns, PETSC_COMM_WORLD, ierr, "Inconsistent data : the matrix has %d lines but %d velocity lines, %d pressure lines and %d remaining lines declared : n_u+n_p +n_neither_U_nor_P=%d, is not equal to the number of lines %d\n", ncolumns, n_u,n_p,n_neither_U_nor_P,n_u+n_p +n_neither_U_nor_P,ncolumns);

	PetscPrintf(PETSC_COMM_WORLD,"-nU and -nP not set (isU and isP set ?) so possibly non contiguous velocity and pressure lines);	
	PetscPrintf(PETSC_COMM_SELF,"Process %d local rows : irow_min = %d, irow_max = %d\n", rank, irow_min, irow_max);
	}

	PetscPrintf(PETSC_COMM_WORLD,"The global matrix has %d lines : %d velocity lines, %d pressure and %d remaining lines\n", n, n_u,n_p,n_neither_U_nor_P);
	
	PetscPrintf(PETSC_COMM_WORLD,"Extraction of the 4 blocks M,G,D,C,R :\n M G *\n D C *\n * * *\n");
	
	MatCreateSubMatrix(A_input,is_U, is_U,MAT_INITIAL_MATRIX,&M);
	MatCreateSubMatrix(A_input,is_U, is_P,MAT_INITIAL_MATRIX,&G);
	MatCreateSubMatrix(A_input,is_P, is_U,MAT_INITIAL_MATRIX,&D);
	MatCreateSubMatrix(A_input,is_P, is_P,MAT_INITIAL_MATRIX,&C);
	PetscPrintf(PETSC_COMM_WORLD,"... end of extraction\n");

	//#Display some informations about the four blocs
	int size1, size2;
	MatGetSize(M, &size1,&size2);
	PetscPrintf(PETSC_COMM_WORLD,"Size of M : %d,%d\n", size1,size2);
	MatGetSize(C, &size1,&size2);
	PetscPrintf(PETSC_COMM_WORLD,"Size of C : %d,%d\n", size1,size2);
	MatGetSize(G, &size1,&size2);
	PetscPrintf(PETSC_COMM_WORLD,"Size of G : %d,%d\n", size1,size2);
	MatGetSize(D, &size1,&size2);
	PetscPrintf(PETSC_COMM_WORLD,"Size of D : %d,%d\n", size1,size2);

//##### Definition of the right hand side to test the preconditioner
	Vec b_input, X_hat, X_anal;
	PetscScalar y[nb_pressure_lines];//To store the values
	PetscInt  i_p[nb_pressure_lines];//To store the indices

	PetscPrintf(PETSC_COMM_WORLD,"Creation of the RHS, exact and numerical solution vectors...\n");
	MatCreateVecs( A_input,&b_input,&X_anal );// parallel distribution of vectors should optimise the computation A_input*X_anal=b_input
	VecDuplicate(X_anal, &X_hat);// X_hat will store the numerical solution of the transformed system

	VecSet(X_anal,0.0);

	for (int i = min_pressure_lines;i<min(irow_max, n_neither_U_nor_P);i++){
		y[i-n_u]=1.0/i;//valeur second membre à imposer ici
		i_p[i-n_u]=i;
	}
	
	VecSetValues(X_anal,nb_pressure_lines,i_p,y,INSERT_VALUES);
	VecAssemblyBegin(X_anal );
	VecAssemblyEnd(  X_anal );
	VecNormalize(    X_anal, NULL);

	MatMult( A_input, X_anal, b_input);
	PetscPrintf(PETSC_COMM_WORLD,"... vectors created \n");	
	MatDestroy(&A_input);//Early destruction since A_input is a sequential matrix stored on process 0

//##### Create the matrix used for the change of variable
	Mat Id_M, Id_C, Id_R, D_M_inv_G, Pmat, invPmat, PMat_array[9], invPMat_array[9];
	Vec v_M, v_C, v_R;
	Vec v_redistributed;//different distribution of coefficients among the processors for v_M
	VecScatter scat;//tool to redistribute a vector on the processors
	IS is_to, is_from;
	
	MatCreateVecs(M,NULL,&v_M);//v_M has the parallel distribution of M
	MatCreateVecs(C,NULL,&v_C);//v_C has the parallel distribution of C
	MatCreateVecs(R,NULL,&v_R);//v_R has the parallel distribution of R

	VecSet(v_M,1.0);
	VecSet(v_C,1.0);
	VecSet(v_R,1.0);

	MatDuplicate(M, MAT_DO_NOT_COPY_VALUES, &Id_M);
	MatDuplicate(C, MAT_DO_NOT_COPY_VALUES, &Id_C);
	MatDuplicate(R, MAT_DO_NOT_COPY_VALUES, &Id_R);

	MatEliminateZeros(Id_M, PETSC_TRUE);
	MatEliminateZeros(Id_C, PETSC_TRUE);
	MatEliminateZeros(Id_R, PETSC_TRUE);
	
	MatDiagonalSet(Id_M, v_M, INSERT_VALUES);
	MatDiagonalSet(Id_C, v_C, INSERT_VALUES);
	MatDiagonalSet(Id_R, v_R, INSERT_VALUES);

	MatGetDiagonal(M,v_M);
	VecReciprocal(v_M);
	
	// Creation of D_M_inv_G = D_M_inv*G
	MatDuplicate(G,MAT_COPY_VALUES,&D_M_inv_G);//D_M_inv_G contains G
	MatCreateVecs(D_M_inv_G,NULL,&v_redistributed);//v_redistributed has the parallel distribution of D_M_inv_G
	PetscInt col_min, col_max;
	VecGetOwnershipRange(v_M,&col_min,&col_max);
	ISCreateStride(PETSC_COMM_WORLD, col_max-col_min, col_min, 1, &is_from);
	VecGetOwnershipRange(v_redistributed,&col_min,&col_max);
	ISCreateStride(PETSC_COMM_WORLD, col_max-col_min, col_min, 1, &is_to);
	VecScatterCreate(v_M,is_from,v_redistributed,is_to,&scat);
	VecScatterBegin(scat, v_M, v_redistributed,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(  scat, v_M, v_redistributed,INSERT_VALUES,SCATTER_FORWARD);
	MatDiagonalScale( D_M_inv_G, v_redistributed, NULL);//D_M_inv_G contains D_M_inv*G

	PMat_array[0]=Id_M;
	PMat_array[1]=D_M_inv_G;
	PMat_array[2]=NULL;
	PMat_array[3]=NULL;
	PMat_array[4]=Id_C;
	PMat_array[5]=NULL;
	PMat_array[6]=NULL;
	PMat_array[7]=NULL;
	PMat_array[8]=Id_R;

	invPMat_array[0]=Id_M;
	invPMat_array[1]=-D_M_inv_G;
	invPMat_array[2]=NULL;
	invPMat_array[3]=NULL;
	invPMat_array[4]=Id_C;
	invPMat_array[5]=NULL;
	invPMat_array[6]=NULL;
	invPMat_array[7]=NULL;
	invPMat_array[8]=Id_R;

	MatCreateNest(PETSC_COMM_WORLD,3,NULL,3,NULL,   PMat_array,&   Pmat);
	MatCreateNest(PETSC_COMM_WORLD,3,NULL,3,NULL,invPMat_array,&invPmat);

//##### Application of the transformation A -> A_hat
	// Declaration
	Mat Mat_array[4];
	Mat A_hat, C_hat, G_hat;
	Mat diag_2M;//Will store 2*diagonal part of M (to approximate the Schur complement)
	
	//Creation of 2*diag(M)
	MatGetDiagonal(M,v_M);
	MatDuplicate(M, MAT_DO_NOT_COPY_VALUES, &diag_2M);
	MatEliminateZeros(diag_2M, PETSC_TRUE);
	MatDiagonalSet(diag_2M, v_M, INSERT_VALUES);
	MatScale(diag_2M,2);//store 2*diagonal part of M
	//Create the matrix 2*diag(M). Why not use MatCreateDiagonal ??? Problem of conversion from MATCONSTANTDIAGONAL to MATAIJ
	//MatCreateConstantDiagonal(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n_u, n_u, 2, &diag_2M);
	//MatDiagonalScale(diag_2M, v_M, NULL);//store 2*diagonal part of M
	/*  Problem of conversion from MATDIAGONAL to MATAIJ
	MatCreateDiagonal(v_M,&diag_2M);
	MatScale(diag_2M,2);//store 2*diagonal part of M
	PetscPrintf(PETSC_COMM_WORLD,"Printing matrix diag_2M before conversion \n");
	MatView( diag_2M, PETSC_VIEWER_STDOUT_WORLD);
	MatConvert(diag_2M,  MATAIJ, MAT_INITIAL_MATRIX, &diag_2Maij);
	PetscPrintf(PETSC_COMM_WORLD,"Printing matrix diag_2M after conversion \n");
	MatView( diag_2Maij, PETSC_VIEWER_STDOUT_WORLD);
	*/

	// Creation of C_hat
	MatMatMult(D,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C_hat);//C_hat contains D*D_M_inv*G
	MatAXPY(C_hat,1.0,C,SUBSET_NONZERO_PATTERN);//C_hat contains C + D*D_M_inv*G

	// Creation of G_hat
	MatMatMult(M,D_M_inv_G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&G_hat);//G_hat contains M*D_M_inv*G
	MatAYPX(G_hat,-1.0,G,UNKNOWN_NONZERO_PATTERN);//G_hat contains G - M*D_M_inv*G

	// Creation of -D
	MatScale(D,-1.0);

	//Creation of global matrices using MatCreateNest
	if( size==1 )
	{
	Mat_array[3]=C_hat;//Top left block of A_hat
	Mat_array[2]=D;//Top right block of A_hat
	Mat_array[1]=G_hat;//Bottom left block of A_hat
	Mat_array[0]=M;//Bottom left block of A_hat
	// Creation of A_hat = transformed+ reordered A_input
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,Mat_array,&A_hat);

	// Creation of Pmat
	Mat_array[0]=diag_2M;
	Mat_array[1]=NULL;//Cancel top right block
	MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,Mat_array,&Pmat);
	}
	else// bug ordonancement matnest en parallèle ?
	{
	Mat_array[0]=C_hat;//Top left block of A_hat
	Mat_array[1]=D;//Top right block of A_hat
	Mat_array[2]=G_hat;//Bottom left block of A_hat
	Mat_array[3]=M;//Bottom left block of A_hat
	//Creation IS pour que la creation du matnest se passe bien
	IS IS_array[2];
	IS_array[0]=is_P;
	IS_array[1]=is_U;	
	MatCreateNest(PETSC_COMM_WORLD,2,IS_array,2,IS_array,Mat_array,&A_hat);
	// Creation of Pmat
	Mat_array[3]=diag_2M;
	Mat_array[2]=NULL;//Cancel top right block
	MatCreateNest(PETSC_COMM_WORLD,2,IS_array,2,IS_array,Mat_array,&Pmat);
	}

//##### Call KSP solver and monitor convergence
	double residu, abstol, rtol=1e-7, dtol;
	int iter, iter1, iter2, numberMaxOfIter;
	int nblocks=2;
	KSP ksp;
	KSP * kspArray;
	KSPType type, type1, type2;
	PC pc, pc1, pc2;
	PCType pctype, pctype1, pctype2;
	
	PetscPrintf(PETSC_COMM_WORLD,"Definition of the KSP solver to test the preconditioner...\n");
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetType(ksp,KSPFGMRES);
	KSPSetOperators(ksp,A_hat,Pmat);
	KSPSetTolerances(ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT, PETSC_DEFAULT);
	KSPGetPC(ksp,&pc);

	PetscPrintf(PETSC_COMM_WORLD,"Setting the preconditioner...\n");
	PCSetType(pc,PCFIELDSPLIT);
	PCFieldSplitSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
	PCFieldSplitSetIS(pc, "0",is_U);
	PCFieldSplitSetIS(pc, "1",is_P);
	PCFieldSplitGetSubKSP( pc, &nblocks, &kspArray);
	KSPSetType( kspArray[0], KSPGMRES);
	KSPSetType( kspArray[1], KSPGMRES);
	KSPGetPC(kspArray[0], &pc1);
	KSPGetPC(kspArray[1], &pc2);
	PCSetType( pc1, PCBJACOBI);
	PCSetType( pc2, PCGAMG);
	PCGAMGSetType( pc2, PCGAMGAGG);

/*		
		KSP * subKSP;
		PC subpc;
		int nlocal;//nb local blocs (should equal 1)
		
		//PCSetUp(pc2);
		KSPSetUp(kspArray[1]);//to set the block Jacobi data structures (including creation of an internal KSP context for each block)
		PCBJacobiGetSubKSP( pc2,&nlocal,NULL,&subKSP);
		PetscPrintf(PETSC_COMM_SELF,"Number of local jacobi blocks : %d\n", nlocal);

		KSPSetType(subKSP[0], KSPPREONLY);//local block solver is same as global
		KSPGetPC(subKSP[0],&subpc);
		PCSetType(subpc,PCLU);

		//PetscOptionsSetValue(NULL,"-fieldsplit_0_ksp_type","gmres");	
		//PetscOptionsSetValue(NULL,"-fieldsplit_0_pc_type","gamg");
		//PetscOptionsSetValue(NULL,"-fieldsplit_0_pc_gamg_type","agg");
		//PetscOptionsSetValue(NULL,"-sub_pc_type ","lu");
		//PetscOptionsSetValue(NULL,"-sub_ksp_type ","preonly");	
*/		

	PCSetFromOptions(pc);
	PCSetUp(pc);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	PetscPrintf(PETSC_COMM_WORLD,"Solving the linear system...\n");
	KSPSolve(ksp,b_input,X_hat);

	//Extract informations about the convergence
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);
	KSPGetIterationNumber(ksp,&iter);
	
	if (reason>0)
		PetscPrintf(PETSC_COMM_WORLD, "\nLinear system converged in %d iterations \n", iter);
	else
		PetscPrintf(PETSC_COMM_WORLD, "\n!!!!!!!!!!!!!!!!!! Linear system diverged  after %d iterations !!!!!!!!!!!!!!\n", iter);
		
	KSPGetResidualNorm( ksp, &residu);
	KSPGetTolerances( ksp, &rtol, &abstol, &dtol, &numberMaxOfIter);
	PCFieldSplitGetSubKSP( pc, &nblocks, &kspArray);
	KSPGetType( ksp, &type);
	KSPGetType( kspArray[0], &type1);
	KSPGetType( kspArray[1], &type2);
	KSPGetIterationNumber(kspArray[0],&iter1);
	KSPGetIterationNumber(kspArray[1],&iter2);
	KSPGetPC(kspArray[0],&pc1);
	KSPGetPC(kspArray[1],&pc2);
	PCGetType( pc, &pctype);
	PCGetType( pc1, &pctype1);
	PCGetType( pc2, &pctype2);
	PetscFree(kspArray);

	PetscPrintf(PETSC_COMM_WORLD, "\n############ : monitoring of the linear solver \n");
	PetscPrintf(PETSC_COMM_WORLD, "Linear solver name: %s, preconditioner %s, %d iterations \n", type, pctype, iter);
	PetscPrintf(PETSC_COMM_WORLD, "    sub solver 1 name : %s, preconditioner %s, %d iterations \n", type1, pctype1, iter1);
	PetscPrintf(PETSC_COMM_WORLD, "    sub solver 2 name: %s, preconditioner %s, %d iterations \n", type2, pctype2, iter2);

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
	Vec X_output, X_output_array[2];
	
	VecGetSubVector( X_hat, is_P, &X_hat_p);
	VecGetSubVector( X_hat, is_U, &X_hat_u);

	VecDuplicate(X_hat_u,&X_u);
	VecDuplicate(X_hat_p,&X_p);
	VecCopy(X_hat_p,X_p);
	MatMult( G, X_hat_p, X_u);
	VecPointwiseMult(X_u,X_u,v);
	VecAYPX( X_u, -1, X_hat_u);

	X_output_array[0] = X_u;
	X_output_array[1] = X_p;

	//VecCreateNest( PETSC_COMM_WORLD, 2, NULL, X_array, &X_output);//This generate an error message : "Nest vector argument 3 not setup "
	VecConcatenate(2, X_output_array, &X_output, NULL);

//##### Compute the error and check it is small
	Vec X_anal_p, X_anal_u;//Pressure and velocity components of the analytic solution
	double error, error_p, error_u;
	
	VecGetSubVector( X_anal, is_P, &X_anal_p);
	VecGetSubVector( X_anal, is_U, &X_anal_u);

	VecAXPY(  X_p, -1, X_anal_p);
	VecNorm(  X_p, NORM_2, &error_p);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error on pressure p : ||X_anal_p - X_num_p|| = %e\n", error_p);
	VecAXPY(  X_u, -1, X_anal_u);
	VecNorm(  X_u, NORM_2, &error_u);
	PetscPrintf(PETSC_COMM_WORLD,"L2 Error on velocity u : ||X_anal_u - X_num_u|| = %e \n", error_u);

	VecAXPY(X_output, -1, X_anal);
	VecNorm( X_output, NORM_2, &error);
	PetscPrintf(PETSC_COMM_WORLD,"L2 total Error : ||X_anal - X_num|| = %e, (remember ||X_anal||=1)\n", error);

	PetscCheck( error < 1.e-5, PETSC_COMM_WORLD, ierr, "Linear system did not return accurate solution. Error is too high\n");
	
//##### Cleaning of the code
	MatDestroy(&M);
	MatDestroy(&G_hat);
	MatDestroy(&C_hat);
	MatDestroy(&D);	
	MatDestroy(&G);
	MatDestroy(&C);
	MatDestroy(&diag_2M);
	MatDestroy(&A_hat);
	MatDestroy(&Pmat);
	MatDestroy(&Id_M);
	MatDestroy(&Id_C);
	MatDestroy(&Id_R);
	
	VecDestroy(&b_input);
	VecDestroy(&X_hat);
	VecDestroy(&X_anal);
	VecDestroy(&v);
	VecDestroy(&v_M);
	VecDestroy(&v_C);
	VecDestroy(&v_R);

	ISDestroy(&is_U);
	ISDestroy(&is_P);

	KSPDestroy(&ksp);
	VecScatterDestroy(&scat);

	PetscFinalize();
	return ierr;
}
