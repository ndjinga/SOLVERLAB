project(SaddlePointPreconditioner)
cmake_minimum_required(VERSION 3.1)

#Case where the tests are run out of SOLVERLAB
if( NOT DEFINED PETSC_FOUND )
  message ( STATUS "Trying to find PETSc at PETSC_DIR = ${PETSC_DIR}" )
  list (APPEND CMAKE_MODULE_PATH "${CMAKE_FILES}")
  find_package (PETSc 3.4 REQUIRED)
  message ( STATUS "PETSc found. Version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" )
  INCLUDE_DIRECTORIES( ${PETSC_INCLUDES} )
  enable_testing ()
endif( NOT DEFINED PETSC_FOUND )

# copy tests matrices into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Matrices DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target (saddle    COMMAND ctest -R Saddle)# may be replace ctest -R with ctest -L

###### Initial version of the saddle point preconditoner (no Shat)
add_executable(SaddlePointLinearSolver.exe ./SaddlePointLinearSolver.c)
target_link_libraries(SaddlePointLinearSolver.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver.exe DESTINATION share/examples)

    SET(MATRIX_FILE  Matrices/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu -on_error_attach_debugger  )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  Matrices/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu -log_view)     
set_property (TEST SaddlePointLinearSolver_2DCartesian_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DTetra_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DHexa_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DKershaw_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Assembly_PolyMAC_1.petsc  )
    SET(nu  20344  )
    SET(np  9552  )
add_test(NAME SaddlePointLinearSolver_3DAssembly_schur.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DAssembly_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

###### Second version of the saddle point preconditoner (with Shat)
add_executable(SaddlePointLinearSolver_Shat.exe ./SaddlePointLinearSolver_Shat.c)
target_link_libraries(SaddlePointLinearSolver_Shat.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver_Shat.exe DESTINATION share/examples)

    SET(MATRIX_FILE  Matrices/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu -on_error_attach_debugger)     
set_property (TEST SaddlePointLinearSolver_2DKershaw_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  Matrices/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu -log_view)     
set_property (TEST SaddlePointLinearSolver_2DCartesian_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DTetra_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DHexa_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DKershaw_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Assembly_PolyMAC_1.petsc  )
    SET(nu  20344  )
    SET(np  9552  )
add_test(NAME SaddlePointLinearSolver_3DAssembly_Shat_schur.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DAssembly_Shat_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     


###### Parallel (third) version of the saddle point preconditoner (use schur option since preconditioner is lower triangular)
add_executable(SaddlePointLinearSolver_par.exe ./SaddlePointLinearSolver_par.c)
target_link_libraries(SaddlePointLinearSolver_par.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver_par.exe DESTINATION share/examples)

    SET(MATRIX_FILE  Matrices/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu -on_error_attach_debugger)     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    if( SOLVERLAB_WITH_MPI )#Tests parallèles
add_test(NAME SaddlePointLinearSolver_2DKershaw_par_schur_MPI.exe     COMMAND "${MPIEXEC}" "-n" "2" "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -on_error_attach_debugger)     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur_MPI.exe PROPERTY LABELS "C++ PETSC multigrid")     
    endif( SOLVERLAB_WITH_MPI )

    SET(MATRIX_FILE  Matrices/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_1_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Assembly_PolyMAC_1.petsc  )
    SET(nu  20344  )
    SET(np  9552  )
add_test(NAME SaddlePointLinearSolver_3DAssembly_par_schur.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DAssembly_par_schur.exe PROPERTY LABELS "C++ PETSC multigrid")     


###### Fourth version of the saddle point preconditoner (use pc_composite_multiplicative for a triangular preconditioner)
    SET(MATRIX_FILE  Matrices/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    if( SOLVERLAB_WITH_MPI )#Tests parallèles
add_test(NAME SaddlePointLinearSolver_2DKershaw_par_pc_composite_MPI.exe     COMMAND "${MPIEXEC}" "-n" "2" "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np}  -on_error_attach_debugger)   
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite_MPI.exe PROPERTY LABELS "C++ PETSC multigrid")     
    endif( SOLVERLAB_WITH_MPI )

    SET(MATRIX_FILE  Matrices/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  Matrices/Assembly_PolyMAC_1.petsc  )
    SET(nu  20344  )
    SET(np  9552  )
add_test(NAME SaddlePointLinearSolver_3DAssembly_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type multiplicative -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_3DAssembly_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     


