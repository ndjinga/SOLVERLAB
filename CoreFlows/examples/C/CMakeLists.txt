project(testC)

INCLUDE_DIRECTORIES(
  ${CoreFlows_INCLUDES}											    #
)


if(CMAKE_COMPILER_IS_GNUCXX)
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
    include(CodeCoverage)
    setup_target_for_coverage(cov ctest coverage)
    endif()
endif()


##################################### test generation with ctest

# this function creates a target and a ctest test
function(CreateTestExec SourceTestFile libList)
     get_filename_component( FILE_BASENAME ${SourceTestFile} NAME_WE) # <path>/testxxx.c --> testxxx
     set( EXECNAME "${FILE_BASENAME}.exe" )                     # testxxx          --> testxxx.exe
     add_executable(${EXECNAME} ${SourceTestFile})                    # compilation of the testxxx.exe 
     set_target_properties(${EXECNAME} PROPERTIES COMPILE_FLAGS "")
     target_link_libraries(${EXECNAME} ${libList})              # provide required lib for testxxx.exe 
     add_test(${FILE_BASENAME} ${EXECNAME} "./${EXECNAME}")     # adding a ctest Test
endfunction(CreateTestExec)

# this function creates a target and a ctest test
# and also create install rules for copying the example
# in the install dir
function(CreateTestExecAndInstall SourceTestFile libList testLabels)
     get_filename_component( FILE_BASENAME ${SourceTestFile} NAME_WE) # <path>/testxxx.c --> testxxx
     set( EXECNAME "${FILE_BASENAME}.exe" )                     # testxxx          --> testxxx.exe
     add_executable(${EXECNAME} ${SourceTestFile})                    # compilation of the testxxx.exe 
     set_target_properties(${EXECNAME} PROPERTIES COMPILE_FLAGS "")
     target_link_libraries(${EXECNAME} ${libList})              # provide required lib for testxxx.exe 
     add_test(NAME ${EXECNAME} COMMAND "./${EXECNAME}")     # adding a ctest Test
     set_property (TEST ${EXECNAME} PROPERTY LABELS ${testLabels})     
     install(TARGETS ${EXECNAME} DESTINATION share/examples)
endfunction(CreateTestExecAndInstall)


set( libs_for_tests CoreFlowsLibs ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

# copy tests resources (med files etc.) into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

CreateTestExecAndInstall(CoupledTransportDiffusionEquations_1DHeatedChannel.cxx  "${libs_for_tests}" "C++ fv coupled transport diffusion 1D heating structured coreflows")

CreateTestExecAndInstall(TransportEquation_1DHeatedChannel.cxx  "${libs_for_tests}" "C++ fv transport 1D heating structured coreflows" )

CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod.cxx  "${libs_for_tests}" "C++ fv diffusion 1D explicit heating structured coreflows" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod_Implicit.cxx  "${libs_for_tests}" "C++ fv diffusion 1D implicit heating structured coreflows" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod_FE.cxx  "${libs_for_tests}" "C++ fe diffusion 1D explicit heating structured coreflows" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod_FE_Implicit.cxx  "${libs_for_tests}" "C++ fe diffusion 1D implicit heating structured coreflows" )
add_executable(DiffusionEquation_2D_SQUARE_Implicit.exe DiffusionEquation_2D_SQUARE_Implicit.cxx)
target_link_libraries(DiffusionEquation_2D_SQUARE_Implicit.exe CoreFlowsLibs )
install(TARGETS DiffusionEquation_2D_SQUARE_Implicit.exe DESTINATION share/examples)
#Generate several tests from the same executable ./DiffusionEquation_2D_SQUARE_Implicit.exe using various parameters (finite elemets vs volumes, Dirichlet vs Neumann)
add_test(NAME DiffusionEquation_2D_SQUARE_Implicit_FE_Neumann.exe         COMMAND ./DiffusionEquation_2D_SQUARE_Implicit.exe FE NEUMANN )     
  set_property (TEST DiffusionEquation_2D_SQUARE_Implicit_FE_Neumann.exe PROPERTY LABELS "C++ diffusion heating 2D fe Neumann squares SQUARE coreflows")     
add_test(NAME DiffusionEquation_2D_SQUARE_Implicit_FE_Dirichlet.exe         COMMAND ./DiffusionEquation_2D_SQUARE_Implicit.exe FE Dirichlet )     
  set_property (TEST DiffusionEquation_2D_SQUARE_Implicit_FE_Dirichlet.exe PROPERTY LABELS "C++ diffusion heating 2D fe Dirichlet squares SQUARE coreflows")     
add_test(NAME DiffusionEquation_2D_SQUARE_Implicit_FV_Neumann.exe         COMMAND ./DiffusionEquation_2D_SQUARE_Implicit.exe FV Neumann )     
  set_property (TEST DiffusionEquation_2D_SQUARE_Implicit_FV_Neumann.exe PROPERTY LABELS "C++ diffusion heating 2D fv Neumann squares SQUARE coreflows")     
add_test(NAME DiffusionEquation_2D_SQUARE_Implicit_FV_Dirichlet.exe         COMMAND ./DiffusionEquation_2D_SQUARE_Implicit.exe FV DIRICHLET )     
  set_property (TEST DiffusionEquation_2D_SQUARE_Implicit_FV_Dirichlet.exe PROPERTY LABELS "C++ diffusion heating 2D fv Dirichlet squares SQUARE coreflows")     

CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles.cxx  "${libs_for_tests}" "C++ fe stationary diffusion 2D dirichlet structured triangles coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" "C++ fe stationary diffusion 2D neumann structured triangles coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_UnstructuredTriangles.cxx  "${libs_for_tests}" "C++ fe stationary diffusion 2D dirichlet unstructured triangles coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles.cxx  "${libs_for_tests}" "C++ fv stationary diffusion 2D dirichlet structured triangles coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" "C++ fv stationary diffusion 2D neumann structured triangles coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredSquares.cxx  "${libs_for_tests}" "C++ fv stationary diffusion 2D dirichlet structured squares coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DEF_StructuredTetrahedra.cxx  "${libs_for_tests}" "C++ fe stationary diffusion 3D dirichlet structured tetrahedra coreflows" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DFV_StructuredTetrahedra.cxx  "${libs_for_tests}" "C++ fv stationary diffusion 3D dirichlet structured tetrahedra coreflows" )

CreateTestExecAndInstall(testEOS.cxx  "${libs_for_tests}" "C++ eos coreflows")

CreateTestExecAndInstall(SinglePhase_1DDepressurisation.cxx  "${libs_for_tests}" "C++ fv singlephase 1D depressurisation structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DHeatedChannel.cxx  "${libs_for_tests}" "C++ fv singlephase 1D heating structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DPorosityJump.cxx  "${libs_for_tests}" "C++ fv singlephase 1D porosity structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem.cxx  "${libs_for_tests}" "C++ fv singlephase 1D riemann structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem_1bar.cxx  "${libs_for_tests}" "C++ fv singlephase 1D riemann structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem_Implicit.cxx  "${libs_for_tests}" "C++ fv singlephase 1D riemann implicit structured coreflows" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem_Implicit_LineSearch.cxx  "${libs_for_tests}" "C++ fv singlephase 1D riemann implicit linesearch structured coreflows" )
CreateTestExecAndInstall(SinglePhase_2DHeatDrivenCavity.cxx  "${libs_for_tests}" "C++ fv singlephase 2D heating structured squares coreflows" )
CreateTestExecAndInstall(SinglePhase_2DHeatDrivenCavity_unstructured.cxx  "${libs_for_tests}" "C++ fv singlephase 2D heating unstructured squares coreflows" )
CreateTestExecAndInstall(SinglePhase_2DHeatedChannelInclined.cxx  "${libs_for_tests}" "C++ fv singlephase 2D heating structured squares inclined coreflows" )
CreateTestExecAndInstall(SinglePhase_2DLidDrivenCavity.cxx  "${libs_for_tests}" "C++ fv singlephase 2D viscosity structured squares coreflows" )
CreateTestExecAndInstall(SinglePhase_2DLidDrivenCavity_unstructured.cxx  "${libs_for_tests}" "C++ fv singlephase 2D viscosity unstructured triangles coreflows" )
CreateTestExecAndInstall(SinglePhase_2DSphericalExplosion_unstructured.cxx  "${libs_for_tests}" "C++ fv singlephase 2D riemann unstructured triangles coreflows" )
CreateTestExecAndInstall(SinglePhase_3DSphericalExplosion_unstructured.cxx  "${libs_for_tests}" "C++ fv singlephase 3D riemann unstructured triangles coreflows" )
CreateTestExecAndInstall(SinglePhase_2DWallHeatedChannel_ChangeSect.cxx  "${libs_for_tests}" "C++ fv singlephase 2D heating structured squares coreflows" )
CreateTestExecAndInstall(SinglePhase_2DWallHeatedChannel.cxx  "${libs_for_tests}" "C++ fv singlephase 2D heating structured squares coreflows" )
CreateTestExecAndInstall(SinglePhase_3DHeatDrivenCavity.cxx  "${libs_for_tests}" "C++ fv singlephase 3D heating structured cubes coreflows" )
CreateTestExecAndInstall(SinglePhase_HeatedWire_2Branches.cxx  "${libs_for_tests}" "C++ fv singlephase 1D branches heating unstructured coreflows" )

CreateTestExecAndInstall(IsothermalSinglePhase_1DRiemannProblem.cxx  "${libs_for_tests}" "C++ fv isothermalsinglephase 1D riemann explicit upwindstructured coreflows" )
CreateTestExecAndInstall(IsothermalSinglePhase_1DRiemannProblem_Implicit.cxx  "${libs_for_tests}" "C++ fv isothermalsinglephase 1D implicit upwind riemann structured coreflows" )
CreateTestExecAndInstall(IsothermalSinglePhase_1DChannelGravity_Incompressible.cxx "${libs_for_tests}" "C++ fv isothermalsinglephase 1D incompressible implicit staggered gravity structured coreflows" )
CreateTestExecAndInstall(IsothermalSinglePhase_2DLidDrivenCavity_Incompressible.cxx "${libs_for_tests}" "C++ fv isothermalsinglephase 1D incompressible  implicit staggered viscosity structured coreflows" )
CreateTestExecAndInstall(IsothermalSinglePhase_2DPoiseuilleFlow.cxx "${libs_for_tests}" "C++ fv isothermalsinglephase 1D incompressible  implicit staggered viscosity structured coreflows" )

CreateTestExecAndInstall(DriftModel_1DBoilingAssembly.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D heating boiling structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DBoilingChannel.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D heating boiling structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DChannelGravity.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D gravity structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DDepressurisation.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D depressurisation structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DPorosityJump.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D porosity structured" )
CreateTestExecAndInstall(DriftModel_1DPressureLoss.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D pressureloss structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DRiemannProblem.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D riemann structured coreflows" )
CreateTestExecAndInstall(DriftModel_1DVidangeReservoir.cxx  "${libs_for_tests}" "C++ fv driftmodel 1D gravity drainage structured coreflows" )
CreateTestExecAndInstall(DriftModel_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" "C++ fv driftmodel 2D heating boiling structured inclined coreflows" )
CreateTestExecAndInstall(DriftModel_2DInclinedChannelGravity.cxx  "${libs_for_tests}" "C++ fv driftmodel 2D gravity structured inclined coreflows" )
CreateTestExecAndInstall(DriftModel_2DInclinedChannelGravityBarriers.cxx  "${libs_for_tests}" "C++ fv driftmodel 2D gravity barrier structured inclined coreflows" )
CreateTestExecAndInstall(DriftModel_3DCanalCloison.cxx  "${libs_for_tests}" "C++ fv driftmodel 3D barrier structured coreflows" )

CreateTestExecAndInstall(FiveEqsTwoFluid_1DBoilingChannel.cxx  "${libs_for_tests}" "C++ fv fiveqs 1D heating boiling structured coreflows" )
CreateTestExecAndInstall(FiveEqsTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" "C++ fv fiveqs 1D depressurisation structured coreflows" )
CreateTestExecAndInstall(FiveEqsTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" "C++ fv fiveqs 1D riemann structured coreflows" )
CreateTestExecAndInstall(FiveEqsTwoFluid_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" "C++ fv fiveqs 2D heating boiling structured inclined coreflows" )
CreateTestExecAndInstall(FiveEqsTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" "C++ fv fiveqs 2D gravity sedimentation structured inclined coreflows" )

CreateTestExecAndInstall(IsothermalTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" "C++ fv isothermaltwofluid 1D depressurisation structured coreflows" )
CreateTestExecAndInstall(IsothermalTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" "C++ fv isothermaltwofluid 1D riemann structured coreflows" )
#CreateTestExecAndInstall(IsothermalTwoFluid_1DSedimentation.cxx  "${libs_for_tests}" "C++ fv isothermaltwofluid 1D gravity sedimentation structured coreflows" )
CreateTestExecAndInstall(IsothermalTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" "C++ fv isothermaltwofluid 2D gravity sedimentation structured inclined coreflows" )
CreateTestExecAndInstall(IsothermalTwoFluid_2DVidangeReservoir.cxx  "${libs_for_tests}" "C++ fv isothermaltwofluid 2D gravity drainage structured coreflows" )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/PETSc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

###### Initial version of the saddle point preconditoner (no Shat)
add_executable(SaddlePointLinearSolver.exe PETSc/SaddlePointLinearSolver.c)
target_link_libraries(SaddlePointLinearSolver.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver.exe DESTINATION share/examples)

    SET(MATRIX_FILE  PETSc/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_2DKershaw.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  PETSc/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_2DCartesian.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DTetra.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DHexa.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw.exe         COMMAND "./SaddlePointLinearSolver.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DKershaw.exe PROPERTY LABELS "C++ PETSC multigrid")     

###### Second version of the saddle point preconditoner (with Shat)
add_executable(SaddlePointLinearSolver_Shat.exe PETSc/SaddlePointLinearSolver_Shat.c)
target_link_libraries(SaddlePointLinearSolver_Shat.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver_Shat.exe DESTINATION share/examples)

    SET(MATRIX_FILE  PETSc/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_Shat.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_2DKershaw_Shat.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  PETSc/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_Shat.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_2DCartesian_Shat.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_Shat.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DTetra_Shat.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_Shat.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DHexa_Shat.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_Shat.exe         COMMAND "./SaddlePointLinearSolver_Shat.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_type agg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type lu >> res.txt)     
set_property (TEST SaddlePointLinearSolver_3DKershaw_Shat.exe PROPERTY LABELS "C++ PETSC multigrid")     


###### Parallel (third) version of the saddle point preconditoner (no Shat)
add_executable(SaddlePointLinearSolver_par.exe PETSc/SaddlePointLinearSolver_par.c)
target_link_libraries(SaddlePointLinearSolver_par.exe ${PETSC_LIBRARIES} )
install(TARGETS SaddlePointLinearSolver_par.exe DESTINATION share/examples)

    SET(MATRIX_FILE  PETSc/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_par.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type upper -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  PETSc/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_par.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type upper -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_par.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type upper -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_par.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type upper -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_par.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type upper -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par.exe PROPERTY LABELS "C++ PETSC multigrid")     


###### Fourth version of the saddle point preconditoner (no Schur)
    SET(MATRIX_FILE  PETSc/Kershaw2D_PolyMAC_1.petsc  )
    SET(nu  612  )
    SET(np  289  )
add_test(NAME SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     
 
    SET(MATRIX_FILE  PETSc/Cartesian_PolyMAC_1.petsc  )
    SET(nu  40  )
    SET(np  16  )
add_test(NAME SaddlePointLinearSolver_2DCartesian_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Tetra_PolyMAC_1.petsc  )
    SET(nu  500  )
    SET(np  215  )
add_test(NAME SaddlePointLinearSolver_3DTetra_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Hexa_PolyMAC_1.petsc  )
    SET(nu  36  )
    SET(np  8  )
add_test(NAME SaddlePointLinearSolver_3DHexa_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     

    SET(MATRIX_FILE  PETSc/Kershaw_PolyMAC_1.petsc  )
    SET(nu  1728  )
    SET(np  512  )
add_test(NAME SaddlePointLinearSolver_3DKershaw_par_pc_composite.exe         COMMAND "./SaddlePointLinearSolver_par.exe" -f0 ${MATRIX_FILE} -nU ${nu} -nP ${np} -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type gamg -fieldsplit_1_pc_gamg_type agg -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type lu )     
set_property (TEST SaddlePointLinearSolver_2DKershaw_par_pc_composite.exe PROPERTY LABELS "C++ PETSC multigrid")     


add_executable(WaveSystem_FV_SphericalExplosion_CDMATH.exe WaveSystem_FV_SphericalExplosion_CDMATH.cxx)
target_link_libraries(WaveSystem_FV_SphericalExplosion_CDMATH.exe CoreFlowsLibs )
install(TARGETS WaveSystem_FV_SphericalExplosion_CDMATH.exe DESTINATION share/examples)

#Generate several tests from the same executable WaveSystem_FV_SphericalExplosion_CDMATH.exe (change mesh name)
add_test(NAME WaveSystem_2DFV_SphericalExplosion_CDMATH_SQUARE.exe         COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe")     
set_property (TEST WaveSystem_2DFV_SphericalExplosion_CDMATH_SQUARE.exe PROPERTY LABELS "C++ wave 2D fv riemann squares SQUARE cdmath")     
add_test(NAME WaveSystem_2DFV_SphericalExplosion_CDMATH_HEXAGON.exe        COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshHexagonWithTriangles10.med)
set_property (TEST WaveSystem_2DFV_SphericalExplosion_CDMATH_HEXAGON.exe PROPERTY LABELS "C++ wave 2D fv riemann triangles HEXAGON cdmath")     
add_test(NAME WaveSystem_3DFV_SphericalExplosion_CDMATH_CUBE.exe           COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshCube.med)     
set_property (TEST WaveSystem_3DFV_SphericalExplosion_CDMATH_CUBE.exe PROPERTY LABELS "C++ wave 3D fv riemann tetrahedra CUBE cdmath")     
add_test(NAME WaveSystem_3DFV_SphericalExplosion_CDMATH_TETRAHEDRON.exe    COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshTetrahedron10.med)     
set_property (TEST WaveSystem_3DFV_SphericalExplosion_CDMATH_TETRAHEDRON.exe PROPERTY LABELS "C++ wave 3D fv riemann tetrahedra TETRAHEDRON cdmath")     

if( SOLVERLAB_WITH_MPI )#Tests parallèles
  add_executable(MEDCouplingSendRecvFieldSameMesh_MPI.exe MEDCouplingSendRecvFieldSameMesh_MPI.cxx)
  target_link_libraries(MEDCouplingSendRecvFieldSameMesh_MPI.exe ${MPI_LIBRARY} medcoupling medloader paramedmem)
  install(TARGETS MEDCouplingSendRecvFieldSameMesh_MPI.exe DESTINATION share/examples)

  add_test(NAME MEDCouplingSendRecvFieldSameMesh_MPI_4Procs.exe                      COMMAND "${MPIEXEC}" "-n" "4" "./MEDCouplingSendRecvFieldSameMesh_MPI.exe" ) 
  set_property (TEST MEDCouplingSendRecvFieldSameMesh_MPI_4Procs.exe PROPERTY LABELS "C++ medcoupling mpi 4procs cdmath")     

  add_executable(MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe MEDCouplingSendRecvFieldDifferentMeshes_MPI.cxx)
  target_link_libraries(MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe ${MPI_LIBRARY} medcoupling medloader paramedmem)
  install(TARGETS MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe DESTINATION share/examples)

  add_test(NAME MEDCouplingSendRecvFieldDifferentMeshes_MPI_4Procs.exe                      COMMAND "${MPIEXEC}" "-n" "4" "./MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe" ) 
  set_property (TEST MEDCouplingSendRecvFieldDifferentMeshes_MPI_4Procs.exe PROPERTY LABELS "C++ medcoupling mpi 4procs cdmath")     

  add_executable(DiffusionEquation_1DHeatedRod_FE_MPI.exe DiffusionEquation_1DHeatedRod_FE_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(DiffusionEquation_1DHeatedRod_FE_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS DiffusionEquation_1DHeatedRod_FE_MPI.exe DESTINATION share/examples)

  add_test(NAME DiffusionEquation_1DHeatedRod_FE_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./DiffusionEquation_1DHeatedRod_FE_MPI.exe")     
  set_property (TEST MEDCouplingSendRecvFieldDifferentMeshes_MPI_4Procs.exe PROPERTY LABELS "C++ medcoupling mpi 2procs cdmath")     

  add_executable(DiffusionEquation_1DHeatedRod_FV_MPI.exe DiffusionEquation_1DHeatedRod_FV_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(DiffusionEquation_1DHeatedRod_FV_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS DiffusionEquation_1DHeatedRod_FV_MPI.exe DESTINATION share/examples)

  add_test(NAME DiffusionEquation_1DHeatedRod_FV_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./DiffusionEquation_1DHeatedRod_FV_MPI.exe")     
  set_property (TEST DiffusionEquation_1DHeatedRod_FV_MPI_2Procs.exe PROPERTY LABELS "C++ fv diffusion 1D heating structured mpi 2procs coreflows")     

  add_executable(TransportEquation_1DHeatedChannel_MPI.exe TransportEquation_1DHeatedChannel_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(TransportEquation_1DHeatedChannel_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS TransportEquation_1DHeatedChannel_MPI.exe DESTINATION share/examples)

  add_test(NAME TransportEquation_1DHeatedChannel_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./TransportEquation_1DHeatedChannel_MPI.exe")     
  set_property (TEST TransportEquation_1DHeatedChannel_MPI_2Procs.exe PROPERTY LABELS "C++ fv transport 1D heating structured mpi 2procs coreflows")     

  add_executable(WaveSystem_FV_SphericalExplosion_MPI.exe WaveSystem_FV_SphericalExplosion_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(WaveSystem_FV_SphericalExplosion_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS WaveSystem_FV_SphericalExplosion_MPI.exe DESTINATION share/examples)

  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_SQUARE.exe         COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe")     
  set_property (TEST WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_SQUARE.exe PROPERTY LABELS "C++ wave 2D fv riemann mpi 1Proc squares SQUARE cdmath")     
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_HEXAGON.exe        COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshHexagonWithTriangles10.med) 
  set_property (TEST WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_HEXAGON.exe PROPERTY LABELS "C++ wave 2D fv riemann mpi 1Proc triangles HEXAGON cdmath")     
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_CUBE.exe           COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshCube.med) 
  set_property (TEST WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_CUBE.exe PROPERTY LABELS "C++ wave 3D fv riemann mpi 1Proc tetrahedra CUBE cdmath")     
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_TETRAHEDRON.exe    COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshTetrahedron10.med) 
  set_property (TEST WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_TETRAHEDRON.exe PROPERTY LABELS "C++ wave 3D fv riemann mpi 1Proc tetrahedra TETRAHEDRON cdmath")     
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_SQUARE.exe         COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe")
  set_property (TEST WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_SQUARE.exe PROPERTY LABELS "C++ wave 2D fv riemann mpi 2Proc squares SQUARE cdmath")     
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_HEXAGON.exe        COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshHexagonWithTriangles10.med)
  set_property (TEST WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_HEXAGON.exe PROPERTY LABELS "C++ wave 2D fv riemann mpi 2Proc triangles HEXAGON cdmath")     
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_CUBE.exe           COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshCube.med) 
  set_property (TEST WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_CUBE.exe PROPERTY LABELS "C++ wave 3D fv riemann mpi 2Proc tetrahedra CUBE cdmath")     
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_TETRAHEDRON.exe    COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshTetrahedron10.med) 
  set_property (TEST WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_TETRAHEDRON.exe PROPERTY LABELS "C++ wave 3D fv riemann mpi 2Proc tetrahedra TETRAHEDRON cdmath")

endif( SOLVERLAB_WITH_MPI )

