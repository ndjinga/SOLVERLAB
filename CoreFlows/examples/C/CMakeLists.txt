project(testC)

INCLUDE_DIRECTORIES(
  ${CoreFlows_INCLUDES}											    #
  ${PETSC_INCLUDES}
  ${SLEPC_INCLUDES}
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
function(CreateTestExecAndInstall SourceTestFile libList)
     get_filename_component( FILE_BASENAME ${SourceTestFile} NAME_WE) # <path>/testxxx.c --> testxxx
     set( EXECNAME "${FILE_BASENAME}.exe" )                     # testxxx          --> testxxx.exe
     add_executable(${EXECNAME} ${SourceTestFile})                    # compilation of the testxxx.exe 
     set_target_properties(${EXECNAME} PROPERTIES COMPILE_FLAGS "")
     target_link_libraries(${EXECNAME} ${libList})              # provide required lib for testxxx.exe 
     add_test(NAME ${EXECNAME} COMMAND "./${EXECNAME}")     # adding a ctest Test
     install(TARGETS ${EXECNAME} DESTINATION share/examples)
endfunction(CreateTestExecAndInstall)


set( libs_for_tests CoreFlowsLibs )

# copy tests resources (med files etc.) into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

CreateTestExecAndInstall(CoupledTransportDiffusionEquations_1DHeatedChannel.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(TransportEquation_1DHeatedChannel.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod_FE.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_UnstructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredSquares.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DEF_StructuredTetrahedra.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DFV_StructuredTetrahedra.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(testEOS.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(SinglePhase_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DPorosityJump.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem_Implicit.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem_Implicit_LineSearch.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DHeatDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DHeatDrivenCavity_unstructured.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DHeatedChannelInclined.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DLidDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DLidDrivenCavity_unstructured.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DSphericalExplosion_unstructured.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_3DSphericalExplosion_unstructured.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DWallHeatedChannel_ChangeSect.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_2DWallHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_3DHeatDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_HeatedWire_2Branches.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(DriftModel_1DBoilingAssembly.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DChannelGravity.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DPorosityJump.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DPressureLoss.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_1DVidangeReservoir.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_2DInclinedChannelGravity.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_2DInclinedChannelGravityBarriers.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DriftModel_3DCanalCloison.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(FiveEqsTwoFluid_1DBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(FiveEqsTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(FiveEqsTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(FiveEqsTwoFluid_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(FiveEqsTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" )

CreateTestExecAndInstall(IsothermalTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(IsothermalTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" )
#CreateTestExecAndInstall(IsothermalTwoFluid_1DSedimentation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(IsothermalTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(IsothermalTwoFluid_2DVidangeReservoir.cxx  "${libs_for_tests}" )

 
add_executable(WaveSystem_FV_SphericalExplosion_CDMATH.exe WaveSystem_FV_SphericalExplosion_CDMATH.cxx)
target_link_libraries(WaveSystem_FV_SphericalExplosion_CDMATH.exe CoreFlowsLibs )
install(TARGETS WaveSystem_FV_SphericalExplosion_CDMATH.exe DESTINATION share/examples)

add_test(NAME WaveSystem_2DFV_SphericalExplosion_CDMATH_SQUARE         COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe")     
add_test(NAME WaveSystem_2DFV_SphericalExplosion_CDMATH_HEXAGON        COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshHexagonWithTriangles10.med)
add_test(NAME WaveSystem_3DFV_SphericalExplosion_CDMATH_CUBE           COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshCube.med)     
add_test(NAME WaveSystem_3DFV_SphericalExplosion_CDMATH_TETRAHEDRON    COMMAND "./WaveSystem_FV_SphericalExplosion_CDMATH.exe" resources/meshTetrahedron10.med)     

if( SOLVERLAB_WITH_MPI )#Tests parallèles
  add_executable(MEDCouplingSendRecvFieldSameMesh_MPI.exe MEDCouplingSendRecvFieldSameMesh_MPI.cxx)
  target_link_libraries(MEDCouplingSendRecvFieldSameMesh_MPI.exe ${MPI_LIBRARY} medcoupling medloader paramedmem)
  install(TARGETS MEDCouplingSendRecvFieldSameMesh_MPI.exe DESTINATION share/examples)

  add_test(NAME MEDCouplingSendRecvFieldSameMesh_MPI_4Procs.exe                      COMMAND "${MPIEXEC}" "-n" "4" "./MEDCouplingSendRecvFieldSameMesh_MPI.exe" ) 

  add_executable(MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe MEDCouplingSendRecvFieldDifferentMeshes_MPI.cxx)
  target_link_libraries(MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe ${MPI_LIBRARY} medcoupling medloader paramedmem)
  install(TARGETS MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe DESTINATION share/examples)

  add_test(NAME MEDCouplingSendRecvFieldDifferentMeshes_MPI_4Procs.exe                      COMMAND "${MPIEXEC}" "-n" "4" "./MEDCouplingSendRecvFieldDifferentMeshes_MPI.exe" ) 

  add_executable(DiffusionEquation_1DHeatedRod_FE_MPI.exe DiffusionEquation_1DHeatedRod_FE_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(DiffusionEquation_1DHeatedRod_FE_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS DiffusionEquation_1DHeatedRod_FE_MPI.exe DESTINATION share/examples)

  add_test(NAME DiffusionEquation_1DHeatedRod_FE_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./DiffusionEquation_1DHeatedRod_FE_MPI.exe")     

  add_executable(DiffusionEquation_1DHeatedRod_FV_MPI.exe DiffusionEquation_1DHeatedRod_FV_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(DiffusionEquation_1DHeatedRod_FV_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS DiffusionEquation_1DHeatedRod_FV_MPI.exe DESTINATION share/examples)

  add_test(NAME DiffusionEquation_1DHeatedRod_FV_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./DiffusionEquation_1DHeatedRod_FV_MPI.exe")     

  add_executable(TransportEquation_1DHeatedChannel_MPI.exe TransportEquation_1DHeatedChannel_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(TransportEquation_1DHeatedChannel_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS TransportEquation_1DHeatedChannel_MPI.exe DESTINATION share/examples)

  add_test(NAME TransportEquation_1DHeatedChannel_MPI_2Procs.exe         COMMAND "${MPIEXEC}" "-n" "2" "./TransportEquation_1DHeatedChannel_MPI.exe")     

  add_executable(WaveSystem_FV_SphericalExplosion_MPI.exe WaveSystem_FV_SphericalExplosion_MPI.cxx)                    # compilation of the testxxx.exe 
  target_link_libraries(WaveSystem_FV_SphericalExplosion_MPI.exe CoreFlowsLibs ${MPI_LIBRARY})              # provide required lib for testxxx.exe 
  install(TARGETS WaveSystem_FV_SphericalExplosion_MPI.exe DESTINATION share/examples)

  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_SQUARE.exe         COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe")     
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_1Proc_HEXAGON.exe        COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshHexagonWithTriangles10.med) 
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_CUBE.exe           COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshCube.med) 
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_1Proc_TETRAHEDRON.exe    COMMAND "${MPIEXEC}" "-n" "1" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshTetrahedron10.med) 
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_SQUARE.exe         COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe")
  add_test(NAME WaveSystem_2DFV_SphericalExplosion_MPI_2Procs_HEXAGON.exe        COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshHexagonWithTriangles10.med)
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_CUBE.exe           COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshCube.med) 
  add_test(NAME WaveSystem_3DFV_SphericalExplosion_MPI_2Procs_TETRAHEDRON.exe    COMMAND "${MPIEXEC}" "-n" "2" "./WaveSystem_FV_SphericalExplosion_MPI.exe" resources/meshTetrahedron10.med) 

endif( SOLVERLAB_WITH_MPI )

