project(testC)

INCLUDE_DIRECTORIES(
  ${CoreFlows_INCLUDES}											    #
)


SET(_extra_lib_CoreFlows CoreFlows ${PETSC_LIBRARIES} ${CDMATH_LIBRARIES})


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


set( libs_for_tests ${_extra_lib_CoreFlows} )

# copy tests resources (med files etc.) into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

CreateTestExecAndInstall(CoupledTransportDiffusionEquations_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(DiffusionEquation_1DHeatedRod_FE.cxx  "${libs_for_tests}" )
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
CreateTestExecAndInstall(SinglePhase_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DPorosityJump.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(SinglePhase_1DRiemannProblem.cxx  "${libs_for_tests}" )
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
CreateTestExecAndInstall(TransportEquation_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DEF_UnstructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredTriangles_Neumann.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_2DFV_StructuredSquares.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DEF_StructuredTetrahedra.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(StationaryDiffusionEquation_3DFV_StructuredTetrahedra.cxx  "${libs_for_tests}" )
CreateTestExecAndInstall(testEOS.cxx  "${libs_for_tests}" )

