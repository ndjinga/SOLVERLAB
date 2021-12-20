# Graphic interface---------------------------------------------------------------------------------------------#
														#
if  ( (EXISTS $ENV{CONFIGURATION_ROOT_DIR}) AND (EXISTS $ENV{KERNEL_ROOT_DIR}) AND (EXISTS $ENV{GUI_ROOT_DIR}) )

	# Common CMake macros
	# ===================
	SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake configuration files")
	IF(EXISTS ${CONFIGURATION_ROOT_DIR})
	  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
	  INCLUDE(SalomeMacros NO_POLICY_SCOPE)
	ELSE()
	  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
	ENDIF()
						    								#
	# Find Salome KERNEL 
	# ==================
	SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
	IF( EXISTS ${KERNEL_ROOT_DIR} )
	  FIND_PACKAGE(SalomeKERNEL REQUIRED)
	ELSE( EXISTS ${KERNEL_ROOT_DIR} )
	  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR!")
	ENDIF( EXISTS ${KERNEL_ROOT_DIR} )
	    
	MESSAGE(STATUS "SalomeKERNEL_VERSION is ${SalomeKERNEL_VERSION}")

	SET(_all_SCRIPTS
	    SinglePhase/SOLVERLAB_2DSphericalExplosion_HEXAGON.py
	    )
	
	SALOME_INSTALL_SCRIPTS("${_all_SCRIPTS}" ${SALOME_INSTALL_SCRIPT_PYTHON}/tests)
	
	SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PROJECT_BINARY_DIR}/src) 
	SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
	
	ADD_TEST(PythonTests ${PYTHON_EXECUTABLE} ${SALOME_INSTALL_SCRIPT_PYTHON}/tests/SOLVERLAB_2DSphericalExplosion_HEXAGON.py ${CMAKE_INSTALL_PREFIX}/share/meshes/2DHexagonWithTriangles/meshHexagonWithTriangles200.med)
	SET_TESTS_PROPERTIES(PythonTests PROPERTIES ENVIRONMENT "${tests_env}")
						    								#
endif( (EXISTS $ENV{CONFIGURATION_ROOT_DIR}) AND (EXISTS $ENV{KERNEL_ROOT_DIR}) AND (EXISTS $ENV{GUI_ROOT_DIR}) )                                                                                                        #
#---------------------------------------------------------------------------------------------------------------#


#SET(TESTS_LIBRARY_PATH ${CMAKE_BINARY_DIR}/Models/src:${CDMATH_DIR}/lib:${MEDCOUPLING_LIBRARIES}:${MEDFILE_C_LIBRARIES}:${PETSC_DIR}/${PETSC_ARCH}/lib:${PV_LIB_DIR} )
#SET(TESTS_PYTHON_PATH  ${CMAKE_BINARY_DIR}/Models/src:${CMAKE_BINARY_DIR}/swig:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/cdmath:${CDMATH_DIR}/bin/cdmath:${CDMATH_DIR}/bin/cdmath/postprocessing:${MEDCOUPLING_LIBRARIES}:${MEDFILE_C_LIBRARIES}:${PETSC_DIR}/${PETSC_ARCH}/lib:${PV_PYTHON_DIR} )

############# python example tests

function(CreatePythonTest pythonFile)
    add_test(${pythonFile} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${pythonFile})
#    SET_PROPERTY(TEST ${pythonFile}        PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH} )
#    SET_PROPERTY(TEST ${pythonFile} APPEND PROPERTY ENVIRONMENT      PYTHONPATH=${TESTS_PYTHON_PATH}  )
endfunction(CreatePythonTest)

# copy all *.py tests files before build
file(GLOB pythonTestFiles "*.py" )
file(COPY ${pythonTestFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
)

CreatePythonTest(DiffusionEquation/DiffusionEquation_1DHeatedRod.py)
CreatePythonTest(DiffusionEquation/DiffusionEquation_1DHeatedRod_FE.py)
add_test(DiffusionEquation/DiffusionEquation_2DSpherical_FE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionEquation/DiffusionEquation_2DSpherical.py 1)
add_test(DiffusionEquation/DiffusionEquation_2DSpherical_FV ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionEquation/DiffusionEquation_2DSpherical.py 0)

CreatePythonTest(DriftModel/DriftModel_1DBoilingAssembly.py)
CreatePythonTest(DriftModel/DriftModel_1DBoilingChannel.py)
CreatePythonTest(DriftModel/DriftModel_1DChannelGravity.py)
CreatePythonTest(DriftModel/DriftModel_1DDepressurisation.py)
CreatePythonTest(DriftModel/DriftModel_1DPorosityJump.py)
CreatePythonTest(DriftModel/DriftModel_1DPressureLoss.py)
CreatePythonTest(DriftModel/DriftModel_1DRiemannProblem.py)
CreatePythonTest(DriftModel/DriftModel_1DVidangeReservoir.py)
CreatePythonTest(DriftModel/DriftModel_2BranchesBoilingChannels.py)
CreatePythonTest(DriftModel/DriftModel_2DInclinedChannelGravity.py)
CreatePythonTest(DriftModel/DriftModel_2DInclinedChannelGravityTriangles.py)
CreatePythonTest(DriftModel/DriftModel_2DInclinedChannelGravityBarriers.py)
CreatePythonTest(DriftModel/DriftModel_2DBoilingChannelBarrier.py)
CreatePythonTest(DriftModel/DriftModel_2DInclinedBoilingChannelBarrier.py)
CreatePythonTest(DriftModel/DriftModel_2DInclinedBoilingChannel.py)
CreatePythonTest(DriftModel/DriftModel_2DPorosityJump.py)
CreatePythonTest(DriftModel/DriftModel_2DPressureLoss.py)
CreatePythonTest(DriftModel/DriftModel_2DVidangeReservoir.py)
CreatePythonTest(DriftModel/DriftModel_2DVidangeReservoirUnstructured.py)
CreatePythonTest(DriftModel/DriftModel_3DBoilingChannelBarrier.py)

CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_1DBoilingAssembly.py)
CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_1DBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_2DInclinedBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_2DInclinedSedimentation.py)
CreatePythonTest(FiveEqsTwoFluid/FiveEqsTwoFluid_2DVidangeReservoir.py)

CreatePythonTest(IsothermalTwoFluid/IsothermalTwoFluid_1DSedimentation.py)
CreatePythonTest(IsothermalTwoFluid/IsothermalTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(IsothermalTwoFluid/IsothermalTwoFluid_2DVidangeReservoir.py)

CreatePythonTest(SinglePhase/SinglePhase_1DDepressurisation.py)
CreatePythonTest(SinglePhase/SinglePhase_1DHeatedAssembly.py)
CreatePythonTest(SinglePhase/SinglePhase_1DHeatedChannel.py)
CreatePythonTest(SinglePhase/SinglePhase_1DHeatedChannel_Implicit.py)
CreatePythonTest(SinglePhase/SinglePhase_1DRiemannProblem.py)
CreatePythonTest(SinglePhase/SinglePhase_1DRiemannProblem_Implicit.py)
CreatePythonTest(SinglePhase/SinglePhase_1DRiemannProblem_Implicit_LineSearch.py)
CreatePythonTest(SinglePhase/SinglePhase_1DWaterHammer.py)
CreatePythonTest(SinglePhase/SinglePhase_2BranchesHeatedChannels.py)
CreatePythonTest(SinglePhase/SinglePhase_2DVidangeReservoir.py)
CreatePythonTest(SinglePhase/SinglePhase_2DLidDrivenCavity.py)
CreatePythonTest(SinglePhase/SinglePhase_2DLidDrivenCavity_unstructured.py)
CreatePythonTest(SinglePhase/SinglePhase_2DPoiseuilleFlow.py)
CreatePythonTest(SinglePhase/SinglePhase_2DPoiseuilleFlow_restart.py)
CreatePythonTest(SinglePhase/SinglePhase_2DPoiseuilleFlow_outputFields.py)
CreatePythonTest(SinglePhase/SinglePhase_2DAirfoil.py)
CreatePythonTest(SinglePhase/SinglePhase_2DSphericalExplosion_unstructured.py)
CreatePythonTest(SinglePhase/SinglePhase_2DSphericalExplosion_HEXAGON.py)
CreatePythonTest(SinglePhase/SinglePhase_3DSphericalExplosion_TETRAHEDRON.py)
CreatePythonTest(SinglePhase/SinglePhase_2DHeatedChannelInclined.py)
CreatePythonTest(SinglePhase/SinglePhase_2DThermalDiffusion.py)
CreatePythonTest(SinglePhase/SinglePhase_2DWallHeatedChannel_ChangeSect.py)
CreatePythonTest(SinglePhase/SinglePhase_3DSphericalExplosion_unstructured.py)
CreatePythonTest(SinglePhase/SinglePhase_3DHeatDrivenCavity.py)
CreatePythonTest(SinglePhase/SinglePhase_3DVortexTube_NoCone_NoViscosity.py)
CreatePythonTest(SinglePhase/SinglePhase_3DVortexTube_WithCone_NoViscosity.py)

CreatePythonTest(TransportEquation/TransportEquation_1DHeatedChannel.py)
add_test(TransportEquation/TransportEquation_2DSpherical ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TransportEquation/TransportEquation_2DSpherical.py )

CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DEF.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DEF_Neumann.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DFV_StructuredTriangles.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DFV_EquilateralTriangles.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DFV_StructuredSquares.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_2DFV_StructuredSquares_Neumann.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DEF.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DFV_StructuredCubes.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DFV_StructuredTetrahedra.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DEF_RoomCooling.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DVF_RoomCooling_StructuredCubes.py)
CreatePythonTest(StationaryDiffusionEquation/StationaryDiffusionEquation_3DVF_RoomCooling_UnstructuredTetras.py)

if( SOLVERLAB_WITH_MPI )
  add_test(NAME MPI4PY_2Procs                                        COMMAND "${MPIEXEC}" "-n" "2" "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/MPI4PY/testMPI4PY.py")     
  add_test(NAME StationaryDiffusion_TwoIndependentSimulations_2Procs COMMAND "${MPIEXEC}" "-n" "2" "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/MPI4PY/testTwoSimulations.py")     
  add_test(NAME SendRecvMEDCouplingFieldSameMesh_2Procs              COMMAND "${MPIEXEC}" "-n" "2" "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/MPI4PY/testSendRecvFieldSameMesh.py")     
  add_test(NAME SendRecvMEDCouplingFieldDifferentMeshes_2Procs       COMMAND "${MPIEXEC}" "-n" "2" "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/MPI4PY/testSendRecvFieldDifferentMeshes.py")     
  add_test(NAME SendRecvMEDCouplingFieldSubComm_4Procs               COMMAND "${MPIEXEC}" "-n" "4" "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/MPI4PY/testSendRecvFieldSubComm.py")     
endif( SOLVERLAB_WITH_MPI )

############# python convergence tests

function(CreatePythonConvergenceTest model convergenceTest pythonSolver MESHES_AND_PICTURES_PATH)
    FILE(GLOB MY_MESHES_AND_PICTURES  MESHES_AND_PICTURES_PATH )
    file(COPY convergence_${model}_${convergenceTest}.py ../${pythonSolver} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE )

    add_test(convergence_${model}_${convergenceTest} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/convergence_${model}_${convergenceTest}.py)
    #SET_PROPERTY(TEST convergence_${model}_${convergenceTest}        PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH} )
    #SET_PROPERTY(TEST convergence_${model}_${convergenceTest} APPEND PROPERTY ENVIRONMENT      PYTHONPATH=${TESTS_PYTHON_PATH} )

endfunction(CreatePythonConvergenceTest)

add_subdirectory (Convergence/StationaryDiffusion)
add_subdirectory (Convergence/SinglePhase)

