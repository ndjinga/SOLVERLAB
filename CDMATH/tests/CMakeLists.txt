# Nom du projet tests
PROJECT(tests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (CDMATH_WITH_TESTS)
  INCLUDE_DIRECTORIES(
    ${MEDCOUPLING_INCLUDE_DIR}
    ${BASE_DIR}/inc
    ${MESH_DIR}/inc
    cdmath
    )

  SET(MAIN_SOURCES
    MainTests.cxx
    )

  ADD_EXECUTABLE(TestCDMATH ${MAIN_SOURCES})

  if  (CDMATH_WITH_PETSC)
      INCLUDE_DIRECTORIES(
        ${LINEARSOLVER_DIR}/inc
        ${PETSC_INCLUDES}
        ${SLEPC_INCLUDES}
        )
      SET(_extra_lib_TESTS linearsolver ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
      SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS -DCDMATH_WITH_PETSC)
      SET(TESTS_LIBRARY_PATH ${CMAKE_BINARY_DIR}/base/lib:${CMAKE_BINARY_DIR}/mesh/lib:${CMAKE_BINARY_DIR}/linearsolver/lib:${MEDCOUPLING_ROOT_DIR}/lib:${MEDFILE_ROOT_DIR}/lib:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT}/lib ) 
  else (CDMATH_WITH_PETSC)
      SET(_extra_lib_TESTS)
      SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS "")
      SET(TESTS_LIBRARY_PATH ${CMAKE_BINARY_DIR}/base/lib:${CMAKE_BINARY_DIR}/mesh/lib:${MEDCOUPLING_ROOT_DIR}/lib:${MEDFILE_ROOT_DIR}/lib:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT}/lib ) 
  endif(CDMATH_WITH_PETSC)

  if   (CDMATH_WITH_PYTHON)
    SET(TESTS_PYTHON_PATH ${CMAKE_BINARY_DIR}/swig:${CMAKE_BINARY_DIR}/postprocessing ) 
  endif(CDMATH_WITH_PYTHON)

  TARGET_LINK_LIBRARIES(TestCDMATH TestMesh base mesh ${CPPUNIT_LIBS} ${_extra_lib_TESTS})

  ADD_TEST(UnitTestCDMATH TestCDMATH)
  SET_PROPERTY(TEST UnitTestCDMATH PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}")
  INSTALL(TARGETS TestCDMATH DESTINATION tests/unit_tests)

  ADD_SUBDIRECTORY(cdmath)

  IF (CDMATH_WITH_PYTHON)
    ADD_SUBDIRECTORY(swig)
  ENDIF (CDMATH_WITH_PYTHON)

ENDIF (CDMATH_WITH_TESTS)

ADD_SUBDIRECTORY(ressources)

FILE(COPY doc/ DESTINATION convergence_tests)
INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/convergence_tests DESTINATION doc/)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
  # Examples tests
  set (examples_dir ${CDMATH_SOURCE_DIR}/tests/examples)
  add_subdirectory (${examples_dir})
  add_subdirectory (validation)

  add_custom_target (tests_CDMATH COMMAND ctest -O tests.log)
  add_custom_target (check COMMAND ctest -E 'validation|Example')# may be replace ctest -E with ctest -LE
  add_custom_target (examples COMMAND ctest -R Example)# may be replace ctest -R with ctest -L
  add_custom_target (validation COMMAND ctest -R validation)# may be replace ctest -R with ctest -L

  add_custom_target (spectrum COMMAND ctest -R Spectrum)# may be replace ctest -R with ctest -L
  add_custom_target (beltrami COMMAND ctest -R Beltrami)# may be replace ctest -R with ctest -L
  add_custom_target (poisson COMMAND ctest -R Poisson)# may be replace ctest -R with ctest -L
  add_custom_target (diffusion COMMAND ctest -R Diffusion)# may be replace ctest -R with ctest -L
  add_custom_target (transport COMMAND ctest -R Transport)# may be replace ctest -R with ctest -L
  add_custom_target (wave COMMAND ctest -R WaveSystem)# may be replace ctest -R with ctest -L
  add_custom_target (upwind COMMAND ctest -R Upwind)# may be replace ctest -R with ctest -L
  add_custom_target (pstag COMMAND ctest -R PStag)# may be replace ctest -R with ctest -L
  add_custom_target (staggered COMMAND ctest -R Staggered)# may be replace ctest -R with ctest -L
  add_custom_target (centered COMMAND ctest -R Centered)# may be replace ctest -R with ctest -L
  add_custom_target (triangles COMMAND ctest -R triangles)# may be replace ctest -R with ctest -L
  add_custom_target (checkerboard COMMAND ctest -R checkerboard)# may be replace ctest -R with ctest -L
  add_custom_target (squares COMMAND ctest -R squares)# may be replace ctest -R with ctest -L
  add_custom_target (scaling COMMAND ctest -R '_scaling' -E no_scaling)# may be replace ctest -R with ctest -L
  add_custom_target (noscaling COMMAND ctest -R no_scaling)# may be replace ctest -R with ctest -L
  add_custom_target (2D COMMAND ctest -R 2D)# may be replace ctest -R with ctest -L
  add_custom_target (3D COMMAND ctest -R 3D)# may be replace ctest -R with ctest -L
  add_custom_target (FE COMMAND ctest -R FE)# may be replace ctest -R with ctest -L
  add_custom_target (FV COMMAND ctest -R FV)# may be replace ctest -R with ctest -L
  add_custom_target (shock COMMAND ctest -R Shock)# may be replace ctest -R with ctest -L
  add_custom_target (vortex COMMAND ctest -R Vortex)# may be replace ctest -R with ctest -L
  add_custom_target (riemann COMMAND ctest -R RiemannProblem)# may be replace ctest -R with ctest -L
endif ()

