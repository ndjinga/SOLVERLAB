# Nom du projet tests
PROJECT(tests)

IF (CDMATH_WITH_TESTS)
  INCLUDE_DIRECTORIES(
    ${MEDCOUPLING_INCLUDE_DIR}
    ${BASE_DIR}/inc
    ${MESH_DIR}/inc
    cdmath
    ${CPPUNIT_INCLUDE_DIR}
    )

  SET(MAIN_SOURCES
    MainTests.cxx
    )

  ADD_EXECUTABLE(UnitTestsCDMATH ${MAIN_SOURCES})

  if  (CDMATH_WITH_PETSC)
      INCLUDE_DIRECTORIES(
        ${LINEARSOLVER_DIR}/inc
        ${PETSC_INCLUDES}
        ${SLEPC_INCLUDES}
        )
      SET(_extra_lib_TESTS linearsolver ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
      SET_TARGET_PROPERTIES(UnitTestsCDMATH PROPERTIES COMPILE_FLAGS -DCDMATH_WITH_PETSC)
      SET(TESTS_LIBRARY_PATH ${CMAKE_BINARY_DIR}/CDMATH/base/lib:${CMAKE_BINARY_DIR}/CDMATH/mesh/lib:${CMAKE_BINARY_DIR}/CDMATH/linearsolver/lib:${MEDCOUPLING_ROOT_DIR}/lib:${MEDFILE_ROOT_DIR}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT}/lib:${Python_ROOT_DIR}/lib/:${PV_LIB_DIR} ) 
  else (CDMATH_WITH_PETSC)
      SET(_extra_lib_TESTS)
      SET_TARGET_PROPERTIES(UnitTestsCDMATH PROPERTIES COMPILE_FLAGS "")
      SET(TESTS_LIBRARY_PATH ${CMAKE_BINARY_DIR}/CDMATH/base/lib:${CMAKE_BINARY_DIR}/CDMATH/mesh/lib:${MEDCOUPLING_ROOT_DIR}/lib:${MEDFILE_ROOT_DIR}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${HDF5_ROOT}/lib:${Python_ROOT_DIR}/lib/:${PV_LIB_DIR}) 
  endif(CDMATH_WITH_PETSC)

  if   (CDMATH_WITH_PYTHON)
    SET(TESTS_PYTHON_PATH ${CMAKE_BINARY_DIR}/CDMATH/swig:${CMAKE_BINARY_DIR}/CDMATH/postprocessing:${MEDCOUPLING_ROOT_DIR}/lib:${MEDFILE_ROOT_DIR}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PV_PYTHON_DIR} ) 
  endif(CDMATH_WITH_PYTHON)

  TARGET_LINK_LIBRARIES(UnitTestsCDMATH TestMesh base mesh ${CPPUNIT_LIBS} ${_extra_lib_TESTS})

  ADD_TEST(UnitTestsCDMATH.exe UnitTestsCDMATH)
  INSTALL(TARGETS UnitTestsCDMATH DESTINATION bin/unit_tests)

  if   ( ${CPPUNIT_IS_OK} )
    ADD_SUBDIRECTORY(cdmath)
  endif( ${CPPUNIT_IS_OK} )

  IF (CDMATH_WITH_PYTHON)
    ADD_SUBDIRECTORY(swig)
  ENDIF (CDMATH_WITH_PYTHON)

ENDIF (CDMATH_WITH_TESTS)

ADD_SUBDIRECTORY(ressources)
SET(MED_MESHES  ${CMAKE_CURRENT_SOURCE_DIR}/ressources  )

FILE(COPY doc/ DESTINATION share/doc/)

#variables to be used in both example and validation tests
SET(WaveSystemLinearSolvers_script ${CMAKE_CURRENT_SOURCE_DIR}/examples/WaveSystem/WaveSystem_LinearSolvers/WaveSystemLinearSolvers.py)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

  add_subdirectory (examples) #basic test that run a single calculation
  add_subdirectory (validation) # convergence studies (several runs)

  add_custom_target (CDMATH COMMAND ctest -O testsCDMATH.log)

  add_custom_target (cdmath COMMAND ctest -R TestsCDMATH)
  add_custom_target (medcouplingtests COMMAND ctest -R MEDCoupling)
  add_custom_target (GEOMSMESHtests COMMAND ctest -R GEOMSMESH)
  add_custom_target (check COMMAND ctest -E 'validation|Example')# may be replace ctest -E with ctest -LE
  add_custom_target (vtk COMMAND ctest -R VTK)
  add_custom_target (pv COMMAND ctest -R PV)
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

