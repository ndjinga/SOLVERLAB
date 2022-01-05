
SET(MESH_MED
  ../ressources/meshSquare.med
  ../ressources/meshSphere.med
  ../ressources/meshCube.med
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/swig)
ADD_TEST(UnitTestsCDMATH.py ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/mainTests.py)

if (CDMATH_WITH_PETSC)
    ADD_TEST(TestsLinearSolver.py ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/linearSolverTests.py)
    INSTALL(FILES linearSolverTests.py DESTINATION bin/unit_tests)
endif (CDMATH_WITH_PETSC)


if (CDMATH_WITH_POSTPRO)
  ADD_TEST(Tests_VTK_routines ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_VTK_routines.py)

  ADD_TEST(Tests_PV_routines ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_PV_routines.py)

  add_subdirectory (BoySurface)


  ADD_TEST(Tests_Field_Creation_And_Saving_routines ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testFieldCreationAndSave.py )

  INSTALL(FILES ${MESH_MED}  mainTests.py test_VTK_routines.py test_PV_routines.py testFieldCreationAndSave.py DESTINATION bin/unit_tests)
endif (CDMATH_WITH_POSTPRO)

