
SET(SCRIPT
  ../../scripts/Euler1DIsothermal/EulerIsothermal1DSchemeComparison.py
  test_validationEulerIsothermal_1D_mesh_refinement.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationEulerEquations_1D_mesh_refinement ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validationEulerIsothermal_1D_mesh_refinement.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


