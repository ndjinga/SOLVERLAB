
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
    test_validation2DPoissonVF_right_triangles.py
  )
file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_right_triangles ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_right_triangles.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


