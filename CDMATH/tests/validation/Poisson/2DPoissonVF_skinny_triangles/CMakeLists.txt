
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
    test_validation2DPoissonVF_skinny_triangles.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DLongTriangles/squareWithLongTriangles_*" )

SET(NOTEBOOK ../../../doc/2DPoissonVF/2DPoissonVF_long_triangles/Convergence_Poisson_FV5_SQUARE_long_triangles.ipynb )

file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_skinny_triangles ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_skinny_triangles.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


