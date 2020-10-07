
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
    test_validation2DPoissonVF_flat_cross_triangles.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DFlatTriangles/squareWithFlatTriangles_*" )

SET(NOTEBOOK ../../../doc/2DPoissonVF/2DPoissonVF_flat_triangles/Convergence_Poisson_FV5_SQUARE_flat_triangles.ipynb )

file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_flat_cross_triangles ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_flat_cross_triangles.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


