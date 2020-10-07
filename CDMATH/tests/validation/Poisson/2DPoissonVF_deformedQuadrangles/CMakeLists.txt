
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
  test_validation2DPoissonVF_deformedQuadrangles.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DDeformedQuadrangles/squareWithDeformedQuadrangles_*" )

SET(NOTEBOOK ../../../doc/2DPoissonVF/2DPoissonVF_deformedQuadrangles/Convergence_Poisson_FV5_SQUARE_deformedQuadrangles.ipynb )

file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_deformed_quadrangles ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_deformedQuadrangles.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


