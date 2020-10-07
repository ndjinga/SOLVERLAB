
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
  test_validation2DPoissonVF_hexagons.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DHexagons/squareWithHexagons_*" )

SET(NOTEBOOK ../../../doc/2DPoissonVF/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb )

file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_hexagons ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_hexagons.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


