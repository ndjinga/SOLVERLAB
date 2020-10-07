
SET(SCRIPT
  ../../scripts/Poisson2DVF_StiffBC_SQUARE/FiniteVolumes2DPoissonStiffBC_SQUARE.py
  test_validation2DPoissonVF_StiffBCSquare_squares.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DCartesien/squareWithSquares_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_StiffBCSquare_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DPoissonVF_StiffBCSquare_squares.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


