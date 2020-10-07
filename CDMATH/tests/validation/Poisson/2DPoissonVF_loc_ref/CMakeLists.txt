
SET(SCRIPT
  ../../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
    test_validation2DPoissonVF_loc_ref.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_*" )

SET(NOTEBOOK ../../../doc/2DPoissonVF/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb )

file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFV_loc_ref ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_loc_ref.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


