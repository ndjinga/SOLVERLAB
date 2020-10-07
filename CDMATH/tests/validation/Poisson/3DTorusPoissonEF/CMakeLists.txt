
SET(SCRIPT
  ../../scripts/Poisson3DTorusEF/FiniteElements3DPoissonTorus.py
  test_validation3DTorusPoissonEF.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/3DTorus/meshTorus_*" )
file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_3DFE_TORUS ${PYTHON_EXECUTABLE} test_validation3DTorusPoissonEF.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


