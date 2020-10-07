
SET(SCRIPT
  ../../scripts/Poisson3DCubeSkinEF/FiniteElements3DPoissonCubeSkin.py
  test_validation3DCubeSkinPoissonEF.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/3DCubeSkin/meshCubeSkin_*" )
file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_3DFE_CUBESKIN ${PYTHON_EXECUTABLE} test_validation3DCubeSkinPoissonEF.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


