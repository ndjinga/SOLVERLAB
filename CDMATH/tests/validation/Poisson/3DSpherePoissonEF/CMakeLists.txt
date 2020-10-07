
SET(SCRIPT
  ../../scripts/Poisson3DSphereEF/FiniteElements3DPoissonSphere.py
  test_validation3DSpherePoissonEF.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/3DSphere/meshSphere_*" )
file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_3DFE_SPHERE ${PYTHON_EXECUTABLE} test_validation3DSpherePoissonEF.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


