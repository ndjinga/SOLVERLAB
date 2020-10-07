
SET(SCRIPT
  ../../scripts/Poisson3DVF/FiniteVolumes3DPoisson_CUBE.py
  test_validation3DPoissonVF_tetrahedra.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/3DTetrahedra/meshCubeTetrahedra_*" )
file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_3DFV_tetrahedra ${PYTHON_EXECUTABLE} test_validation3DPoissonVF_tetrahedra.py)


endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


