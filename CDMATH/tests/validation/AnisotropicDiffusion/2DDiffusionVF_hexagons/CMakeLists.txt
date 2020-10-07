
SET(SCRIPT
  ../../scripts/AnisotropicDiffusion2DVF/FiniteVolumes2DDiffusion_SQUARE.py
  test_validation2DDiffusionVF_hexagons.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DHexagons/squareWithHexagons_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationDiffusion_2DFV_hexagons ${PYTHON_EXECUTABLE} test_validation2DDiffusionVF_hexagons.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


