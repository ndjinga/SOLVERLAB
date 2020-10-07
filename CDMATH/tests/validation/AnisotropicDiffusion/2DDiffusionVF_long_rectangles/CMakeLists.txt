
SET(SCRIPT
  ../../scripts/AnisotropicDiffusion2DVF/FiniteVolumes2DDiffusion_SQUARE.py
    test_validation2DDiffusionVF_long_rectangles.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DLongRectangles/squareWithLongRectangles_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationDiffusion_2DFV_long_rectangles ${PYTHON_EXECUTABLE} test_validation2DDiffusionVF_long_rectangles.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


