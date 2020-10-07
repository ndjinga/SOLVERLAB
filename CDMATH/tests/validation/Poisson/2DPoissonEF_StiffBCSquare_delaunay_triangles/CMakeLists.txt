
SET(SCRIPT
  ../../scripts/Poisson2DEF_StiffBC_SQUARE/FiniteElements2DPoissonStiffBC_SQUARE.py
  test_validation2DPoissonEF_StiffBCSquare_Delaunay_triangles.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DTriangles/squareWithTriangles_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFE_StiffBCSquare_Delaunay_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DPoissonEF_StiffBCSquare_Delaunay_triangles.py)



endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


