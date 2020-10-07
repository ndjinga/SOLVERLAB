
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py)

    SET(MESH_FILE  ../../ressources/meshSquare.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/squareWithLocRefSquares.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_loc_ref ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/squareWithCheckerboardSquares.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/squareWithHexagons.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/squareWithBrickWall.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/squareWithDeformedQuadrangles.med  )

    ADD_TEST(ExamplePoisson_2DVF_SQUARE_deformed_quadrangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_SQUARE.py ${MESH_FILE})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


