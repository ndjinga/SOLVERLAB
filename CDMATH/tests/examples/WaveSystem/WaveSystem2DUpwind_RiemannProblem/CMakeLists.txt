
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(IMPLICIT_SCHEME  0 )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithBrickWall.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithCheckerboardSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithDeformedQuadrangles.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_deformedQuadrangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(IMPLICIT_SCHEME  1 )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithBrickWall.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithCheckerboardSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithDeformedQuadrangles.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_deformedQuadrangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ${MED_MESHES}/squareWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


