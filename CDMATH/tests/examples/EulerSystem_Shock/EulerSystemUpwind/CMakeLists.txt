
SET(MESH_MED
  ../../../ressources/squareWithTriangles.med
  ../../../ressources/meshCube.med
  ../../../ressources/squareWithSquares.med
  ../../../ressources/cubeWithCubes.med
  ../../../ressources/diskWithTriangles.med
  ../../../ressources/diskWithSquares.med
  ../../../ressources/diskWithSpiderWeb.med
  ../../../ressources/diskWithHexagons.med
  ../../../ressources/ballWithTetrahedra.med
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/tests/examples/EulerSystemUpwind)
install(FILES ${MESH_MED} DESTINATION share/examples/EulerSystemUpwind)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(IMPLICIT_SCHEME  0 )

    SET(MESH_FILE  ../../../ressources/meshSquare.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/meshCube.med  )

    ADD_TEST(ExampleEulerSystem_3DShock_UpwindExplicit_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/cubeWithCubes.med  )

    ADD_TEST(ExampleEulerSystem_3DShock_UpwindExplicit_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithSquares.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/meshHexagonWithTriangles.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindExplicit_HEXAGON_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(IMPLICIT_SCHEME  1 )

    SET(MESH_FILE  ../../../ressources/meshSquare.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/meshCube.med  )

    ADD_TEST(ExampleEulerSystem_3DShock_UpwindImplicit_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/cubeWithCubes.med  )

    ADD_TEST(ExampleEulerSystem_3DShock_UpwindImplicit_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithSquares.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

    SET(MESH_FILE  ../../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExampleEulerSystem_2DShock_UpwindImplicit_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


