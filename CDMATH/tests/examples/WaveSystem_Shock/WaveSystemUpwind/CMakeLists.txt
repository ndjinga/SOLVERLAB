
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

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/tests/examples/WaveSystemUpwind)
install(FILES ${MESH_MED} DESTINATION share/examples/WaveSystemUpwind)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(IMPLICIT_SCHEME  0 )

    SET(MESH_FILE  ../../../ressources/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_triangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_squares        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_SQUARE_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_tetrahedra        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_tetrahedra APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_cubes        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindExplicit_CUBE_cubes APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_triangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_squares        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_spiderWeb        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_spiderWeb APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindExplicit_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_hexagons        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindExplicit_DISK_hexagons APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(IMPLICIT_SCHEME  1 )

    SET(MESH_FILE  ../../../ressources/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_triangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_squares        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_SQUARE_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_tetrahedra        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_tetrahedra APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_cubes        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_3DShock_UpwindImplicit_CUBE_cubes APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_triangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_squares        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_spiderWeb        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_spiderWeb APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_UpwindImplicit_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_hexagons        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DShock_UpwindImplicit_DISK_hexagons APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


