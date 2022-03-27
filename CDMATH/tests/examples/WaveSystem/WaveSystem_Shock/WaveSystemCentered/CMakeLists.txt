
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_Centered_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_3DShock_Centered_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/diskWithTriangles.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/diskWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/diskWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ${MED_MESHES}/meshHexagonWithTriangles.med  )

    ADD_TEST(ExampleWaveSystem_2DShock_Centered_HEXAGON_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


