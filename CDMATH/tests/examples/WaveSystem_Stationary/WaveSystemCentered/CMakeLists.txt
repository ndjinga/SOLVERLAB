
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

file(COPY ${MESH_MED} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${MESH_MED} DESTINATION share/examples/WaveSystemCentered)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ../../../ressources/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_3DVortex_Centered_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_3DVortex_Centered_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_DISK_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/diskWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_DISK_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_DISK_spiderWeb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

    SET(MESH_FILE  ../../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DVortex_Centered_DISK_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemCentered.py  ${MESH_FILE})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


