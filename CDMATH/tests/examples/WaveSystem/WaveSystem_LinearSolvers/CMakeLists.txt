
file(COPY ${MESH_MED} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${MESH_MED} DESTINATION share/examples/WaveSystemLinearSolvers)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(NUM_METHOD  centered  )

    SET(LIN_SOLVER  gmres  )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_GMRES_Centered_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_GMRES_Centered_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_GMRES_Centered_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_GMRES_Centered_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(LIN_SOLVER  cgne  )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_CGNE_Centered_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_CGNE_Centered_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_CGNE_Centered_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_CGNE_Centered_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(LIN_SOLVER  lu  )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_LU_Centered_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_LU_Centered_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/meshCube.med  )

    ADD_TEST(ExampleWaveSystem_LU_Centered_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(MESH_FILE  ${MED_MESHES}/cubeWithCubes.med  )

    ADD_TEST(ExampleWaveSystem_LU_Centered_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemLinearSolvers.py  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


