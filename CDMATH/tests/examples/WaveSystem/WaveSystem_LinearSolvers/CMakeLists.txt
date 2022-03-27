
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    SET(LIN_SOLVER  gmres  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(ExampleWaveSystem_GMRES_Centered_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(ExampleWaveSystem_GMRES_Upwind_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  PStag  )

    ADD_TEST(ExampleWaveSystem_GMRES_PStag_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(LIN_SOLVER  cgne  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(ExampleWaveSystem_CGNE_Centered_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(ExampleWaveSystem_CGNE_Upwind_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  PStag  )

    ADD_TEST(ExampleWaveSystem_CGNE_PStag_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(LIN_SOLVER  lu  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(ExampleWaveSystem_LU_Centered_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(ExampleWaveSystem_LU_Upwind_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

    SET(NUM_METHOD  PStag  )

    ADD_TEST(ExampleWaveSystem_LU_PStag_triangles ${PYTHON_EXECUTABLE} ${WaveSystemLinearSolvers_script}  ${MESH_FILE} ${NUM_METHOD} ${LIN_SOLVER})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


