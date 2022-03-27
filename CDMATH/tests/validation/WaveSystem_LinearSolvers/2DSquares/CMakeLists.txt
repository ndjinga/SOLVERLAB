
file(COPY ${WaveSystemLinearSolvers_script} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${WaveSystemLinearSolvers_script} DESTINATION share/validation/test_validationWaveSystemLinearSolvers_Squares)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC )

    SET(CFL  10 )

    SET(LIN_SOLVER  gmres  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(validationWaveSystem_Centered_GMRES_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(validationWaveSystem_Upwind_GMRES_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

    SET(LIN_SOLVER  cgne  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(validationWaveSystem_Centered_CGNE_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(validationWaveSystem_Upwind_CGNE_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

    SET(LIN_SOLVER  lu  )

    SET(NUM_METHOD  centered  )

    ADD_TEST(validationWaveSystem_Centered_LU_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

    SET(NUM_METHOD  upwind  )

    ADD_TEST(validationWaveSystem_Upwind_LU_Squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemLinearSolvers_Squares.py ${CFL} ${NUM_METHOD} ${LIN_SOLVER} ${MED_MESHES})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC )


