
SET(SCRIPT
  ../../scripts/WaveSystemCentered/WaveSystemCentered.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation2DWaveSystemCentered_squares)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_2DFV_Centered_squares_no_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemCenteredSquares.py ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_Centered_squares_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemCenteredSquares.py ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


