
SET(SCRIPT
  ../../scripts/WaveSystemStaggered/WaveSystemStaggered.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation2DWaveSystemStaggered_squares)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_2DFV_Staggered_squares_no_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemStaggeredSquares.py ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_Staggered_squares_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemStaggeredSquares.py ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


