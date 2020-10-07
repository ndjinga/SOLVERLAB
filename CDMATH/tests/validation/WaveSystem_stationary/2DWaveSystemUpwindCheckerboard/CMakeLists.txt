
SET(SCRIPT
  ../../scripts/WaveSystemUpwind/WaveSystemUpwind.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation2DWaveSystemUpwindCheckerboard)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(BCTYPE  Periodic  )

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_Upwind_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemUpwindCheckerboard.py ${BCTYPE} ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


