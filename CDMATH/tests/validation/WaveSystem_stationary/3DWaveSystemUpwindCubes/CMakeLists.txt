
SET(SCRIPT
  ../../scripts/WaveSystemUpwind/WaveSystemUpwind.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation3DWaveSystemUpwindCubes)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(BCTYPE  Periodic  )

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_3DFV_Upwind_cubes_no_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindCubes.py ${BCTYPE} ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_3DFV_Upwind_cubes_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindCubes.py ${BCTYPE} ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


