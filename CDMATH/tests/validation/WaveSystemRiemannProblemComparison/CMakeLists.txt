
SET(SCRIPT
  ..//scripts/WaveSystemRiemannProblemComparison/WaveSystemRiemannProblemComparison.py
  test_validationWaveSystemRiemannProblem_1D_Increase_Sound_Speed.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationWaveSystemRiemannProblemComparison_Increase_Sound_Speed ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validationWaveSystemRiemannProblem_1D_Increase_Sound_Speed.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


