
SET(SCRIPT
  ../SinglePhase_1DRiemannProblem.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validationSinglePhase_1DRiemannProblem_UpwindExplicit)

if (CDMATH_WITH_PYTHON )

    SET(ISEXPLICIT  1 )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(validationSinglePhase_1DRiemannProblem_UpwindExplicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validationSinglePhase_1DRiemannProblem_UpwindExplicit.py  ${CFL} ${ISEXPLICIT} )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(validationSinglePhase_1DRiemannProblem_UpwindExplicit_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validationSinglePhase_1DRiemannProblem_UpwindExplicit.py  ${CFL} ${ISEXPLICIT} )

endif (CDMATH_WITH_PYTHON )


