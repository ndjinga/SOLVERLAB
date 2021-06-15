
SET(SCRIPT
  ./SinglePhase_1DRiemannProblem.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON )

    SET(SCHEME Upwind )

    install(FILES ${SCRIPT} DESTINATION share/convergence/test_convergenceSinglePhase_1DRiemannProblem_Upwind)

    SET(ISEXPLICIT  1 )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_UpwindExplicit_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME})

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_UpwindExplicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME} )

    SET(ISEXPLICIT  0 )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_UpwindImplicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME} )

    SET(CFL  10  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_UpwindImplicit_CFL10 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME} )

    SET(SCHEME Centered )

    install(FILES ${SCRIPT} DESTINATION share/convergence/test_convergenceSinglePhase_1DRiemannProblem_Centered)

    SET(ISEXPLICIT  0 )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_CenteredImplicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME} )

    SET(CFL  10  )#Courant Friedrichs Lewy number

    ADD_TEST(convergenceSinglePhase_1DRiemannProblem_CenteredImplicit_CFL10 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergenceSinglePhase_1DRiemannProblem.py  ${CFL} ${ISEXPLICIT} ${SCHEME} )

endif (CDMATH_WITH_PYTHON )


