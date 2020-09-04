
if (CDMATH_WITH_PYTHON )

    SET(NX 50 )#Number of cells

    SET(SMOOTHNESS  0 )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_CenteredExplicit_Stiff_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredExplicit.py ${NX} ${CFL} ${SMOOTHNESS} )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_CenteredExplicit_Stiff_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredExplicit.py ${NX} ${CFL} ${SMOOTHNESS} )

    SET(SMOOTHNESS  1 )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_CenteredExplicit_Smooth_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredExplicit.py ${NX} ${CFL} ${SMOOTHNESS} )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_CenteredExplicit_Smooth_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredExplicit.py ${NX} ${CFL} ${SMOOTHNESS} )

endif (CDMATH_WITH_PYTHON )


