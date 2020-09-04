
SET(SCRIPT
  ./TransportEquation1DUpwindExplicit.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DTransportEquationUpwindExplicit)

if (CDMATH_WITH_PYTHON )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_UpwindExplicit_Stiff_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )

    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_UpwindExplicit_Smooth_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_UpwindExplicit_Stiff_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )


    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_UpwindExplicit_Smooth_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )

endif (CDMATH_WITH_PYTHON )


