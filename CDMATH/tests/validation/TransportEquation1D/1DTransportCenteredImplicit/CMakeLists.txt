
SET(SCRIPT
  ./TransportEquation1DCenteredImplicit.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DTransportEquationCenteredImplicit)

if (CDMATH_WITH_PYTHON )

    SET(CFL  1  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Stiff_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )

    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Smooth_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )

    SET(CFL  2  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Stiff_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )


    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Smooth_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )

endif (CDMATH_WITH_PYTHON )


