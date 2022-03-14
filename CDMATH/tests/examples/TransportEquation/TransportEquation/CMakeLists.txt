
if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleTransportEquation_2DUpwind_Explicit ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TransportEquationUpwind.py)
    ADD_TEST(ExampleSpectrumTransportEquation_2DUpwind ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumTransportEquation.py)
    ADD_TEST(ExampleSpectrumTransportEquation_2DCentered ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumTransportEquation.py "")

endif (CDMATH_WITH_PYTHON )


