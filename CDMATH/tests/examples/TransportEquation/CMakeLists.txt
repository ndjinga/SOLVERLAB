
if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleTransportEquation_2DUpwind_Explicit ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TransportEquationUpwind.py)

endif (CDMATH_WITH_PYTHON )


