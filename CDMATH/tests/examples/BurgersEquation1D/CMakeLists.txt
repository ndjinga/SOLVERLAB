
if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleBurgersEquation_1DFV_shock_formation ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DBurgersEquation_shock_formation.py )
    ADD_TEST(ExampleBurgersEquation_1DFV_rarefaction_wave ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DBurgersEquation_rarefaction_wave.py )

endif (CDMATH_WITH_PYTHON )


