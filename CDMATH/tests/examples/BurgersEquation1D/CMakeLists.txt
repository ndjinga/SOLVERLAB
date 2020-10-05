
if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleBurgersEquation_1DFV_shock_formation ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DBurgersEquation_shock_formation.py )
    SET_PROPERTY(TEST ExampleBurgersEquation_1DFV_shock_formation        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleBurgersEquation_1DFV_shock_formation APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    ADD_TEST(ExampleBurgersEquation_1DFV_rarefaction_wave ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DBurgersEquation_rarefaction_wave.py )
    SET_PROPERTY(TEST ExampleBurgersEquation_1DFV_rarefaction_wave        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleBurgersEquation_1DFV_rarefaction_wave APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON )


