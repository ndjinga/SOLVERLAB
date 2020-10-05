
if (CDMATH_WITH_PYTHON )

    SET(NX 100 )#Number of cells

    SET(CFL  0.9999  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL1         PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL1  APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL0.5         PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL0.5  APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(CFL  2  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL2         PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Upwind_Explicit_CFL2  APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON )


