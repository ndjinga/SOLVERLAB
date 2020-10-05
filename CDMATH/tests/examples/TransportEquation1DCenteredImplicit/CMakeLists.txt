
if (CDMATH_WITH_PYTHON )

    SET(NX 50 )#Number of cells

    SET(CFL  1  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Centered_Implicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredImplicit.py ${NX} ${CFL})
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Centered_Implicit_CFL1         PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Centered_Implicit_CFL1  APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(CFL  2  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Centered_Implicit_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationCenteredImplicit.py ${NX} ${CFL})
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Centered_Implicit_CFL2         PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleTransportEquation_1DFV_Centered_Implicit_CFL2  APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")


endif (CDMATH_WITH_PYTHON )


