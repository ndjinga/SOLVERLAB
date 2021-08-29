
if (CDMATH_WITH_PYTHON )

    SET(NX 50 )#Number of cells

    SET(CFL  1  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Implicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindImplicit.py ${NX} ${CFL})

    SET(CFL  2  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Implicit_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindImplicit.py ${NX} ${CFL})


endif (CDMATH_WITH_PYTHON )


