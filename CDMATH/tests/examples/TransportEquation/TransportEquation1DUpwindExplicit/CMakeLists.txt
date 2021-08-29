
if (CDMATH_WITH_PYTHON )

    SET(NX 100 )#Number of cells

    SET(CFL  0.9999  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})

    SET(CFL  2  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleTransportEquation_1DFV_Upwind_Explicit_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DTransportEquationUpwindExplicit.py ${NX} ${CFL})

endif (CDMATH_WITH_PYTHON )


