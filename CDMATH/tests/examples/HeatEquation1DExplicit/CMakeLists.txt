
if (CDMATH_WITH_PYTHON )

    SET(NX 100 )#Number of cells

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleHeatEquation_1DFV_Explicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/HeatEquation1DExplicit.py ${NX} ${CFL})

    SET(CFL  2  )#Courant Friedrichs Lewy number

    ADD_TEST(ExampleHeatEquation_1DFV_Explicit_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/HeatEquation1DExplicit.py ${NX} ${CFL})

endif (CDMATH_WITH_PYTHON )


