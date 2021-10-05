
SET(SCRIPT
  ../../scripts/1DEulerEquations/EulerEquations1D.py
  ../../scripts/ReferenceSolutions/exact_rs_stiffenedgas.py
  ./test_validation1DEulerEquationsToro.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DEulerEquationsToro)

if (CDMATH_WITH_PYTHON )

    SET(SCHEME  Roe  )#Roe numerical method

    SET(ISIMPLICIT  0 )#Explicit scheme

    ADD_TEST(validation1DEulerEquationsToro_Roe_Explicit ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DEulerEquationsToro.py  ${SCHEME} ${ISIMPLICIT} )

    SET(ISIMPLICIT  1 )#Implicit scheme

    ADD_TEST(validation1DEulerEquationsToro_Roe_Implicit ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DEulerEquationsToro.py  ${SCHEME} ${ISIMPLICIT} )

    SET(SCHEME  Stag  )#Conservative staggered numerical method

    SET(ISIMPLICIT  1 )#Implicit scheme

    ADD_TEST(validation1DEulerEquationsToro_Stag_Implicit ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DEulerEquationsToro.py  ${SCHEME} ${ISIMPLICIT} )


endif (CDMATH_WITH_PYTHON )

install(      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/    DESTINATION share/validation/test_validation1DEulerEquationsToro
              FILES_MATCHING PATTERN "*.py"
                             PATTERN "*.png"
)

