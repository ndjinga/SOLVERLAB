
SET(SCRIPT
  ../scripts/ReferenceSolutions/exact_rs_stiffenedgas.py
  ../scripts/ReferenceSolutions/exact_rs_stiffenedgas_isentropic.py
 ExactRiemannProblemSolutions.py
 ExactRiemannProblemSolutions_isentropic.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationEuler1D_ReferenceSolutions_RiemannProblem            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ExactRiemannProblemSolutions.py )
    ADD_TEST(validationEuler1D_ReferenceSolutions_RiemannProblem_isentropic ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ExactRiemannProblemSolutions_isentropic.py )

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

install(      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/      DESTINATION share/validation/ReferenceSolutions 
              FILES_MATCHING PATTERN "*.dat"
                             PATTERN "*.sh"
                             PATTERN "*.py"
                             PATTERN "*.png"
)

