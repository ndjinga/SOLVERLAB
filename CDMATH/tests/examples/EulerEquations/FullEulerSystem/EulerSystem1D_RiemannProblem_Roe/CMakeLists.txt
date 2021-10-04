
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(ISIMPLICIT  0 )#Explicit scheme
    ADD_TEST(ExampleEulerSystem_1DRiemammProblem_ExplicitRoe ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerEquations_RiemannProblem_Roe.py  ${ISIMPLICIT} )

    SET(ISIMPLICIT  1 )#Implicit scheme
    ADD_TEST(ExampleEulerSystem_1DRiemammProblem_ImplicitRoe ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerEquations_RiemannProblem_Roe.py  ${ISIMPLICIT} )

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


