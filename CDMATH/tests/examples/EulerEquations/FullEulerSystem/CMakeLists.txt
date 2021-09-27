file(GLOB FullEulerEquations_EXAMPLES_TO_INSTALL 
EulerSystem1D_HeatedChannel_Roe EulerSystem1D_HeatedChannel_MAC EulerSystem1D_RiemannProblem_Roe EulerSystem1D_RiemannProblem_MAC
)

install(DIRECTORY ${FullEulerEquations_EXAMPLES_TO_INSTALL} DESTINATION share/examples/EulerEquations/FullEulerEquations)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
    ADD_SUBDIRECTORY(EulerSystem1D_RiemannProblem_Roe)
    ADD_SUBDIRECTORY(EulerSystem1D_RiemannProblem_MAC)
    ADD_SUBDIRECTORY(EulerSystem1D_HeatedChannel_Roe)
    ADD_SUBDIRECTORY(EulerSystem1D_HeatedChannel_MAC)
ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


