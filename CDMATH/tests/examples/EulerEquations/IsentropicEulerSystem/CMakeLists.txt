file(GLOB IsentropicEulerEquations_EXAMPLES_TO_INSTALL 
EulerSystem1D/EulerSystem1DUpwind EulerSystem1D/EulerSystem1DUpwindEntrCorr EulerSystem1D/EulerSystem1DConservativeStaggered EulerSystem_Shock/EulerSystemStaggered EulerSystem_Shock/EulerSystemUpwind
)

install(DIRECTORY ${IsentropicEulerEquations_EXAMPLES_TO_INSTALL} DESTINATION share/examples/EulerEquations/IsentropicEulerEquations)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
    ADD_SUBDIRECTORY(EulerSystem1D/EulerSystem1DUpwind)
    ADD_SUBDIRECTORY(EulerSystem1D/EulerSystem1DUpwindEntrCorr)
    ADD_SUBDIRECTORY(EulerSystem1D/EulerSystem1DConservativeStaggered)
    ADD_SUBDIRECTORY(EulerSystem_Shock/EulerSystemUpwind)
    ADD_SUBDIRECTORY(EulerSystem_Shock/EulerSystemStaggered)
ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


