file(GLOB TransportEquation_EXAMPLES_TO_INSTALL 
  TransportEquation1DUpwindExplicit TransportEquation1DUpwindImplicit TransportEquation1DCenteredImplicit TransportEquation1DCenteredExplicit # 1D Transport equation 
  TransportEquation
)

install(DIRECTORY ${TransportEquation_EXAMPLES_TO_INSTALL} DESTINATION share/examples/TransportEquation)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_SUBDIRECTORY(TransportEquation1DUpwindExplicit)
    ADD_SUBDIRECTORY(TransportEquation1DUpwindImplicit)
    ADD_SUBDIRECTORY(TransportEquation1DCenteredImplicit)
    ADD_SUBDIRECTORY(TransportEquation1DCenteredExplicit)
    ADD_SUBDIRECTORY(TransportEquation)

ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


