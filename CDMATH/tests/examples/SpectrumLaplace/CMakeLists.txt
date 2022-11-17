file(GLOB SpectrumLaplace_EXAMPLES_TO_INSTALL 
  SpectrumLaplace2DEF  SpectrumLaplace2DVF  SpectrumLaplaceBeltrami3DEF SpectrumLaplace3DEF_Neumann #Spectrum of Laplace operator
)

install(DIRECTORY ${SpectrumLaplace_EXAMPLES_TO_INSTALL} DESTINATION share/examples/SpectrumLaplace)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_SUBDIRECTORY(SpectrumLaplace2DEF)
    ADD_SUBDIRECTORY(SpectrumLaplace2DVF)
    ADD_SUBDIRECTORY(SpectrumLaplaceBeltrami3DEF)
    ADD_SUBDIRECTORY(SpectrumLaplace3DEF_Neumann)

ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


