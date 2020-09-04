
INCLUDE_DIRECTORIES(
  ${CoreFlows_INCLUDES}											    #
  )

SET(src_models_CXX
    DriftModel.cxx
    IsothermalTwoFluid.cxx
    SinglePhase.cxx
    FiveEqsTwoFluid.cxx
    ProblemCoreFlows.cxx
    TransportEquation.cxx
    DiffusionEquation.cxx
    StationaryDiffusionEquation.cxx
    Fluide.cxx
    ProblemFluid.cxx
    utilitaire_algebre.cxx
  )

ADD_LIBRARY(CoreFlows SHARED ${src_models_CXX})
target_link_libraries(CoreFlows ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

INSTALL(TARGETS CoreFlows DESTINATION lib)
