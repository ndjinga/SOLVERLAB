
INCLUDE_DIRECTORIES(
  ${CoreFlows_INCLUDES}											    #
  ${PETSC_INCLUDES}
  ${SLEPC_INCLUDES}
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
#    LinearElasticityModel.cxx
    ProblemFluid.cxx
    utilitaire_algebre.cxx
    EOS/StiffenedGas.cxx
#    EOS/IAPWS97.cxx
  )

ADD_LIBRARY(CoreFlowsLibs SHARED ${src_models_CXX})
target_link_libraries(CoreFlowsLibs ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

INSTALL(TARGETS CoreFlowsLibs DESTINATION lib)
