project(Models C CXX)

INCLUDE_DIRECTORIES(
  ${EOS_INCLUDES}
  ${MODELS_INCLUDES}
  ${CDMATH_INCLUDES}
  ${PETSC_INCLUDES}
  )

SET(src_models_CXX
    DriftModel.cxx
    IsothermalTwoFluid.cxx
    IsothermalSinglePhase.cxx
    SinglePhase.cxx
    FiveEqsTwoFluid.cxx
    ProblemCoreFlows.cxx
    TransportEquation.cxx
    DiffusionEquation.cxx
    StationaryDiffusionEquation.cxx
#    LinearElasticityModel.cxx
    ProblemFluid.cxx
    utilitaire_algebre.cxx
  )

ADD_LIBRARY(CoreFlowsLibs SHARED ${src_models_CXX})
target_link_libraries(CoreFlowsLibs ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} EosLibs)
INSTALL(TARGETS CoreFlowsLibs DESTINATION lib)

ADD_SUBDIRECTORY(EOS)

