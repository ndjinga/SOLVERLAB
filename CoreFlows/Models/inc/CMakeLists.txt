SET(include_models_HXX
    DriftModel.hxx
    IsothermalTwoFluid.hxx
    IsothermalSinglePhase.hxx
    SinglePhase.hxx
    FiveEqsTwoFluid.hxx
    ProblemCoreFlows.hxx
    TransportEquation.hxx
    DiffusionEquation.hxx
    StationaryDiffusionEquation.hxx
    LinearElasticityModel.hxx
    ProblemFluid.hxx
    utilitaire_algebre.h
  )

ADD_SUBDIRECTORY(EOS)

INSTALL(FILES ${include_models_HXX} DESTINATION include)
