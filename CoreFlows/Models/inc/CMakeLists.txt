SET(include_models_HXX
    DriftModel.hxx
    IsothermalTwoFluid.hxx
    SinglePhase.hxx
    FiveEqsTwoFluid.hxx
    ProblemCoreFlows.hxx
    TransportEquation.hxx
    DiffusionEquation.hxx
    StationaryDiffusionEquation.hxx
    LinearElasticityModel.hxx
    Fluide.h
    ProblemFluid.hxx
    utilitaire_algebre.h
    IAPWS-IF97.h
  )
INSTALL(FILES ${include_models_HXX} DESTINATION include)
