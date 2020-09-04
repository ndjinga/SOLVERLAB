SET(include_models_HXX
    DriftModel.hxx
    IsothermalTwoFluid.hxx
    SinglePhase.hxx
    FiveEqsTwoFluid.hxx
    ProblemCoreFlows.hxx
    TransportEquation.hxx
    DiffusionEquation.hxx
    StationaryDiffusionEquation.hxx
    Fluide.h
    ProblemFluid.hxx
    utilitaire_algebre.h
  )
INSTALL(FILES ${include_models_HXX} DESTINATION include)
