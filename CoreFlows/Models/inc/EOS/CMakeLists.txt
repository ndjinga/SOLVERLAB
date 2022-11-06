SET(include_eos_HXX
    EosException.hxx
    Fluide.h
    StiffenedGas.hxx
#   IAPWS97.hxx
)

INSTALL(FILES ${include_eos_HXX} DESTINATION include)
