project(EOS C CXX)

INCLUDE_DIRECTORIES(
  ${EOS_INCLUDES}
  )

SET(src_eos_CXX
    EosException.cxx
    StiffenedGas.cxx
#   IAPWS97.cxx
  )

ADD_LIBRARY(EosLibs SHARED ${src_eos_CXX})
target_link_libraries( EosLibs )
INSTALL(TARGETS EosLibs DESTINATION lib)

