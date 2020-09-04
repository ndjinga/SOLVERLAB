INCLUDE_DIRECTORIES(
  ${SLEPC_INCLUDES}
)

FILE(GLOB include_linearsolver_HXX . *.hxx  )

INSTALL(FILES ${include_linearsolver_HXX} DESTINATION include)
