
SET(MESH_MED
  ./Mesh_RadiatorAndWindow.med
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${MESH_MED} DESTINATION share/examples/Poisson3DEF_RadiatorAndWindow)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(ExamplePoisson_3DEF_CUBE_RadiatorAndWindow ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteElements3DPoisson_CUBE_RadiatorAndWindow.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


