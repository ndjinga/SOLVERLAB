
SET(MESH_MED
  ${MED_MESHES}/meshCube.med
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${MESH_MED} DESTINATION share/examples/Poisson3DEF)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(ExamplePoisson_3DEF_CUBE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteElements3DPoisson_CUBE.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


