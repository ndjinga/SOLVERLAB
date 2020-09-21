SET(MESH_MED2
  ../ressources/
  )

file(COPY ${MESH_MED2} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)

install(DIRECTORY ${MESH_MED2} DESTINATION share/meshes)


