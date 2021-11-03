
file(COPY ${MED_MESHES} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)
file(COPY scripts/mesh_ref_1.typ2 scripts/mesh_tetra_0.typ3 scripts/mesh_tetra_0.msh DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)

if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleMEDCoupling_create_weird_mesh ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/weird_mesh.py )

    ADD_TEST(ExampleMEDCoupling_create_brick_wall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DBrickWall/create_mesh_brickWall.py)

    ADD_TEST(ExampleMEDCoupling_convert_typ2_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/convert_typ_to_med.py mesh_ref_1.typ2)

    ADD_TEST(ExampleMEDCoupling_convert_typ3_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/convert_typ_to_med.py mesh_tetra_0.typ3)

    ADD_TEST(ExampleMEDCoupling_convert_msh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/convert_gmsh_to_med.py mesh_tetra_0.msh)

endif (CDMATH_WITH_PYTHON )

install(DIRECTORY ${MED_MESHES} DESTINATION share/meshes)


