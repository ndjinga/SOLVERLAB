
file(COPY ${MED_MESHES} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)
file(COPY scripts/checkerboard_2x2x2_new.msh scripts/locrafgrid_1_new.msh DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)

if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleMEDCoupling_create_weird_mesh ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/weird_mesh.py )

    ADD_TEST(ExampleMEDCoupling_create_brick_wall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DBrickWall/create_mesh_brickWall.py)

    ADD_TEST(ExampleMEDCoupling_convert_locraf_msh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/convert_msh_to_med.py locrafgrid_1_new.msh)

    ADD_TEST(ExampleMEDCoupling_convert_checkerboard_msh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/convert_msh_to_med.py checkerboard_2x2x2_new.msh)

endif (CDMATH_WITH_PYTHON )

install(DIRECTORY ${MED_MESHES} DESTINATION share/meshes)


