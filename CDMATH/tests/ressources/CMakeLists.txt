
file(COPY ${MED_MESHES} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)
file(COPY 
ConversionScripts/FVCA6/checkerboard_2x2x2_new.msh ConversionScripts/FVCA6/locrafgrid_1_new.msh ConversionScripts/FVCA6/RandMesh4.msh ConversionScripts/FVCA6/tet.00.msh ConversionScripts/FVCA6/WellMesh_1.msh
ConversionScripts/FVCA5/mesh1_1.typ2 ConversionScripts/FVCA5/mesh2_1.typ2 ConversionScripts/FVCA5/mesh3_1.typ2 ConversionScripts/FVCA5/mesh4_1.typ2 ConversionScripts/FVCA5/mesh4_2.typ2
DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)

if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleMEDCoupling_create_weird_mesh ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/weird_mesh.py )

    ADD_TEST(ExampleMEDCoupling_create_brick_wall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DBrickWall/create_mesh_brickWall.py)

    ADD_TEST(ExampleMEDCoupling_convert_LocRef_3Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_msh_to_med.py locrafgrid_1_new.msh)

    ADD_TEST(ExampleMEDCoupling_convert_Checkerboard_3Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_msh_to_med.py checkerboard_2x2x2_new.msh)

    ADD_TEST(ExampleMEDCoupling_convert_RandMesh_3Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_msh_to_med.py RandMesh4.msh)

    ADD_TEST(ExampleMEDCoupling_convert_Tetra_3Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_msh_to_med.py tet.00.msh)

    ADD_TEST(ExampleMEDCoupling_convert_Well_3Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_msh_to_med.py WellMesh_1.msh)

    ADD_TEST(ExampleMEDCoupling_convert_Triangle_2Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_2Dmsh_to_med.py mesh1_1.typ2)

    ADD_TEST(ExampleMEDCoupling_convert_Cart_2Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_2Dmsh_to_med.py mesh2_1.typ2)

    ADD_TEST(ExampleMEDCoupling_convert_LocRef_2Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_2Dmsh_to_med.py mesh3_1.typ2)

    ADD_TEST(ExampleMEDCoupling_convert_Kershaw1_2Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_2Dmsh_to_med.py mesh4_1.typ2)

    ADD_TEST(ExampleMEDCoupling_convert_Kershaw2_2Dmsh_to_med ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/convert_2Dmsh_to_med.py mesh4_2.typ2)

endif (CDMATH_WITH_PYTHON )

install(DIRECTORY ${MED_MESHES} DESTINATION share/meshes)


