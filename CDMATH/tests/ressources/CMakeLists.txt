
file(COPY ${MED_MESHES} DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)
file(COPY 
ConversionScripts/FVCA6/checkerboard_2x2x2_new.msh ConversionScripts/FVCA6/locrafgrid_1_new.msh ConversionScripts/FVCA6/RandMesh4.msh ConversionScripts/FVCA6/tet.00.msh ConversionScripts/FVCA6/WellMesh_1.msh
ConversionScripts/FVCA5/mesh1_1.typ2 ConversionScripts/FVCA5/mesh2_1.typ2 ConversionScripts/FVCA5/mesh3_1.typ2 ConversionScripts/FVCA5/mesh4_1.typ2 ConversionScripts/FVCA5/mesh4_2.typ2
DESTINATION ${CMAKE_BINARY_DIR}/CDMATH/tests/ressources)

if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleMEDCoupling_create_2D_brick_wall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DBrickWall/create_mesh_brickWall.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_cartesian ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DCartesien/create_mesh_of_squares_with_squares.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_Checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DCheckerboard/create_checkerboard_2d_quad_meshes_amr.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_CrossTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DCrossTriangles/create_mesh_of_square_with_cross_triangles.py)
    # test 2DdiskWithHexagons generate error pb with slice
    ADD_TEST(ExampleMEDCoupling_create_2D_DiskWithHexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DDiskWithHexagons/create_mesh_diskWithHexagons.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_disk_spider_web ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DDiskSpiderWeb/create_mesh_diskSpiderWeb.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_diskWithSquares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DdiskWithSquares/create_mesh_diskWithSquares.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_EquilateralTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DEquilateralTriangles/create_mesh_squareWithEquilateralTriangles.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_FlatCrossTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DFlatCrossTriangles/create_mesh_of_square_with_cross_triangles.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_Hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DHexagons/create_mesh_squareWithHexagons.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_LocRefinedSquares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DLocRefinedSquares/create_locally_refined_2d_quad_meshes_amr.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_LongRectangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DLongRectangles/create_mesh_of_rectangle_with_rectangles.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_RightTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DRightTriangles/create_mesh_right_triangles.py)

    ADD_TEST(ExampleMEDCoupling_create_2D_SkinnyTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DSkinnyTriangles/create_mesh_of_rectangle_with_right_triangles.py)

    ADD_TEST(ExampleMEDCoupling_create_3D_Hexahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/3DHexahedra/create_mesh_of_cube_with_cuboids.py)

    ADD_TEST(ExampleMEDCoupling_create_3DLocRefinedCubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/3DLocRefinedCubes/create_locally_refined_3d_quad_meshes_amr.py)

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

    ADD_TEST(ExampleMEDCoupling_create_weird_mesh ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConversionScripts/weird_mesh.py )

    IF(EXISTS ${KERNEL_ROOT_DIR})
        ADD_TEST(ExampleGEOMSMESH_create_2D_NACA_profile ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DNACA_profile/naca_profile.py)

        ADD_TEST(ExampleGEOMSMESH_create_2D_HexagonWithTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2DHexagonWithTriangles/HexagonMesh.py)

        ADD_TEST(ExampleGEOMSMESH_create_3DTetrahedron ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/3DTetrahedron/Tetrahedron.py)
    ENDIF(EXISTS ${KERNEL_ROOT_DIR})
endif (CDMATH_WITH_PYTHON )

install(DIRECTORY ${MED_MESHES} DESTINATION share/meshes)


