
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    #Cartesian cells split into right triangles
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_rightTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py)
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_rightTriangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_rightTriangles APPEND PROPERTY ENVIRONMENT      "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../ressources/meshSquare.med  )
    SET(MESH_NAME  "DelaunayTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_DelaunayTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_DelaunayTriangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_DelaunayTriangles APPEND PROPERTY ENVIRONMENT      "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../ressources/2DFlatCrossTriangles/squareWithFlatCrossTriangles_0.med  )
    SET(MESH_NAME  "flatCrossTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_flatCrossTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_flatCrossTriangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_flatCrossTriangles APPEND PROPERTY ENVIRONMENT      "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

    SET(MESH_FILE  ../../ressources/2DSkinnyTriangles/squareWithSkinnyTriangles_0.med  )
    SET(MESH_NAME  "skinnyTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_skinnyTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_skinnyTriangles        PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${TESTS_LIBRARY_PATH}:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleSpectrumLaplace_2DEF_SQUARE_skinnyTriangles APPEND PROPERTY ENVIRONMENT      "PYTHONPATH=${TESTS_PYTHON_PATH}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


