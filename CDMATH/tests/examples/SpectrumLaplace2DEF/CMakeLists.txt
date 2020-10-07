
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    #Cartesian cells split into right triangles
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_rightTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py)

    SET(MESH_FILE  ../../ressources/meshSquare.med  )
    SET(MESH_NAME  "DelaunayTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_DelaunayTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})

    SET(MESH_FILE  ../../ressources/2DFlatCrossTriangles/squareWithFlatCrossTriangles_0.med  )
    SET(MESH_NAME  "flatCrossTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_flatCrossTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})

    SET(MESH_FILE  ../../ressources/2DSkinnyTriangles/squareWithSkinnyTriangles_0.med  )
    SET(MESH_NAME  "skinnyTriangles" )
    ADD_TEST(ExampleSpectrumLaplace_2DEF_SQUARE_skinnyTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace2DEF_SQUARE.py ${MESH_FILE} ${MESH_NAME})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


