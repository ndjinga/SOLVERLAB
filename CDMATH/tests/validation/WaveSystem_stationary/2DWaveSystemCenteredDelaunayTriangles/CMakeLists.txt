
SET(SCRIPT
  ../../scripts/WaveSystemCentered/WaveSystemCentered.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DTriangles/squareWithTriangles_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION share/validation/test_validation2DWaveSystemCentered_Delaunay_triangles)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_2DFV_Centered_Delaunay_triangles_no_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemCenteredDelaunayTriangles.py ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_Centered_Delaunay_triangles_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemCenteredDelaunayTriangles.py ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


