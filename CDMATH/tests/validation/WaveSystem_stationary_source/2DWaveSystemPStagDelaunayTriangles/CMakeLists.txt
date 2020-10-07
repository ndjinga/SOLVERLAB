
SET(SCRIPT
  ../../scripts/WaveSystemPStag/WaveSystemPStag.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation2DWaveSystemPStag_Delaunay_triangles_source)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_2DFV_PStag_Delaunay_triangles_no_scaling_source ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemPStagDelaunayTriangles_source.py ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_PStag_Delaunay_triangles_scaling_source ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemPStagDelaunayTriangles_source.py ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


