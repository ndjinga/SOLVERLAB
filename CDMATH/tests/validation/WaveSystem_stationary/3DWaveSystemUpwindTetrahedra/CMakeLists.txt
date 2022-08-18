
SET(SCRIPT
  ../../scripts/WaveSystemUpwind/WaveSystemUpwind.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/3DTetrahedra/meshCubeTetrahedra_*" )
file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION share/validation/test_validation3DWaveSystemUpwindTetrahedra)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(BCTYPE  Wall  )

    SET(SCALING  0  )

    ADD_TEST(validationWaveSystem_3DFV_Upwind_tetrahedra_no_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindTetrahedra.py ${BCTYPE} ${SCALING})

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_3DFV_Upwind_tetrahedra_scaling ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindTetrahedra.py ${BCTYPE} ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


