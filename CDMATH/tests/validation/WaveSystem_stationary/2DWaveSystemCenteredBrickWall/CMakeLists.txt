
SET(SCRIPT
  ../../scripts/WaveSystemCentered/WaveSystemCentered.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DBrickWall/squareWithBrickWall_*" )

FILE(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION share/validation/test_validation2DWaveSystemCentered_brickwall)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(SCALING  2  )

    ADD_TEST(validationWaveSystem_2DFV_Centered_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DWaveSystemCenteredBrickWall.py ${SCALING})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


