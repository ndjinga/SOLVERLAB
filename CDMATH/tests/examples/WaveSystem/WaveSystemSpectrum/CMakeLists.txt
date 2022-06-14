
if (CDMATH_WITH_PYTHON )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleSpectrumWaveSystem_2DUpwind   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemSpectrum.py  ${MESH_FILE} )
    ADD_TEST(ExampleSpectrumWaveSystem_2DCentered ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemSpectrum.py  ${MESH_FILE} "Centred")
    ADD_TEST(ExampleSpectrumWaveSystem_2DPStag    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemSpectrum.py  ${MESH_FILE} "PStag")

endif (CDMATH_WITH_PYTHON )


