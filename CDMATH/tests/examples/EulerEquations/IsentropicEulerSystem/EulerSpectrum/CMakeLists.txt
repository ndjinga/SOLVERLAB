
if (CDMATH_WITH_PYTHON )

    SET(MESH_FILE  ${MED_MESHES}/meshSquare.med  )

    ADD_TEST(ExampleSpectrumIsentropicEulerSystem_2DUpwind   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSpectrum.py  ${MESH_FILE} )
    ADD_TEST(ExampleSpectrumIsentropicEulerSystem_2DCentered ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSpectrum.py  ${MESH_FILE} "Centred")
    ADD_TEST(ExampleSpectrumIsentropicEulerSystem_2DPStag    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/EulerSpectrum.py  ${MESH_FILE} "PStag")

endif (CDMATH_WITH_PYTHON )


