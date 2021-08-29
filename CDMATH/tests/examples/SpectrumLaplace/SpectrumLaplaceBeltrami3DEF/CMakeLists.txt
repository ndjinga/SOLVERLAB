
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ${MED_MESHES}/meshSphere.med  )
    SET(MESH_NAME  "Sphere" )

    ADD_TEST(ExampleSpectrumLaplaceBeltrami_3DFE_SPHERE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumFiniteElements3DLaplace-Beltrami.py ${MESH_FILE} ${MESH_NAME} )

    SET(MESH_FILE  ${MED_MESHES}/meshTorus.med  )
    SET(MESH_NAME  "Torus" )

    ADD_TEST(ExampleSpectrumLaplaceBeltrami_3DFE_TORUS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumFiniteElements3DLaplace-Beltrami.py ${MESH_FILE} ${MESH_NAME} )

    SET(MESH_FILE  ${MED_MESHES}/meshCubeSkin.med  )
    SET(MESH_NAME  "CubeSkin" )

    ADD_TEST(ExampleSpectrumLaplaceBeltrami_3DFE_CUBESKIN ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumFiniteElements3DLaplace-Beltrami.py ${MESH_FILE} ${MESH_NAME} )

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


