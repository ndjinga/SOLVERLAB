
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ../../ressources/diskWithSquares.med  )

    ADD_TEST(ExamplePoisson_2DVF_DISK_StiffBC_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_DISK_StiffBC.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/diskWithTriangles.med  )

    ADD_TEST(ExamplePoisson_2DVF_DISK_StiffBC_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_DISK_StiffBC.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/diskWithHexagons.med  )

    ADD_TEST(ExamplePoisson_2DVF_DISK_StiffBC_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_DISK_StiffBC.py ${MESH_FILE})

    SET(MESH_FILE  ../../ressources/diskWithSpiderWeb.med  )

    ADD_TEST(ExamplePoisson_2DVF_DISK_StiffBC_spiderweb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes2DPoisson_DISK_StiffBC.py ${MESH_FILE})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


