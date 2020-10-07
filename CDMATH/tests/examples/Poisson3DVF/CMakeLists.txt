
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(ExamplePoisson_3DVF_CUBE_cubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes3DPoisson_CUBE.py)

    SET(MESH_FILE  ../../ressources/meshCube.med  )

    ADD_TEST(ExamplePoisson_3DVF_CUBE_tetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes3DPoisson_CUBE.py)

    SET(MESH_FILE  ../../ressources/cubeWithLocRefCubes.med  )

    ADD_TEST(ExamplePoisson_3DVF_CUBE_loc_ref ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes3DPoisson_CUBE.py)

    SET(MESH_FILE  ../../ressources/3DCheckerboard/checkerboard_4x4x4.med  )

    ADD_TEST(ExamplePoisson_3DVF_CUBE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolumes3DPoisson_CUBE.py)

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


