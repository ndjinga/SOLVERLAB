
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(MESH_FILE  ${MED_MESHES}/ballWithTetrahedra.med  )
    SET(MESH_NAME  "Tetrahedra" )
    ADD_TEST(ExampleSpectrumLaplace_3DEF_BALL_Neumann ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpectrumLaplace3DEF.py ${MESH_FILE} ${MESH_NAME})

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


