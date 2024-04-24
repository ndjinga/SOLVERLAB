
if (COREFLOWS_WITH_PYTHON )
    IF(EXISTS ${KERNEL_ROOT_DIR})
        ADD_TEST(ExampleGEOMSMESH_2D_SinusoidalBump_triangle ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SinusoidalBump/mesh_salome_triangle.py)

        ADD_TEST(ExampleGEOMSMESH_2D_SinusoidalBump_hexaedre ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SinusoidalBump/mesh_salome_hexaedre.py)

        ADD_TEST(ExampleGEOMSMESH_2D_BackwardFacingStep_triangle ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/BackwardFacingStep/backward_facing_step_triangles.py)

        ADD_TEST(ExampleGEOMSMESH_2D_BackwardFacingStep_hexaedre ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/BackwardFacingStep/backward_facing_step_hexaedra.py)

    ENDIF(EXISTS ${KERNEL_ROOT_DIR})
endif (COREFLOWS_WITH_PYTHON )


