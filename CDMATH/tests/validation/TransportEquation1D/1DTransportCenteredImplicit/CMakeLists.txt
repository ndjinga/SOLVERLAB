
SET(SCRIPT
  ./TransportEquation1DCenteredImplicit.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DTransportEquationCenteredImplicit)

if (CDMATH_WITH_PYTHON )

    SET(CFL  1  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Stiff_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Stiff_CFL1 PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Stiff_CFL1 APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Smooth_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Smooth_CFL1 PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Smooth_CFL1 APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

    SET(CFL  2  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Stiff_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Stiff_CFL2 PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Stiff_CFL2 APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")


    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquation_CenteredImplicit_Smooth_CFL2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationCenteredImplicit.py  ${CFL} ${SMOOTHNESS} )
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Smooth_CFL2 PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validation1DTransportEquation_CenteredImplicit_Smooth_CFL2 APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON )


