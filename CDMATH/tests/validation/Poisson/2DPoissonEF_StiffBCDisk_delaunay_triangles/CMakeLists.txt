
SET(SCRIPT
  ../../scripts/Poisson2DEF_StiffBC_DISK/FiniteElements2DPoissonStiffBC_DISK.py
  test_validation2DPoissonEF_StiffBCDisk_Delaunay_triangles.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  "../../../ressources/2DdiskWithTriangles/diskWithTriangles_*" )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_TEST(validationPoisson_2DFE_StiffBCDisk_Delaunay_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DPoissonEF_StiffBCDisk_Delaunay_triangles.py)

    SET_PROPERTY(TEST validationPoisson_2DFE_StiffBCDisk_Delaunay_triangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")

    SET_PROPERTY(TEST validationPoisson_2DFE_StiffBCDisk_Delaunay_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


