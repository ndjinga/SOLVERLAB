INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${MEDCOUPLING_INCLUDE_DIR}
  ${MESH_DIR}/inc 
  ${BASE_DIR}/inc
  )

if (CDMATH_WITH_PETSC)
    INCLUDE_DIRECTORIES(
      ${PETSC_INCLUDES}
      ${LINEARSOLVER_DIR}/inc
      )
else (CDMATH_WITH_PETSC)
endif (CDMATH_WITH_PETSC)


SET(MESH_MED
  ../ressources/meshSquare.med
  ../ressources/meshSphere.med
  ../ressources/meshCube.med
	)

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/tests) 

if (CDMATH_WITH_PETSC)
    FILE(GLOB TestMesh_SOURCES . *.cxx)
    ADD_LIBRARY(TestMesh SHARED ${TestMesh_SOURCES})
    TARGET_LINK_LIBRARIES(TestMesh base mesh linearsolver medloader medcoupling med ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
else (CDMATH_WITH_PETSC)  
    FILE(GLOB TestMesh_SOURCES . CellTests.cxx FaceTests.cxx FieldTests.cxx MatrixTests.cxx MeshTests.cxx NodeTests.cxx PointTests.cxx VectorTests.cxx)
    ADD_LIBRARY(TestMesh SHARED ${TestMesh_SOURCES})
    TARGET_LINK_LIBRARIES(TestMesh base mesh medloader medcoupling med)
endif (CDMATH_WITH_PETSC) 

INSTALL(TARGETS TestMesh DESTINATION bin/unit_tests)
