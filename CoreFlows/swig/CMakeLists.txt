INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(CoreFlows.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(CoreFlows.i PROPERTIES C OFF)
SET_SOURCE_FILES_PROPERTIES(CoreFlows.i PROPERTIES SWIG_DEFINITIONS "-shadow")

INCLUDE_DIRECTORIES(
  ${Python_INCLUDE_DIRS}
  ${CoreFlows_INCLUDES}									    #
  ${PETSC_INCLUDES}
  ${SLEPC_INCLUDES}
  ${SOLVERLAB_SOURCE_DIR}/CDMATH/mesh/inc # To include Mesh.hxx for the enum EntityType
)

IF(    ${SOLVERLAB_WITH_MPI} )
  if(EXISTS ${MPI4PY_ROOT_DIR}/mpi4py/include/mpi4py/mpi4py.i)
    INCLUDE_DIRECTORIES( ${MPI4PY_ROOT_DIR}/mpi4py/include
                         ${MPI4PY_ROOT_DIR}/mpi4py/include/mpi4py )
  elseif (EXISTS ${MPI4PY_ROOT_DIR}/include/mpi4py/mpi4py.i )
    INCLUDE_DIRECTORIES( ${MPI4PY_ROOT_DIR}/include
                         ${MPI4PY_ROOT_DIR}/include/mpi4py )
  endif()
ENDIF( ${SOLVERLAB_WITH_MPI} )

if(     EXISTS ${PETSC4PY_ROOT_DIR}/petsc4py/include/petsc4py/petsc4py.i)
  INCLUDE_DIRECTORIES( ${PETSC4PY_ROOT_DIR}/petsc4py/include
                       ${PETSC4PY_ROOT_DIR}/petsc4py/include/petsc4py )
elseif (EXISTS ${PETSC4PY_ROOT_DIR}/include/petsc4py/petsc4py.i )
  INCLUDE_DIRECTORIES( ${PETSC4PY_ROOT_DIR}/include
                       ${PETSC4PY_ROOT_DIR}/include/petsc4py )
endif(  EXISTS ${PETSC4PY_ROOT_DIR}/petsc4py/include/petsc4py/petsc4py.i)

if(     EXISTS ${SLEPC4PY_ROOT_DIR}/slepc4py/include/slepc4py/slepc4py.i)
  INCLUDE_DIRECTORIES( ${SLEPC4PY_ROOT_DIR}/slepc4py/include
                       ${SLEPC4PY_ROOT_DIR}/slepc4py/include/slepc4py )
elseif (EXISTS ${SLEPC4PY_ROOT_DIR}/include/slepc4py/slepc4py.i )
  INCLUDE_DIRECTORIES( ${SLEPC4PY_ROOT_DIR}/include
                       ${SLEPC4PY_ROOT_DIR}/include/slepc4py )
endif(  EXISTS ${SLEPC4PY_ROOT_DIR}/slepc4py/include/slepc4py/slepc4py.i)

IF   (${CMAKE_VERSION} VERSION_LESS "3.8.0") 
  SWIG_ADD_MODULE(CoreFlows python CoreFlows.i)
ELSE (${CMAKE_VERSION} VERSION_LESS "3.8.0") 
  SWIG_ADD_LIBRARY(CoreFlows TYPE SHARED LANGUAGE python SOURCES CoreFlows.i)
ENDIF(${CMAKE_VERSION} VERSION_LESS "3.8.0")

SWIG_LINK_LIBRARIES(CoreFlows ${Python_LIBRARIES} ${CDMATH_LIBRARIES} CoreFlowsLibs ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

INSTALL(TARGETS _CoreFlows DESTINATION lib/coreflows)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoreFlows.py DESTINATION bin/coreflows)
INSTALL(FILES ${CoreFlows_SOURCE_DIR}/mainCoreFlows.py DESTINATION bin/)

if(    SOLVERLAB_WITH_PYTHON )
  INCLUDE_DIRECTORIES(
    ${Python_INCLUDE_DIRS}
    ${PROJECT_BINARY_DIR}/../CDMATH/swig											    #    
  )

  file(COPY solverlab.i DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/solverlab.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/solverlab.i PROPERTIES C OFF)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/solverlab.i PROPERTIES SWIG_DEFINITIONS "-shadow")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/solverlab.i PROPERTIES DEPENDS _cdmath)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/solverlab.i PROPERTIES DEPENDS _CoreFlows)

  SWIG_ADD_LIBRARY(solverlab TYPE SHARED LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/solverlab.i )
  SWIG_LINK_LIBRARIES(solverlab ${Python_LIBRARIES} ${CDMATH_LIBRARIES} CoreFlowsLibs ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})

  INSTALL(TARGETS _solverlab DESTINATION lib)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/solverlab.py DESTINATION bin)
endif( SOLVERLAB_WITH_PYTHON )


