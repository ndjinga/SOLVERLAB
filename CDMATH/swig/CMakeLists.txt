INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(
  ${Python_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${MEDCOUPLING_INCLUDE_DIR}
  ${BASE_DIR}/inc
  ${MESH_DIR}/inc
)

if (CDMATH_WITH_PETSC)
    INCLUDE_DIRECTORIES(
      ${PETSC_INCLUDES}
      ${SLEPC_INCLUDES}
      ${LINEARSOLVER_DIR}/inc
    )
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

    file(COPY cdmath_with_petsc.i DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cdmath_with_petsc.i ${CMAKE_CURRENT_BINARY_DIR}/cdmath.i)
else (CDMATH_WITH_PETSC)
    file(COPY cdmath_without_petsc.i DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cdmath_without_petsc.i ${CMAKE_CURRENT_BINARY_DIR}/cdmath.i)
endif (CDMATH_WITH_PETSC)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/cdmath.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/cdmath.i PROPERTIES C OFF)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/cdmath.i PROPERTIES SWIG_DEFINITIONS "-shadow")

#Trying to remove warnings due to cmake policies cmp0078 and cmp0086
IF   (${CMAKE_VERSION} VERSION_GREATER "3.13.0")
  cmake_policy(SET CMP0078 NEW)
  IF   (${CMAKE_VERSION} VERSION_GREATER "3.14.0")
    cmake_policy(SET CMP0086 NEW)
  ENDIF(${CMAKE_VERSION} VERSION_GREATER "3.14.0")
ENDIF(${CMAKE_VERSION} VERSION_GREATER "3.13.0")

IF   (${CMAKE_VERSION} VERSION_LESS "3.8.0") 
  SWIG_ADD_MODULE(cdmath python ${CMAKE_CURRENT_BINARY_DIR}/cdmath.i)
ELSE (${CMAKE_VERSION} VERSION_LESS "3.8.0") 
  SWIG_ADD_LIBRARY(cdmath TYPE SHARED LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/cdmath.i)
ENDIF(${CMAKE_VERSION} VERSION_LESS "3.8.0")

if (CDMATH_WITH_PETSC) 
  SWIG_LINK_LIBRARIES(cdmath ${Python_LIBRARIES} base mesh ${medloader} ${medcoupling} linearsolver ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
else (CDMATH_WITH_PETSC)
  SWIG_LINK_LIBRARIES(cdmath ${Python_LIBRARIES} base mesh ${medloader} ${medcoupling} )
endif (CDMATH_WITH_PETSC)

INSTALL(TARGETS _cdmath DESTINATION lib/cdmath)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/cdmath.py DESTINATION bin/cdmath)

