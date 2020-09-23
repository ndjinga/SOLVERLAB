cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project (CoreFlows C CXX)
set (COREFLOWS_VERSION_MAJOR 0)
set (COREFLOWS_VERSION_MINOR 1)

option (COREFLOWS_WITH_DOCUMENTATION "Generate documentation for CoreFlows use." ON)
option (COREFLOWS_WITH_PYTHON "Compile Python interface for COREFLOWS." ON)
option (COREFLOWS_WITH_GUI "Compile Graphic user Interface for COREFLOWS." OFF)
option (COREFLOWS_WITH_PACKAGE "Generate RPM, Debian and tarball packages." OFF)

#Path to installed libraries
set (PETSC_DIR            OFF CACHE STRING "PETSc library path" )
set (SLEPC_DIR            OFF CACHE STRING "SLEPc library path" )

enable_testing()											   #
# Debug mode activates testing and profiling --------------------------------------------------------------#
if (CMAKE_BUILD_TYPE STREQUAL Debug)									   #
													   #
 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg --coverage")						   #
 SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg --coverage")					   #
 SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg --coverage")				   #
													   #
endif ()												   #
													   #
#----------------------------------------------------------------------------------------------------------#

# Swig interface -------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_PYTHON)                                                                                  #
   set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")     #
   set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python-devel, numpy, swig")              #
   set (CoreFlows_SWIG_DIR ${CoreFlows_SOURCE_DIR}/swig)                                                    #
endif ()                                                                                                    #
                                                                                                            #
#-----------------------------------------------------------------------------------------------------------#

# Base directories
set (CoreFlows_SRC ${CoreFlows_SOURCE_DIR}/Models ) 
set (CoreFlows_EXAMPLES ${CoreFlows_SOURCE_DIR}/examples)

set( CoreFlows_INCLUDES 
  ${CDMATH_INCLUDES}											    #
  ${PETSC_INCLUDES}											    #
  ${CoreFlows_SRC}/inc											    #    
  )													    #

add_subdirectory (${CoreFlows_SRC})
add_subdirectory (${CoreFlows_EXAMPLES})
if (COREFLOWS_WITH_PYTHON)                                                                                  #
   add_subdirectory (${CoreFlows_SWIG_DIR})                                                                 #
endif ()                                                                                                    #

# Documentation --------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_DOCUMENTATION)                                                                           #
  configure_file (Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)                              #
  FILE(GLOB_RECURSE CoreFlows_Documentation ${CoreFlows_SOURCE_DIR}/Documentation *.dox)                    #
  SET(Documentation_SRC                                                                                     #
    ${PROJECT_BINARY_DIR}/Doxyfile                                                                          #
    ${CoreFlows_Documentation}                                                                              #
      )                                                                                                     #
  add_custom_target (docCoreFlows COMMAND ${DOXYGEN_EXECUTABLE}                                                  #
  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)			                                            #
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/coreflows-dev-doc OPTIONAL)                  #
  INSTALL(FILES ${CoreFlows_SOURCE_DIR}/Documentation/CoreFlows.pdf DESTINATION share/doc OPTIONAL)	    #
endif ()                                                                                                    #
                                                                                                            #
# ------------- --------------------------------------------------------------------------------------------#

# Graphic interface-----------------------------------------------------------------------------------------#
     													    #
if (COREFLOWS_WITH_GUI)											    					  														    #
 set (CoreFlows_GUI_DIR ${CoreFlows_SOURCE_DIR}/gui)                                                        #
 add_subdirectory(${CoreFlows_GUI_DIR})			

endif ()                                                                                                    #
													    #
# ------------- --------------------------------------------------------------------------------------------#

#--------------------- COMPILATION MAIN --------------------------------------------------------------------#
													    #
INCLUDE_DIRECTORIES(						   					    #
  ${CoreFlows_INCLUDES}											    #
  ${PETSC_INCLUDES}
  ${SLEPC_INCLUDES}
  )													    #
													    #
SET(CoreFlowsMain_SOURCES										    #
    ${CoreFlows_SRC}/src/Fluide.cxx									    #
    ${CoreFlows_SRC}/src/DiffusionEquation.cxx								    #
    ${CoreFlows_SRC}/src/StationaryDiffusionEquation.cx							    #
    ${CoreFlows_SRC}/src/ProblemFluid.cxx								    #
    ${CoreFlows_SRC}/src/IsothermalTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/utilitaire_algebre.cxx								    #
    ${CoreFlows_SRC}/src/DriftModel.cxx									    #
    ${CoreFlows_SRC}/src/SinglePhase.cxx								    #
    ${CoreFlows_SRC}/src/ProblemCoreFlows.cxx								    #
    ${CoreFlows_SRC}/src/FiveEqsTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/TransportEquation.cxx								    #
  )													    #
													    #
SET(MAIN_SOURCES											    #
  mainCoreFlows.cxx												    #
  )													    #
SET(_extra_lib_CoreFlows ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} )					    #
ADD_EXECUTABLE(CoreFlowsMainExe ${MAIN_SOURCES})							    #
SET_TARGET_PROPERTIES(CoreFlowsMainExe PROPERTIES COMPILE_FLAGS "")					    #
TARGET_LINK_LIBRARIES(CoreFlowsMainExe CoreFlows  ${_extra_lib_CoreFlows})				    #
INSTALL(TARGETS CoreFlowsMainExe DESTINATION bin)						    #
													    #
configure_file(
    ${PROJECT_SOURCE_DIR}/env_CoreFlows.sh
    @ONLY
)
													    
													    
# ------------- --------------------------------------------------------------------------------------------#


# Packaging
if (COREFLOWS_WITH_PACKAGE)
  include (${CoreFlows_SOURCE_DIR}/CPackLists.txt)
endif ()
