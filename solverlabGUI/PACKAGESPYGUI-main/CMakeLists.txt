# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8 FATAL_ERROR)
INCLUDE(CMakeDependentOption)

SET(MODULE_NAME PACKAGESPY CACHE STRING "name of current module upperCase")
string(TOLOWER ${MODULE_NAME} MODULE_NAME_LC)

PROJECT( Salome${MODULE_NAME} C CXX )

# get PACKAGESPY cmake macro etc
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_files" )

include( TraceNewMacros )
setcolors()

# Ensure a proper linker behavior:
CMAKE_POLICY(SET CMP0003 NEW)

# Versioning
# ===========
# Project name, upper case
STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)

SET(${PROJECT_NAME_UC}_MAJOR_VERSION 9)
SET(${PROJECT_NAME_UC}_MINOR_VERSION 9)
SET(${PROJECT_NAME_UC}_PATCH_VERSION 0)
SET(${PROJECT_NAME_UC}_VERSION
  ${${PROJECT_NAME_UC}_MAJOR_VERSION}.${${PROJECT_NAME_UC}_MINOR_VERSION}.${${PROJECT_NAME_UC}_PATCH_VERSION})
SET(${PROJECT_NAME_UC}_VERSION_DEV 0)



# Common CMake macros
SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake configuration files")
IF(EXISTS ${CONFIGURATION_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  INCLUDE(SalomeMacros)
ELSE()
  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
ENDIF()

# Find KERNEL
# ===========
SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
IF(EXISTS ${KERNEL_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${KERNEL_ROOT_DIR}/salome_adm/cmake_files")
  FIND_PACKAGE(SalomeKERNEL REQUIRED)
  KERNEL_WITH_CORBA() # check whether KERNEL builded with CORBA
ELSE(EXISTS ${KERNEL_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR")
ENDIF(EXISTS ${KERNEL_ROOT_DIR})

# User options
# ============
OPTION(SALOME_BUILD_TESTS "Build SALOME tests" ON)
OPTION(SALOME_BUILD_DOC "Generate SALOME PACKAGESPY documentation" ON)
OPTION(CMAKE_TRACE "Generate cmake trace variables in log" OFF)

IF(SALOME_BUILD_TESTS)
  ENABLE_TESTING()
ENDIF()

# KERNEL optionals: to compile pyc pyo
FIND_PACKAGE(SalomePythonInterp REQUIRED)
FIND_PACKAGE(SalomePythonLibs REQUIRED)

IF(SALOME_BUILD_DOC)
  FIND_PACKAGE(SalomeSphinx)
  SALOME_LOG_OPTIONAL_PACKAGE(Sphinx SALOME_BUILD_DOC)
ENDIF()

# Find GUI
# ===========
SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR} CACHE PATH "Path to the Salome GUI")
IF(EXISTS ${GUI_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${GUI_ROOT_DIR}/adm_local/cmake_files")
  FIND_PACKAGE(SalomeGUI REQUIRED)
  SALOME_GUI_WITH_CORBA() #remove this line if you don't need GUI module builded with CORBA
  #remove unnecessary options from this list (see SALOME_GUI_MODE() macro description)
  SALOME_GUI_MODE(SALOME_USE_SALOMEOBJECT SALOME_USE_OCCVIEWER SALOME_USE_VTKVIEWER SALOME_USE_PYCONSOLE
                  OPTIONAL SALOME_USE_PLOT2DVIEWER SALOME_USE_GLVIEWER SALOME_USE_GRAPHICSVIEW
                  SALOME_USE_QXGRAPHVIEWER SALOME_USE_PVVIEWER SALOME_USE_PYVIEWER
                  )
  ADD_DEFINITIONS(${GUI_DEFINITIONS})
  INCLUDE_DIRECTORIES(${GUI_INCLUDE_DIRS})
ELSE(EXISTS ${GUI_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome GUI, please define GUI_ROOT_DIR")
ENDIF(EXISTS ${GUI_ROOT_DIR})

# Qt
FIND_PACKAGE(SalomeQt5 REQUIRED)


# Directories
# (default values taken from KERNEL)
# ===========
SET(SALOME_SHORT_INSTALL_DOC "${CMAKE_INSTALL_PREFIX}/doc" CACHE PATH "Install path: SALOME documentation root new design")
SET(SALOME_SHORT_HTML_STATIC_PATH "${SALOME_SHORT_INSTALL_DOC}/_static" CACHE PATH "Install path: SALOME html static data")

# Sources
# ========
IF(SALOME_BUILD_DOC)
  ADD_SUBDIRECTORY(doc)
ENDIF()

# Additional files
# ================
SET(_salome_context_API_SCRIPTS
   ${KERNEL_ROOT_DIR}/bin/salome/salomeContext.py
   ${KERNEL_ROOT_DIR}/bin/salome/salomeContextUtils.py
   ${KERNEL_ROOT_DIR}/bin/salome/parseConfigFile.py
)
SALOME_INSTALL_SCRIPTS("${_salome_context_API_SCRIPTS}" ${SALOME_INSTALL_SCRIPT_PYTHON})

#SET(_salome_AllTestLauncher_SCRIPTS
#   ${CMAKE_CURRENT_SOURCE_DIR}/AllTestLauncherPACKAGESPY.sh
#)

#SALOME_INSTALL_SCRIPTS("${_salome_AllTestLauncher_SCRIPTS}" ${CMAKE_INSTALL_PREFIX})

# --- rules ---
#SALOME_INSTALL_SCRIPTS("${_python_files}" ${SALOME_INSTALL_PYTHON})

#copy directory as it
INSTALL(
  DIRECTORY pythonAppliMatix
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  USE_SOURCE_PERMISSIONS
  PATTERN ".git" EXCLUDE
  PATTERN "tests_crescendo21" EXCLUDE
)

INSTALL(
  DIRECTORY splash
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  USE_SOURCE_PERMISSIONS
  PATTERN ".git" EXCLUDE
)

# resume TRACE log
# "GUI_." is regExp
printCmakeTraceVariable( "GUI_." )

# optionally get trace of all cmake main variables in log
if( CMAKE_TRACE )
  printCmakeTrace()
endif()
