# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8 FATAL_ERROR)
INCLUDE(CMakeDependentOption)

SET(MODULE_NAME SOLVERLABGUI CACHE STRING "name of current module upperCase")
string(TOLOWER ${MODULE_NAME} MODULE_NAME_LC)

PROJECT(Salome${MODULE_NAME} C CXX)

# Ensure a proper linker behavior:
CMAKE_POLICY(SET CMP0003 NEW)

# Versioning
# ===========
# Project name, upper case
STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)

SET(${PROJECT_NAME_UC}_MAJOR_VERSION  9)
SET(${PROJECT_NAME_UC}_MINOR_VERSION 10)
SET(${PROJECT_NAME_UC}_PATCH_VERSION  0)
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
OPTION(SALOME_BUILD_DOC "Generate SALOME SOLVERLABGUI documentation" ON)
OPTION(CMAKE_TRACE "Generate cmake trace variables in log" ON)

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

##
## GUI specifics
##

# OCCT
FIND_PACKAGE(SalomeCAS REQUIRED)

# Qt
FIND_PACKAGE(SalomeQt5 REQUIRED)


# Directories
# (default values taken from KERNEL)
# ===========
SET(SALOME_INSTALL_BINS "${SALOME_INSTALL_BINS}" CACHE PATH "Install path: SALOME binaries")
SET(SALOME_INSTALL_LIBS "${SALOME_INSTALL_LIBS}" CACHE PATH "Install path: SALOME libs")
SET(SALOME_INSTALL_HEADERS "${SALOME_INSTALL_HEADERS}" CACHE PATH "Install path: SALOME headers")
SET(SALOME_INSTALL_SCRIPT_SCRIPTS "${SALOME_INSTALL_SCRIPT_SCRIPTS}" CACHE PATH
   "Install path: SALOME scripts")
SET(SALOME_INSTALL_SCRIPT_DATA "${SALOME_INSTALL_SCRIPT_DATA}" CACHE PATH
   "Install path: SALOME script data")
SET(SALOME_INSTALL_SCRIPT_PYTHON "${SALOME_INSTALL_SCRIPT_PYTHON}" CACHE PATH
   "Install path: SALOME Python scripts")
SET(SALOME_INSTALL_PYTHON "${SALOME_INSTALL_PYTHON}" CACHE PATH "Install path: SALOME Python stuff")
SET(SALOME_INSTALL_PYTHON_SHARED "${SALOME_INSTALL_PYTHON_SHARED}" CACHE PATH
   "Install path: SALOME Python shared modules")
SET(SALOME_INSTALL_CMAKE "${SALOME_INSTALL_CMAKE}" CACHE PATH "Install path: SALOME CMake files")
SET(SALOME_INSTALL_CMAKE_LOCAL "${SALOME_INSTALL_CMAKE_LOCAL}" CACHE PATH
    "Install path: local SALOME CMake files")
SET(SALOME_INSTALL_AMCONFIG_LOCAL "${SALOME_INSTALL_AMCONFIG_LOCAL}" CACHE PATH
  "Install path: local SALOME config files (obsolete, to be removed)")

SET(SALOME_INSTALL_RES "${SALOME_INSTALL_RES}" CACHE PATH "Install path: SALOME resources")
SET(SALOME_INSTALL_DOC "${SALOME_INSTALL_DOC}" CACHE PATH "Install path: SALOME documentation share old design")

## SOLVERLABGUI specific
SET(SALOME_SOLVERLABGUI_INSTALL_RES "${SALOME_INSTALL_RES}/solverlabgui" CACHE PATH "Install path: SOLVERLABGUI resources")
SET(SALOME_SHORT_INSTALL_DOC "${CMAKE_INSTALL_PREFIX}/doc" CACHE PATH "Install path: SALOME documentation root new design")
SET(SALOME_SHORT_HTML_STATIC_PATH "${SALOME_SHORT_INSTALL_DOC}/_static" CACHE PATH "Install path: SALOME html static data")

MARK_AS_ADVANCED(SALOME_INSTALL_BINS SALOME_INSTALL_LIBS SALOME_INSTALL_IDLS SALOME_INSTALL_HEADERS)
MARK_AS_ADVANCED(SALOME_INSTALL_SCRIPT_SCRIPTS SALOME_INSTALL_SCRIPT_DATA SALOME_INSTALL_SCRIPT_PYTHON)
MARK_AS_ADVANCED(SALOME_INSTALL_APPLISKEL_SCRIPTS  SALOME_INSTALL_APPLISKEL_PYTHON SALOME_INSTALL_CMAKE SALOME_INSTALL_CMAKE_LOCAL SALOME_INSTALL_RES)
MARK_AS_ADVANCED(SALOME_INSTALL_PYTHON SALOME_INSTALL_PYTHON_SHARED)
MARK_AS_ADVANCED(SALOME_GUI_INSTALL_PARAVIEW_LIBS SALOME_GUI_INSTALL_RES_DATA SALOME_GUI_INSTALL_RES_SCRIPTS SALOME_GUI_INSTALL_PLUGINS)
MARK_AS_ADVANCED(SALOME_SOLVERLABGUI_INSTALL_RES SALOME_INSTALL_DOC)


# Sources
# ========
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(bin)
ADD_SUBDIRECTORY(resources)
IF(SALOME_BUILD_DOC)
  ADD_SUBDIRECTORY(doc)
ENDIF()

# Header configuration
# ====================
SALOME_XVERSION(${PROJECT_NAME})
SALOME_CONFIGURE_FILE(SOLVERLABGUI_version.h.in SOLVERLABGUI_version.h INSTALL ${SALOME_INSTALL_HEADERS})

# Configuration export
# (here only the level 1 prerequisites are exposed)
# ====================
#   (not done for THIS MODULE)

# --- rules ---
SALOME_INSTALL_SCRIPTS("${_python_files}" ${SALOME_INSTALL_PYTHON})

#copy files and rename in bin (as salome matix)
INSTALL(
  PROGRAMS AllTestLauncher.sh
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  RENAME AllSolverlabGuiTestLauncher.sh
)

INSTALL(
  PROGRAMS solverlabGUI
  DESTINATION ${SALOME_INSTALL_BINS}
)

#copy directory ./examples as it
INSTALL(
  DIRECTORY example
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  USE_SOURCE_PERMISSIONS
  PATTERN ".git" EXCLUDE
)

# avoid duplicates with SOLVERLABGUI and PACKAGESPY
# configparserpy debogpy filewatcherpy
# interpreterpy matplotlibpy
# returncodepy ribbonpy
# salomepy settingspy test unittestpy
# widgetpy xyzpy

INSTALL(
  DIRECTORY solverlabpy
            resources
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  USE_SOURCE_PERMISSIONS
  PATTERN ".git" EXCLUDE
)

# to get trace of all cmake main variables in sat log
IF(CMAKE_TRACE)
  message("%% CMAKE_TRACE BEGIN")
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message("%% ${_variableName}=${${_variableName}}")
  endforeach()
  message("%% CMAKE_TRACE END")
ENDIF()
