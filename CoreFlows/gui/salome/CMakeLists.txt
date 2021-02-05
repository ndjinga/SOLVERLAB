# Copyright (C) 2012-2014  CEA/DEN, EDF R&D
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
INCLUDE(CMakeDependentOption)

# Versioning
# ===========
# Project name, upper case
STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)

SET(${PROJECT_NAME_UC}_MAJOR_VERSION 7)
SET(${PROJECT_NAME_UC}_MINOR_VERSION 8)
SET(${PROJECT_NAME_UC}_PATCH_VERSION 0)
SET(${PROJECT_NAME_UC}_VERSION
  ${${PROJECT_NAME_UC}_MAJOR_VERSION}.${${PROJECT_NAME_UC}_MINOR_VERSION}.${${PROJECT_NAME_UC}_PATCH_VERSION})
SET(${PROJECT_NAME_UC}_VERSION_DEV 1)

# User options
# ============
OPTION(SALOME_BUILD_DOC "Generate SALOME CoreFlows documentation" ON)
OPTION(SALOME_BUILD_TESTS "Generate SALOME CoreFlows tests" ON)

# Common CMake macros
# ===================
SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake configuration files")
IF(EXISTS ${CONFIGURATION_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  INCLUDE(SalomeMacros NO_POLICY_SCOPE)
ELSE()
  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
ENDIF()

# Find KERNEL 
# ==============
SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
IF( EXISTS ${KERNEL_ROOT_DIR} )
  FIND_PACKAGE(SalomeKERNEL REQUIRED)
ELSE( EXISTS ${KERNEL_ROOT_DIR} )
  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR!")
ENDIF( EXISTS ${KERNEL_ROOT_DIR} )

# Find SALOME GUI 
# ==============
SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR} CACHE PATH "Path to the Salome GUI")
IF(EXISTS ${GUI_ROOT_DIR})
  FIND_PACKAGE(SalomeGUI)
ELSE(EXISTS ${GUI_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome GUI, please define GUI_ROOT_DIR!")
ENDIF(EXISTS ${GUI_ROOT_DIR})

# Platform setup
# ==============
INCLUDE(SalomeSetupPlatform)   # From KERNEL
    
# Prerequisites
# =============

# Mandatory products
IF(  ${SalomeKERNEL_VERSION} STRLESS "7.8.0")
         FIND_PACKAGE(SalomePython REQUIRED)
ELSE()
         FIND_PACKAGE(SalomePythonInterp REQUIRED)
         FIND_PACKAGE(SalomePythonLibs REQUIRED)
ENDIF()
MESSAGE(STATUS "SalomeKERNEL_VERSION is ${SalomeKERNEL_VERSION}")

# Qt5
FIND_PACKAGE(SalomeQt5 REQUIRED COMPONENTS QtCore QtGui)

IF(SALOME_BUILD_DOC)
#  FIND_PACKAGE(SalomeDoxygen)
#  FIND_PACKAGE(SalomeSphinx)
#  SALOME_LOG_OPTIONAL_PACKAGE(Doxygen SALOME_BUILD_DOC)
#  SALOME_LOG_OPTIONAL_PACKAGE(Sphinx SALOME_BUILD_DOC)
#  ADD_DEFINITIONS(-DDOXYGEN_IS_OK)
ENDIF()

IF(SALOME_BUILD_TESTS)
  ENABLE_TESTING()
ENDIF()

# Detection summary:
SALOME_PACKAGE_REPORT_AND_CHECK()

# Directories
# ===========
SET(SALOME_INSTALL_LIBS "${SALOME_INSTALL_LIBS}" CACHE PATH "Install path: SALOME libs")
SET(SALOME_INSTALL_PYTHON "${SALOME_INSTALL_PYTHON}" CACHE PATH 
    "Install path: SALOME Python scripts")
SET(SALOME_INSTALL_SCRIPT_PYTHON "${SALOME_INSTALL_SCRIPT_PYTHON}" CACHE PATH 
    "Install path: SALOME Python scripts")
SET(SALOME_INSTALL_CMAKE_LOCAL ${SALOME_INSTALL_CMAKE_LOCAL} CACHE PATH 
    "Install path: local SALOME CMake files") 
SET(SALOME_INSTALL_RES "${SALOME_INSTALL_RES}" CACHE PATH "Install path: SALOME resources")
SET(SALOME_INSTALL_DOC "${SALOME_INSTALL_DOC}" CACHE PATH "Install path: SALOME documentation")
SET(SALOME_INSTALL_HEADERS ${SALOME_INSTALL_HEADERS} CACHE PATH "Install path: SALOME headers")
 
# Specific to CoreFlows:   
SET(SALOME_CoreFlows_INSTALL_RES_DATA ${SALOME_INSTALL_RES}/solverlab CACHE PATH "Install path: SALOME SOLVERLAB specific resources")

# Sources 
# ========
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(resources)
IF(SALOME_BUILD_DOC)
#   ADD_SUBDIRECTORY(doc)
ENDIF()

