# Copyright (C) 2007-2019  CEA/DEN, EDF R&D, OPEN CASCADE
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

MESSAGE(STATUS "Check for XDR ...")

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(XDR_INCLUDE_DIRS rpc/xdr.h PATH_SUFFIXES tirpc)
IF(XDR_INCLUDE_DIRS)
  SET(XDR_DEFINITIONS "-DHAS_XDR")
ENDIF(XDR_INCLUDE_DIRS)

IF(WIN32)
  FIND_LIBRARY(XDR_LIBRARIES xdr)                 # To get the .lib file from XDR
  FIND_PATH(XDR_INCLUDE_DIRS2 stdint.h PATH_SUFFIXES src/msvc)  # To get the stdint.h from XDR (needed by types.h)
  IF(XDR_INCLUDE_DIRS)
    IF(XDR_INCLUDE_DIRS2)
      LIST(APPEND XDR_INCLUDE_DIRS "${XDR_INCLUDE_DIRS2}")
    ELSE()
      SET(XDR_INCLUDE_DIRS "${XDR_INCLUDE_DIRS2}")  # Make the detection fail
    ENDIF()
  ENDIF()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(XDR REQUIRED_VARS XDR_INCLUDE_DIRS XDR_LIBRARIES)
ELSE(WIN32)
  FIND_LIBRARY(XDR_LIBRARY NAMES tirpc xdr)
  IF(NOT XDR_LIBRARY)
    MESSAGE(STATUS "Could not find XDR libraries ...")
  ELSE()
    MESSAGE(STATUS "Found XDR libraries ${XDR_LIBRARY} ...")
    SET(XDR_LIBRARIES ${XDR_LIBRARY})
  ENDIF()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(XDR REQUIRED_VARS XDR_INCLUDE_DIRS)
ENDIF(WIN32)
