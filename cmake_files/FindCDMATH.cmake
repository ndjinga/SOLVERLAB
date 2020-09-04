##Copyright (C) arthurtalpaert.
##All rights reserved.
##
##Redistribution and use in source and binary forms, with or without modification,
##are permitted provided that the following conditions are met:
##
##* Redistributions of source code must retain the above copyright notice, this
##  list of conditions and the following disclaimer.
##
##* Redistributions in binary form must reproduce the above copyright notice, this
##  list of conditions and the following disclaimer in the documentation and/or
##  other materials provided with the distribution.
##
##THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
##ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
##WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
##ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
##(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
##LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
##ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
##SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# - Try to find CDMATH
# Once done this will define
#
#  CDMATH_FOUND        - system has CDMATH
#  CDMATH_INCLUDES     - the CDMATH include directories
#  CDMATH_LIBRARIES    - Link these to use CDMATH
#
#  Usage:
#  find_package(CDMATH)
#
# Setting these changes the behavior of the search:
#  CDMATH_DIR - directory in which CDMATH resides
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

find_path (CDMATH_DIR include/CdmathException.hxx
  HINTS ENV CDMATH_DIR
  PATHS
  /usr
  $ENV{HOME}/cdmath
  $ENV{HOME}/workspace/cdmath_install
  DOC "CDMATH Directory")
message (STATUS "Found CDMATH: ${CDMATH_DIR}")

# Include directories
set(MED_INCLUDES $ENV{MEDFILE_INCLUDE_DIRS})
if (NOT (IS_DIRECTORY  ${MED_INCLUDES}) )
  message (SEND_ERROR "MED_INCLUDES can not be used, ${MED_INCLUDES} does not exist.")
endif () 
set(MEDCOUPLING_INCLUDES $ENV{MEDCOUPLING_INCLUDE_DIR})
if (NOT (IS_DIRECTORY  ${MEDCOUPLING_INCLUDES}) )
  message (SEND_ERROR "MEDCOUPLING_INCLUDES can not be used, ${MEDCOUPLING_INCLUDES} does not exist.")
endif () 
# This sets the variable ${CDMATH_INCLUDES}.
set(CDMATH_INCLUDES ${CDMATH_DIR}/include ${MED_INCLUDES} ${MEDCOUPLING_INCLUDES} )
if (NOT (IS_DIRECTORY  ${CDMATH_DIR}/include) )
  message (SEND_ERROR "CDMATH_INCLUDES can not be used, ${CDMATH_DIR}/include does not exist.")
endif () 

# CDMATH libraries against which to link
# This sets the variable ${CDMATH_LIBRARIES}.
set(CDMATH_LIBDIR ${CDMATH_DIR}/lib)
if ( NOT (IS_DIRECTORY  ${CDMATH_LIBDIR}) )
  message (SEND_ERROR "CDMATH_LIBDIR can not be used, ${CDMATH_LIBDIR} does not exist.")
endif () 
find_library (CDMATHBASE_LIB NAMES base PATHS ${CDMATH_LIBDIR})
find_library (CDMATHMESH_LIB NAMES mesh PATHS ${CDMATH_LIBDIR})
find_library (CDMATHLINEARSOLVER_LIB NAMES linearsolver PATHS ${CDMATH_LIBDIR})
find_library (MEDC_LIB NAMES medC PATHS $ENV{MEDFILE_ROOT_DIR}/lib)
find_library (MEDLOADER_LIB NAMES medloader PATHS $ENV{MEDCOUPLING_LIBRARIES})
find_library (MEDCOUPLING_LIB NAMES medcoupling PATHS $ENV{MEDCOUPLING_LIBRARIES})
set (CDMATH_LIBRARIES
	${MEDC_LIB} 
	${MEDLOADER_LIB} 
	${MEDCOUPLING_LIB}
	${CDMATHBASE_LIB} 
	${CDMATHMESH_LIB} 
	${CDMATHLINEARSOLVER_LIB}
	)
