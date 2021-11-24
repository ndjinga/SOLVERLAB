cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project (SOLVERLAB C CXX)
set (SOLVERLAB_VERSION_MAJOR 9)
set (SOLVERLAB_VERSION_MINOR 7)
set (SOLVERLAB_VERSION_SUBMINOR alpha)

option (SOLVERLAB_WITH_DOCUMENTATION "Generate SOLVERLAB documentation" ON)
option (SOLVERLAB_WITH_PYTHON        "Generate Python interface of SOLVERLAB using SWIG" ON)
option (SOLVERLAB_WITH_GUI           "Compile a Graphic user Interface for SOLVERLAB" OFF)
option (SOLVERLAB_WITH_PACKAGE       "Generate RPM, Debian and tarball packages" OFF)
option (SOLVERLAB_WITH_COREFLOWS     "Generate COREFLOWS module" ON)
option (SOLVERLAB_WITH_TESTS         "Generate SOLVERLAB example tests" ON)
option (SOLVERLAB_WITH_MPI         "Generate a parallel version of SOLVERLAB" OFF)

#Path to installed libraries
set (PETSC_DIR            OFF CACHE STRING "PETSc library path" )
set (SLEPC_DIR            OFF CACHE STRING "SLEPc library path" )
set (HDF5_ROOT            OFF CACHE STRING "HDF5 library path" )#For the compilation of MED with a specific version of HDF5
set (MEDFILE_ROOT_DIR     OFF CACHE STRING "MED library path" )
set (MEDCOUPLING_ROOT_DIR OFF CACHE STRING "MEDCoupling library path" )
set (PYTHON_ROOT_DIR      OFF CACHE STRING "Python library path" )
set (PARAVIEW_ROOT_DIR    OFF CACHE STRING "PARAVIEW library path" )
set (MPI4PY_ROOT_DIR      OFF CACHE STRING "MPI4PY library path" )

#url of tarball librairies
set (DOWNLOAD_PETSC         http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.16.0.tar.gz 
     CACHE STRING           "PETSc tarball path/url" )
set (DOWNLOAD_SLEPC         https://slepc.upv.es/download/distrib/slepc-3.16.0.tar.gz
     CACHE STRING           "SLEPC tarball path/url" )
set (DOWNLOAD_F2CBLASLAPACK https://www.mcs.anl.gov/petsc/mirror/externalpackages/f2cblaslapack-3.4.2.q4.tar.gz
     CACHE STRING           "F2CBLASLAPACK tarball path/url" )
set (DOWNLOAD_HDF5          https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.3/src/hdf5-1.10.3.tar.gz  
     CACHE STRING           "HDF5 tarball path/url" )
set (DOWNLOAD_MED           http://files.salome-platform.org/Salome/other/med-4.1.0.tar.gz
     CACHE STRING           "MED tarball path/url")
set (DOWNLOAD_MEDCOUPLING   http://files.salome-platform.org/Salome/other/medCoupling-9.7.0.tar.gz
     CACHE STRING           "MEDCoupling tarball path/url" )
set (DOWNLOAD_OPENMPI       https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz
     CACHE STRING           "OPENMPI tarball path/url" )
set (DOWNLOAD_HYPRE         https://github.com/hypre-space/hypre/archive/v2.23.0.zip
     CACHE STRING           "HYPRE tarball path/url" )
set (DOWNLOAD_MPI4PY         https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.1.2.tar.gz
     CACHE STRING           "MPY4PY tarball path/url" )

set (CMAKE_BUILD_TYPE Release CACHE STRING "Installation mode")

include(ExternalProject)#For PETSc, MED and MEDCoupling

list (APPEND CMAKE_MODULE_PATH "${SOLVERLAB_SOURCE_DIR}/cmake_files")


#####################
# Detection of MPI  #
#####################
if   ( DEFINED ENV{MPI_ROOT_DIR} OR DEFINED MPI_ROOT_DIR )
  if   ( DEFINED  ENV{MPI_ROOT_DIR} )
    set(MPI_HOME $ENV{MPI_ROOT_DIR})
  else ( DEFINED  ENV{MPI_ROOT_DIR} )
    set(MPI_HOME    ${MPI_ROOT_DIR})
  endif( DEFINED  ENV{MPI_ROOT_DIR} )
endif( DEFINED ENV{MPI_ROOT_DIR} OR DEFINED MPI_ROOT_DIR )

# Definition of SOLVERLAB_WITH_MPI
if(    DEFINED ENV{MPI_HOME} OR MPI_HOME )
  set(SOLVERLAB_WITH_MPI ON)
  find_package(MPI)
  if(   ${MPI_LIBRARY} STREQUAL "MPI_LIBRARY-NOTFOUND")
    set ( MPI_LIBRARY ${MPI_HOME}/lib/libmpi.so )
  endif(${MPI_LIBRARY} STREQUAL "MPI_LIBRARY-NOTFOUND")
  message ( STATUS "MPI executable is ${MPIEXEC}" )
else(  DEFINED ENV{MPI_HOME} OR MPI_HOME )
  set(SOLVERLAB_WITH_MPI OFF)
endif( DEFINED ENV{MPI_HOME} OR MPI_HOME )

# Definition of MPI4PY_ROOT_DIR
if(   ${SOLVERLAB_WITH_PYTHON} AND ${SOLVERLAB_WITH_MPI} AND NOT ${MPI4PY_ROOT_DIR})
  message(STATUS "Defining MPI4PY directory" )
  if(    DEFINED ENV{MPI4PY_ROOT_DIR})
    set(MPI4PY_ROOT_DIR $ENV{MPI4PY_ROOT_DIR})
  else(  DEFINED ENV{MPI4PY_ROOT_DIR})
    set(MPI4PY_ROOT_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib/mpi4py)
  endif( DEFINED ENV{MPI4PY_ROOT_DIR})
endif(${SOLVERLAB_WITH_PYTHON} AND ${SOLVERLAB_WITH_MPI} AND NOT ${MPI4PY_ROOT_DIR})

if(    MPI4PY_ROOT_DIR)
  message(STATUS "MPI4PY directory is ${MPI4PY_ROOT_DIR}" )
  SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DMPI4PY_ROOT_DIR=${MPI4PY_ROOT_DIR})
endif( MPI4PY_ROOT_DIR)

############################
# Definition of HDF5_ROOT  #
############################
if   ( DEFINED ENV{HDF5_ROOT_DIR} OR DEFINED HDF5_ROOT_DIR )
  if   ( DEFINED  ENV{HDF5_ROOT_DIR} )
    set(HDF5_ROOT $ENV{HDF5_ROOT_DIR})
  else ( DEFINED  ENV{HDF5_ROOT_DIR} )
    set(HDF5_ROOT    ${HDF5_ROOT_DIR})
  endif( DEFINED  ENV{HDF5_ROOT_DIR} )
endif( DEFINED ENV{HDF5_ROOT_DIR} OR DEFINED HDF5_ROOT_DIR )

######################################################
# Detection or compilation of PETSc+SLEPc+HDF5       #
######################################################
  if ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR)
    if   (NOT PETSC_DIR)
      set(PETSC_DIR $ENV{PETSC_DIR})
      message ( STATUS "Checking variable PETSC_DIR : " $ENV{PETSC_DIR} )
    else (NOT PETSC_DIR)
      message ( STATUS "Checking variable PETSC_DIR : " ${PETSC_DIR} )
    endif(NOT PETSC_DIR)

    find_package (PETSc 3.4 REQUIRED)
    petsc_get_version ()

    message ( STATUS "PETSc found. Version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" )
    message ( STATUS "PETSc compiler is ${PETSC_COMPILER}" )
    message ( STATUS "PETSc executable is ${PETSC_MPIEXEC}")

    set(PETSC_INSTALL ${PETSC_DIR})

    #Define and search slepc variables
    if   ( NOT SLEPC_DIR )
      if   ( DEFINED ENV{SLEPC_DIR} )
        set(SLEPC_DIR $ENV{SLEPC_DIR})
      else ( DEFINED ENV{SLEPC_DIR} )
        set(SLEPC_DIR ${PETSC_DIR}/${PETSC_ARCH})
      endif( DEFINED ENV{SLEPC_DIR} )
    endif( NOT SLEPC_DIR)

    message ( STATUS "Checking variable SLEPC_DIR" )
    if ( IS_DIRECTORY ${SLEPC_DIR}/include AND EXISTS ${SLEPC_DIR}/lib/libslepc.so)
      set(SLEPC_INCLUDES  ${SLEPC_DIR}/include)
      set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
      message( STATUS "SLEPc found at ${SLEPC_DIR}" )
    else()
      message( FATAL_ERROR "SLEPc not found at ${SLEPC_DIR}" )
    endif()

    #define hdf5 variables
    if   ( NOT HDF5_ROOT )#hdf5 is not defined in cmake arguments
      if  ( DEFINED ENV{HDF5_ROOT} )#hdf5 is defined in the environment
        set(HDF5_ROOT $ENV{HDF5_ROOT})
      else( DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment
        #HDF5 to be found in petsc external packages
        set(HDF5_ROOT ${PETSC_DIR}/${PETSC_ARCH})# define hint for hdf5/med installation in petsc external packages
      endif( DEFINED ENV{HDF5_ROOT} )
    endif( NOT HDF5_ROOT )

    find_package(HDF5 COMPONENTS CXX HL)
    if ( NOT HDF5_FOUND)
      set(HDF5_LIBRARY_DIR  ${HDF5_ROOT}/lib)
      set(HDF5_INCLUDE_DIRS ${HDF5_ROOT}/include)
    endif( NOT HDF5_FOUND)
  
  else ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR )
    message ( STATUS "PETSC_DIR not set, searching PETSc and SLEPc in the system" )

    if ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )  #Case fedora/redhat system install
      message ( STATUS "PETSC includes found in /usr/include/petsc/" )
      message ( STATUS "PETSC library  found in /usr/lib64" )
      set(PETSC_DIR /usr/)
      set(PETSC_INCLUDES  /usr/include/petsc /usr/include/petsc/petsc/mpiuni)
      set(PETSC_LIBRARIES /usr/lib64/libpetsc.so)

      file(STRINGS /usr/include/petsc/petscversion.h petscversion_file)
      FOREACH(line ${petscversion_file})
        string(FIND "${line}" "#define PETSC_VERSION_MAJOR      " pos)
        IF(NOT ${pos} EQUAL -1)
          string(LENGTH ${line} line_length_major)#length of the line
          MATH(EXPR petsc_version_major "${line_length_major}-1")#line ends with "x", that counts for 1 character
          string(SUBSTRING ${line} ${petsc_version_major} 1 PETSC_VERSION_MAJOR)
          break()
        ENDIF(NOT ${pos} EQUAL -1)
      ENDFOREACH(line ${petscversion_file})
      FOREACH(line ${petscversion_file})
        string(FIND "${line}" "#define PETSC_VERSION_MINOR      " pos)
        IF(NOT ${pos} EQUAL -1)
          string(LENGTH ${line} line_length_minor)#length of the line
          if(${line_length_minor} EQUAL ${line_length_major})
            MATH(EXPR petsc_version_minor "${line_length_minor}-1")#line ends with "y", that counts for 1 character
            string(SUBSTRING ${line} ${petsc_version_minor} 1 PETSC_VERSION_MINOR)
          break()
          else(${line_length_minor} EQUAL ${line_length_major})
            MATH(EXPR petsc_version_minor "${line_length_minor}-2")#line ends with "yy", that counts for two characters
            string(SUBSTRING ${line} ${petsc_version_minor} 2 PETSC_VERSION_MINOR)
          endif(${line_length_minor} EQUAL ${line_length_major})
        ENDIF(NOT ${pos} EQUAL -1)
      ENDFOREACH(line ${petscversion_file})
      set(PETSC_VERSION ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}) #3.8 for fedora 26 , 3.9 for fedora 29 , 3.10 for fedora 30 , 3.12 for fedora 32 , 3.14 for fedora 34 , 3.15 for fedora 35
      message(STATUS "PETSc detected, version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}" )
      
      set(CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, petsc-devel (>= 3.4)") # This is not fully working yet. 
      set(PETSC_INSTALL ${PETSC_DIR})

      #Define and search slepc variables
      if   ( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )
        message( STATUS "SLEPc includes found in /usr/include/slepc/" )
        message( STATUS "SLEPc library  found in /usr/lib64/slepc/" )
        set(SLEPC_DIR /usr/)
        set(SLEPC_INCLUDES  ${SLEPC_DIR}/include)
        set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
        file(STRINGS /usr/include/slepc/slepcversion.h slepcversion_file)
        FOREACH(line ${slepcversion_file})
          string(FIND "${line}" "#define SLEPC_VERSION_MAJOR      " pos)
          IF(NOT ${pos} EQUAL -1)
            string(LENGTH ${line} line_length_major)#length of the line
            MATH(EXPR slepc_version_major "${line_length_major}-1")#line ends with "x", that counts for 1 character
            string(SUBSTRING ${line} ${slepc_version_major} 1 SLEPC_VERSION_MAJOR)
            break()
          ENDIF(NOT ${pos} EQUAL -1)
        ENDFOREACH(line ${slepcversion_file})
        FOREACH(line ${slepcversion_file})
          string(FIND "${line}" "#define SLEPC_VERSION_MINOR      " pos)
          IF(NOT ${pos} EQUAL -1)
            string(LENGTH ${line} line_length_minor)#length of the line
            if(${line_length_minor} EQUAL ${line_length_major})
              MATH(EXPR slepc_version_minor "${line_length_minor}-1")#line ends with "y", that counts for 1 character
              string(SUBSTRING ${line} ${slepc_version_minor} 1 SLEPC_VERSION_MINOR)
            break()
            else(${line_length_minor} EQUAL ${line_length_major})
              MATH(EXPR slepc_version_minor "${line_length_minor}-2")#line ends with "yy", that counts for two characters
              string(SUBSTRING ${line} ${slepc_version_minor} 2 SLEPC_VERSION_MINOR)
            endif(${line_length_minor} EQUAL ${line_length_major})
          ENDIF(NOT ${pos} EQUAL -1)
        ENDFOREACH(line ${slepcversion_file})
        set(SLEPC_VERSION ${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}) #3.8 for fedora 26 , 3.9 for fedora 29 , 3.10 for fedora 30 , 3.12 for fedora 32 , 3.14 for fedora 34 , 3.15 for fedora 35
        message(STATUS "SLEPc detected, version is ${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}" )
        set (CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, slepc-devel (>= 3.4)") # This is not fully working yet. 
      else ( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )
        message( FATAL_ERROR "SLEPc not found in the system" )
      endif( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )

      #HDF5 to be found in the system
      if   ( NOT HDF5_ROOT AND NOT DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment
        find_package(HDF5 REQUIRED)
        set(CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-dev")
        set(CPACK_RPM_PACKAGE_REQUIRES   "hdf5-devel")
      endif( NOT HDF5_ROOT AND NOT DEFINED ENV{HDF5_ROOT} )

    #elseif ( IS_DIRECTORY "/usr/lib/petsc/include" AND EXISTS "/usr/lib/petsc/lib/libpetsc_real.so" ) #Case ubuntu/debian system install
    #  message ( STATUS "PETSc includes found in /usr/lib/petsc/include")
    #  message ( STATUS "PETSc library  found in /usr/lib/petsc/lib/libpetsc_real.so")
    #  set(PETSC_DIR /usr/lib/petsc/)
    #  set(PETSC_INCLUDES  ${PETSC_DIR}/include ${PETSC_DIR}/include/petsc/mpiuni)
    #  set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc_real.so)

    #  find_package (PETSc 3.4 REQUIRED)
    #  petsc_get_version ()
    #  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev   (>= 3.4)") # This is not fully working yet. 

    #  Define and search slepc variables
    #  if   ( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )
    #    message( STATUS "SLEPc includes found in /usr/lib/slepc/include" )
    #    message( STATUS "SLEPc library  found in /usr/lib/slepc/lib/libslepc_real.so" )
    #    set(SLEPC_DIR /usr/lib/slepc/)
    #    set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
    #    set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc_real.so)
    #    set (CPACK_DEBIAN_PACKAGE_REQUIRES   "${CPACK_DEBIAN_PACKAGE_REQUIRES}, slepc-devel (>= 3.4)") # This is not fully working yet. 
    #  else ( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )
    #    message( FATAL_ERROR "SLEPc not found in the system" )
    #  endif( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )

    #elseif ( IS_DIRECTORY "/usr/local/lib/python2.7/dist-packages/petsc") #Case ubuntu/debian system pip install

    else ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )  # No petsc in system, do compile petsc along with slepc as an external package
      message (STATUS "PETSC not found in the system")
      message (STATUS "PETSC will be downloaded and compiled from ${DOWNLOAD_PETSC}" )
      message (STATUS "SLEPC will be downloaded and compiled from ${DOWNLOAD_SLEPC}" )
      message (STATUS "HDF5  will be downloaded and compiled from ${DOWNLOAD_HDF5}" )

      #extract petsc tarball name
      string(LENGTH "${DOWNLOAD_PETSC}" tarball_url_length)#length of the tarball
      string(FIND   "${DOWNLOAD_PETSC}" "/" start_tarball_name REVERSE )# last occurence of "/"
      MATH(EXPR start_tarball_name  "${start_tarball_name}+1")#start after the occurence of "/"
      MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
      string(SUBSTRING ${DOWNLOAD_PETSC} ${start_tarball_name} ${tarball_name_length} PETSC_TARBALL_NAME)

      set(PETSC_DIR ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${PETSC_TARBALL_NAME}) # Location of the final install 
      set(PETSC_ARCH arch-linux-c-opt)
      set(PETSC_INSTALL ${CMAKE_INSTALL_PREFIX}/lib/petsc)#folder to copy petsc libraries and include files

      ExternalProject_Add (PETSc
      URL               ${DOWNLOAD_PETSC}   
      SOURCE_DIR        ${PETSC_DIR}
      BUILD_IN_SOURCE   TRUE
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${PETSC_INSTALL} --with-debugging=0 --download-f2cblaslapack=${DOWNLOAD_F2CBLASLAPACK} --with-mpi=0 --with-fc=0 --download-slepc=${DOWNLOAD_SLEPC} --download-hdf5=${DOWNLOAD_HDF5}   --with-petsc4py=yes --download-slepc-configure-arguments="--with-slepc4py=yes" # --download-openmpi=${DOWNLOAD_OPENMPI} --download-hypre=${DOWNLOAD_HYPRE} --download-mpi4py=${DOWNLOAD_MPI4PY}
      BUILD_COMMAND     make all
      TEST_COMMAND      make check
      INSTALL_COMMAND   make install
      INSTALL_DIR       ${PETSC_DIR}
      STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
       )

      message( STATUS "PETSc, SLEPc and HDF5 will be installed at ${PETSC_INSTALL}")

      string(LENGTH "${PETSC_TARBALL_NAME}" tarball_name_length)#length of the tarball
      string(FIND   "${PETSC_TARBALL_NAME}" "-" start_version_name REVERSE )# last occurence of "-"
      MATH(EXPR start_version_name  "${start_version_name}+1")#start after the occurence of "/"
      string(SUBSTRING ${PETSC_TARBALL_NAME} ${start_version_name} ${tarball_name_length} PETSC_VERSION)
      message( STATUS "PETSc version ${PETSC_VERSION}  will be installed" )

      set(PETSC_INCLUDES  ${PETSC_DIR}/${PETSC_ARCH}/include ${PETSC_DIR}/include)
      set(PETSC_LIBRARIES ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)

      #extract slepc tarball name 
      string(LENGTH "${DOWNLOAD_SLEPC}" tarball_url_length)#length of the tarball
      string(FIND   "${DOWNLOAD_SLEPC}" "/" start_tarball_name REVERSE )# last occurence of "/"
      MATH(EXPR start_tarball_name  "${start_tarball_name}+1")#start after the occurence of "/"
      MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
      string(SUBSTRING ${DOWNLOAD_SLEPC} ${start_tarball_name} ${tarball_name_length} SLEPC_TARBALL_NAME)

      #extract slepc version
      string(LENGTH "${SLEPC_TARBALL_NAME}" tarball_name_length)#length of the tarball
      string(FIND   "${SLEPC_TARBALL_NAME}" "-" start_version_name REVERSE )# last occurence of "-"
      MATH(EXPR start_version_name  "${start_version_name}+1")#start after the occurence of "/"
      string(SUBSTRING ${SLEPC_TARBALL_NAME} ${start_version_name} ${tarball_name_length} SLEPC_VERSION)
      message( STATUS "SLEPc version ${SLEPC_VERSION}  will be installed" )

      #define slepc variables
      set(SLEPC_DIR ${PETSC_INSTALL})
      set(SLEPC_INCLUDES ${SLEPC_DIR}/include  ${PETSC_DIR}/${PETSC_ARCH}/externalpackages/${SLEPC_TARBALL_NAME}/include )
      set(SLEPC_LIBRARIES ${PETSC_DIR}/${PETSC_ARCH}/externalpackages/${SLEPC_TARBALL_NAME}/installed-arch-linux2-c-opt/lib/libslepc.so)
      message( STATUS "SLEPc includes ${SLEPC_INCLUDES}")

      #define hdf5 variables
      if   ( NOT HDF5_ROOT )#hdf5 is not defined in cmake arguments
        if  ( DEFINED ENV{HDF5_ROOT} )#hdf5 is defined in the environment
          set(HDF5_ROOT $ENV{HDF5_ROOT})
        else( DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment          
          set(HDF5_ROOT ${PETSC_INSTALL})#HDF5 to be found in petsc installation
        endif( DEFINED ENV{HDF5_ROOT} )
      endif( NOT HDF5_ROOT )
      set(HDF5_LIBRARY_DIR ${HDF5_ROOT}/lib)
      set(HDF5_INCLUDE_DIRS ${HDF5_ROOT}/include)

      #define MPI variables
      #set(SOLVERLAB_WITH_MPI ON)
      #set(MPI_HOME    ${PETSC_INSTALL} )
      #set(MPI_LIBRARY ${PETSC_INSTALL}/lib/libmpi.so)
      #set(MPIEXEC     ${PETSC_INSTALL}/bin/mpiexec)
      #set(MPI4PY_ROOT_DIR    ${PETSC_INSTALL}/lib/mpi4py )

    endif ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )
  endif ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR )

if   ( ${PETSC_VERSION} VERSION_GREATER 3.5 )
  add_definitions(-DPETSC_VERSION_GREATER_3_5)
endif( ${PETSC_VERSION} VERSION_GREATER 3.5 )

string(REPLACE ";" ":"  PETSC_INCLUDES_INSTALL "${PETSC_INCLUDES}")# use colon instead of semicolon in environment file env_SOLVERLAB.sh

add_library(petsc   SHARED IMPORTED)
set_property(TARGET petsc  PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
add_library(slepc   SHARED IMPORTED)
set_property(TARGET slepc  PROPERTY IMPORTED_LOCATION ${SLEPC_LIBRARIES})

######################################################
# Detection or compilation of MEDFile                #
######################################################
if( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})

  if   ( NOT MEDFILE_ROOT_DIR )
    set(MEDFILE_ROOT_DIR $ENV{MEDFILE_ROOT_DIR} )
  endif( NOT MEDFILE_ROOT_DIR )

  find_package (MEDFile REQUIRED)

else( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})

  message(STATUS "MED will be downloaded and installed from ${DOWNLOAD_MED}")
  set(MACHINE PCLINUX)
  set(MEDFILE_DEFINITIONS "-D${MACHINE} -DMEDFILE_INSTALL_DOC=OFF")#Comment ne pas compiler la doc ?
  set(MEDFILE_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/lib/med) # Location of the final install
  
  #extraction of the tarball archive name
  string(LENGTH "${DOWNLOAD_MED}" tarball_url_length)#length of the name of the tarball
  string(FIND "${DOWNLOAD_MED}" "/" start_tarball_name REVERSE )# last occurence of "/"
  MATH(EXPR start_tarball_name "${start_tarball_name}+1")#start after the occurence of "/"
  MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
  string(SUBSTRING ${DOWNLOAD_MED} ${start_tarball_name} ${tarball_name_length} MED_TARBALL_NAME)

  ExternalProject_Add (MED
        URL               ${DOWNLOAD_MED}  #location of med tarball
        SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${MED_TARBALL_NAME}
        BUILD_IN_SOURCE   TRUE
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${MEDFILE_ROOT_DIR}  --with-hdf5=${HDF5_ROOT} --with-hdf5-include=${HDF5_ROOT}/include --with-hdf5-lib=${HDF5_LIBRARY_DIR} --with-hdf5-bin=${HDF5_ROOT}/bin  --with-swig=yes --disable-fortran #CC=${PETSC_INSTALL}/bin/mpicc
        BUILD_COMMAND     make
        INSTALL_COMMAND   make install
        INSTALL_DIR       ${MEDFILE_ROOT_DIR}
        STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
        LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
        LOG_UPDATE        TRUE           # Wrap update in script to log output
        LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
        LOG_BUILD         TRUE           # Wrap build in script to log output
        LOG_TEST          TRUE           # Wrap test in script to log output
        LOG_INSTALL       TRUE           # Wrap install in script to log output
  )
  
  set(MEDFILE_INCLUDE_DIRS ${MEDFILE_ROOT_DIR}/include)# Nécessaire pour le medloader et les sous-dossiers mesh, et IJKMesh
  set(MEDFILE_LIBRARIES    ${MEDFILE_ROOT_DIR}/lib)# Nécessaire pour le env_SOLVERLAB.sh
  set(MEDFILE_C_LIBRARIES med medC)# Nécessaire pour le medloader

endif( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})

add_library(med   SHARED IMPORTED)
set_property(TARGET med  PROPERTY IMPORTED_LOCATION ${MEDFILE_ROOT_DIR}/lib/libmed.so)
add_library(medC   SHARED IMPORTED)
set_property(TARGET medC PROPERTY IMPORTED_LOCATION ${MEDFILE_ROOT_DIR}/lib/libmedC.so)

string(REPLACE ";" ":"  MEDFILE_LIBRARIES_INSTALL "${MEDFILE_LIBRARIES}")# use colon instead of semicolon in environment file env_SOLVERLAB.sh

######################################################
# Detection or compilation of MEDCoupling            #
######################################################
if( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR})

  if   ( NOT MEDCOUPLING_ROOT_DIR )
    set(MEDCOUPLING_ROOT_DIR $ENV{MEDCOUPLING_ROOT_DIR} )
  endif( NOT MEDCOUPLING_ROOT_DIR )

  message (STATUS "Seeking MEDCoupling library in ${MEDCOUPLING_ROOT_DIR}")

  set(MEDCoupling_DIR  ${MEDCOUPLING_ROOT_DIR}/cmake_files)
  find_package(MEDCoupling CONFIG REQUIRED)   # CONFIG is important !
  message( STATUS "MEDCoupling found. Version is ${MEDCoupling_VERSION}")

else( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR})
  message(STATUS "MEDCoupling will be downloaded and installed from ${DOWNLOAD_MEDCOUPLING}")
  set(MEDCOUPLING_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/lib/medcoupling) # Location of the final install
  SET(ENV{HDF5_ROOT} ${HDF5_ROOT})#MEDCoupling install process seems to require an environment variable

  #extraction of the tarball archive name
  string(LENGTH "${DOWNLOAD_MEDCOUPLING}" tarball_url_length)#length of the tarball
  string(FIND   "${DOWNLOAD_MEDCOUPLING}" "/" start_tarball_name REVERSE )# last occurence of "/"
  MATH(EXPR start_tarball_name "${start_tarball_name}+1")#start after the occurence of "/"
  MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
  string(SUBSTRING ${DOWNLOAD_MEDCOUPLING} ${start_tarball_name} ${tarball_name_length} MEDCOUPLING_TARBALL_NAME)
  string(SUBSTRING ${MEDCOUPLING_TARBALL_NAME} 12 ${tarball_name_length} MEDCoupling_VERSION)

  message(STATUS "MEDCoupling will be installed at ${MEDCOUPLING_ROOT_DIR}")

  ExternalProject_Add (MEDCoupling
        URL               ${DOWNLOAD_MEDCOUPLING}  #location of medcoupling tarball
        SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${MEDCOUPLING_TARBALL_NAME}
        BUILD_IN_SOURCE   TRUE
        CONFIGURE_COMMAND cmake <SOURCE_DIR>/${MEDCOUPLING_TARBALL_NAME} -DCMAKE_INSTALL_PREFIX=${MEDCOUPLING_ROOT_DIR} -DCONFIGURATION_ROOT_DIR=<SOURCE_DIR>/configuration-${MEDCoupling_VERSION} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ENABLE_PARTITIONER=OFF -DMEDCOUPLING_PARTITIONER_METIS=OFF -DMEDCOUPLING_PARTITIONER_SCOTCH=OFF -DMEDCOUPLING_ENABLE_RENUMBER=OFF -DMEDCOUPLING_BUILD_DOC=OFF -DHDF5_ROOT_DIR=${HDF5_ROOT} -DMEDCOUPLING_USE_64BIT_IDS=ON #-DSALOME_USE_MPI=ON
        BUILD_COMMAND     make
        INSTALL_COMMAND   make install
        INSTALL_DIR       ${MEDCOUPLING_ROOT_DIR}
        STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
        LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
        LOG_UPDATE        TRUE           # Wrap update in script to log output
        LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
        LOG_BUILD         TRUE           # Wrap build in script to log output
        LOG_TEST          TRUE           # Wrap test in script to log output
        LOG_INSTALL       TRUE           # Wrap install in script to log output
  )

  add_library(medloader   SHARED IMPORTED )
  set_property(TARGET medloader   PROPERTY IMPORTED_LOCATION ${MEDCOUPLING_ROOT_DIR}/lib/libmedloader.so)

  SET(MEDCOUPLING_USE_64BIT_IDS  ON)

endif( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR} )

set(MEDCOUPLING_INCLUDE_DIR ${MEDCOUPLING_ROOT_DIR}/include)
set(MEDCOUPLING_LIBRARIES   ${MEDCOUPLING_ROOT_DIR}/lib    )#for environment file env_SOLVERLAB.sh 

add_library(medcoupling SHARED IMPORTED )
set_property(TARGET medcoupling PROPERTY IMPORTED_LOCATION ${MEDCOUPLING_ROOT_DIR}/lib/libmedcoupling.so)

if   ( ${MEDCoupling_VERSION} VERSION_GREATER 9.4 )
  add_definitions(-DMEDCoupling_VERSION_VERSION_GREATER_9_4)
endif( ${MEDCoupling_VERSION} VERSION_GREATER 9.4 )

if   (TARGET MED AND TARGET PETSc)
  ExternalProject_Add_StepDependencies( MED build PETSc )#PETSc doit être compilé avant MED car c'est PETSc qui compile HDF5
endif(TARGET MED AND TARGET PETSc)

if   (TARGET MED AND TARGET MEDCoupling)
  ExternalProject_Add_StepDependencies( MEDCoupling build MED )#MED doit être compilé avant MEDCoupling
endif(TARGET MED AND TARGET MEDCoupling)
 
if   (TARGET PETSc AND TARGET MEDCoupling)
  ExternalProject_Add_StepDependencies( MEDCoupling build PETSc )#PETSc doit être compilé avant MEDCoupling (car il contient hdf5)
endif(TARGET PETSc AND TARGET MEDCoupling)

IF(MEDCOUPLING_USE_64BIT_IDS)
  SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DMEDCOUPLING_USE_64BIT_IDS)
ENDIF(MEDCOUPLING_USE_64BIT_IDS)
 
####################
# Final operations #
####################                                  

# Paraview variables to choose python version and for env_SOLVERLAB.sh  
if  ( PARAVIEW_ROOT_DIR OR DEFINED ENV{PARAVIEW_ROOT_DIR})
  if  ( PARAVIEW_ROOT_DIR )
    set( PARAVIEW_INCLUDE_DIR   ${PARAVIEW_ROOT_DIR}/include)
    set( PARAVIEW_LIBRARIES_DIR ${PARAVIEW_ROOT_DIR}/lib)
  else( PARAVIEW_ROOT_DIR )
    set( PARAVIEW_ROOT_DIR      $ENV{PARAVIEW_ROOT_DIR} )
    set( PARAVIEW_INCLUDE_DIR   $ENV{PARAVIEW_ROOT_DIR}/include)
    set( PARAVIEW_LIBRARIES_DIR $ENV{PARAVIEW_ROOT_DIR}/lib)
  endif( PARAVIEW_ROOT_DIR )
else( PARAVIEW_ROOT_DIR OR DEFINED ENV{PARAVIEW_ROOT_DIR})
  message(STATUS "PARAVIEW_ROOT_DIR not set. Trying to detect paraview in the system.")
  set( PARAVIEW_INCLUDE_DIR /usr/include/paraview/)
  set( PARAVIEW_LIBRARIES_DIR /usr/lib/paraview/:/usr/lib64/paraview/:/usr/lib64/vtk/:/usr/lib/vtk/:/usr/lib/x86_64-linux-gnu)
endif( PARAVIEW_ROOT_DIR OR DEFINED ENV{PARAVIEW_ROOT_DIR})

if  ( DEFINED PARAVIEW_VERSION OR DEFINED ENV{PARAVIEW_VERSION} )
  if   ( NOT DEFINED PARAVIEW_VERSION )
    set( PARAVIEW_VERSION  $ENV{PARAVIEW_VERSION})
  endif( NOT DEFINED PARAVIEW_VERSION )
  message(STATUS "ParaView version is ${PARAVIEW_VERSION}" )
else( DEFINED PARAVIEW_VERSION OR DEFINED ENV{PARAVIEW_VERSION} )
  set( PARAVIEW_SUFFIXES paraview-5.0 paraview-5.1 paraview-5.2 paraview-5.3 paraview-5.4 paraview-5.5 paraview-5.6 paraview-5.7 paraview-5.8 paraview-5.9 )
  find_file( vtkConfig_file vtkPVConfig.h HINTS ${PARAVIEW_INCLUDE_DIR} PATH_SUFFIXES ${PARAVIEW_SUFFIXES})

  #extraction of the paraview version"
  IF   ( ${vtkConfig_file} STREQUAL "vtkConfig_file-NOTFOUND" )
    message(WARNING "Could not find ParaView configuration file vtkPVConfig.h in folder ${PARAVIEW_INCLUDE_DIR} and subdirectories  ${PARAVIEW_SUFFIXES}. Postprocessing may not work")
    set( PARAVIEW_VERSION  "5.8")
  ELSE ( ${vtkConfig_file} STREQUAL "vtkConfig_file-NOTFOUND" )
    message( STATUS "Found vtkPVConfig.h in ${vtkConfig_file}")
    file(STRINGS ${vtkConfig_file} vtkConfig)
    FOREACH(line ${vtkConfig})
      string(FIND "${line}" "#define PARAVIEW_VERSION_FULL " pos)
      IF(NOT ${pos} EQUAL -1)
        string(LENGTH ${line} line_length)#length of the tarball
        MATH(EXPR start_pv_version "${line_length}-6")#line ends with "x.y.z", that counts for 7 characters
        string(SUBSTRING ${line} ${start_pv_version} 5 PARAVIEW_VERSION)
        break()
      ENDIF(NOT ${pos} EQUAL -1)
    ENDFOREACH(line ${vtkConfig})
    message(STATUS "ParaView detected, version is ${PARAVIEW_VERSION}" )
  ENDIF( ${vtkConfig_file} STREQUAL "vtkConfig_file-NOTFOUND" )
endif( DEFINED PARAVIEW_VERSION OR DEFINED ENV{PARAVIEW_VERSION} )

#Choose python 2 or 3  
IF   ( ${PARAVIEW_VERSION} STRGREATER "5.6" OR ${PARAVIEW_VERSION} STREQUAL "5.6" )#Use python 3, use VERSION_GREATER_EQUAL if cmake >=3.7
  SET(PYTHON2OR3 "3")
  message(STATUS "PARAVIEW_VERSION greater than 5.6. We need python3 for postprocessing scripts.")
ELSE ( ${PARAVIEW_VERSION} STRGREATER "5.6" OR ${PARAVIEW_VERSION} STREQUAL "5.6" )#Use python 2
  SET(PYTHON2OR3 "2")
  message(STATUS "PARAVIEW_VERSION smaller than 5.6. We need python2 for postprocessing scripts.")
ENDIF( ${PARAVIEW_VERSION} STRGREATER "5.6" OR ${PARAVIEW_VERSION} STREQUAL "5.6" )

# Find python
IF   ( PYTHON_ROOT_DIR OR DEFINED ENV{PYTHON_ROOT_DIR})
  IF  ( PYTHON_ROOT_DIR )
    set(Python_ROOT_DIR ${PYTHON_ROOT_DIR} ) 
  ELSE( PYTHON_ROOT_DIR )
    set(Python_ROOT_DIR $ENV{PYTHON_ROOT_DIR} ) 
  ENDIF( PYTHON_ROOT_DIR)
ENDIF( PYTHON_ROOT_DIR OR DEFINED ENV{PYTHON_ROOT_DIR})

IF   (${CMAKE_VERSION} VERSION_GREATER "3.12.0")
  find_package(Python ${PYTHON2OR3} REQUIRED COMPONENTS Interpreter Development )
  SET(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
ELSE (${CMAKE_VERSION} VERSION_GREATER "3.12.0")
  find_package(PythonInterp ${PYTHON2OR3} REQUIRED )
  find_package(PythonLibs   ${PYTHON2OR3} REQUIRED )
  set(Python_LIBRARIES     ${PYTHON_LIBRARIES})
  set(Python_INCLUDE_DIRS  ${PYTHON_INCLUDE_DIRS})
  set(Python_VERSION       ${PYTHON_VERSION_STRING})
  set(Python_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
  set(Python_VERSION_MINOR ${PYTHON_VERSION_MINOR})
ENDIF(${CMAKE_VERSION} VERSION_GREATER "3.12.0")

message(STATUS "Python version found is ${Python_VERSION}")

#Matplotlib root dir
IF   ( NOT DEFINED MATPLOTLIB_ROOT_DIR )
    set(MATPLOTLIB_ROOT_DIR $ENV{MATPLOTLIB_ROOT_DIR} )
ENDIF( NOT DEFINED MATPLOTLIB_ROOT_DIR )

#PyQt5 root dir
IF   ( NOT DEFINED PYQT5_ROOT_DIR )
    set(PYQT5_ROOT_DIR $ENV{PYQT5_ROOT_DIR} )
ENDIF( NOT DEFINED PYQT5_ROOT_DIR )

#Paraview library paths
set (PV_LIB_DIR    /usr/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/dist-packages/paraview/:/usr/lib64/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vtk/:${MATPLOTLIB_ROOT_DIR}/lib/:${PARAVIEW_LIBRARIES_DIR})
set (PV_PYTHON_DIR /usr/lib/python${Python_VERSION_MAJOR}/dist-packages/:/usr/lib/python${Python_VERSION_MAJOR}/dist-packages/paraview/:/usr/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/dist-packages/paraview/:/usr/lib64/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vtk/:/usr/lib64/paraview/site-packages/:/usr/lib64/paraview/site-packages/paraview/:/usr/lib64/paraview/site-packages/vtk/:/usr/lib64/paraview/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/:/usr/lib64/paraview/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/paraview/:/usr/lib64/paraview/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vtkmodules:/usr/lib/paraview/site-packages/paraview:/usr/lib/paraview/site-packages/paraview/vtk:${PARAVIEW_ROOT_DIR}/lib/:${PARAVIEW_ROOT_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/:${PARAVIEW_ROOT_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/vtkmodules/:${MATPLOTLIB_ROOT_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages:${PYQT5_ROOT_DIR}:${PYQT5_ROOT_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages)

# Detection of python-paraview
string(REPLACE ":" ";"  PV_PYTHON_DIR_PATHS "${PV_PYTHON_DIR}")# use semicolon instead of colon in PATH search
find_file(python-paraview simple.py PATHS ${PV_PYTHON_DIR_PATHS})
if( ${python-paraview} STREQUAL "python-paraview-NOTFOUND" )
  message(WARNING "Python Paraview library not found. Post processing will fail in SOLVERLAB scripts.")
  set(PRELOAD_NETCDF_ON_UBUNTU20 OFF)
else( ${python-paraview} STREQUAL "python-paraview-NOTFOUND" )
  if(    ${python-paraview} STREQUAL "/usr/lib/python3/dist-packages/paraview/simple.py" )
    message ( STATUS "Using system Paraview python libraries on Ubuntu 20.04 requires the presence and preload of library netcdf (bug in linkage)." )
    if(    EXISTS "/usr/lib/x86_64-linux-gnu/libnetcdf.so" )
      message ( STATUS "Library netcdf found in the system. Paraview python libraries should be functional." )
      set(PRELOAD_NETCDF_ON_UBUNTU20 ON)
    else(  EXISTS "/usr/lib/x86_64-linux-gnu/libnetcdf.so" )
      message ( WARNING "Library netcdf not found in system. Paraview python libraries may not be functional (bug in linkage)." )
      set(PRELOAD_NETCDF_ON_UBUNTU20 OFF)
    endif( EXISTS "/usr/lib/x86_64-linux-gnu/libnetcdf.so" )
  else(  ${python-paraview} STREQUAL "/usr/lib/python3/dist-packages/paraview/simple.py" )
    message(STATUS "Python Paraview library found. Post processing will work in SOLVERLAB scripts.")
    set(PRELOAD_NETCDF_ON_UBUNTU20 OFF)
  endif( ${python-paraview} STREQUAL "/usr/lib/python3/dist-packages/paraview/simple.py" )
endif( ${python-paraview} STREQUAL "python-paraview-NOTFOUND" )

# Find SWIG
IF   ( DEFINED SWIG_ROOT_DIR OR DEFINED ENV{SWIG_ROOT_DIR})
  IF   ( DEFINED SWIG_ROOT_DIR )
    set(SWIG_EXECUTABLE    ${SWIG_ROOT_DIR}/bin/swig ) 
  ELSE ( DEFINED SWIG_ROOT_DIR )
    set(SWIG_EXECUTABLE $ENV{SWIG_ROOT_DIR}/bin/swig ) 
  ENDIF( DEFINED SWIG_ROOT_DIR)
ENDIF( DEFINED SWIG_ROOT_DIR OR DEFINED ENV{SWIG_ROOT_DIR})
find_package (SWIG 3.0 REQUIRED)

# Find Doxygen
if (SOLVERLAB_WITH_DOCUMENTATION)                                                                           #
 IF   ( DEFINED DOXYGEN_ROOT_DIR OR DEFINED ENV{DOXYGEN_ROOT_DIR})
   IF   ( DEFINED DOXYGEN_ROOT_DIR )
     set(DOXYGEN_EXECUTABLE    ${DOXYGEN_ROOT_DIR}/bin/doxygen ) 
   ELSE ( DEFINED DOXYGEN_ROOT_DIR )
     set(DOXYGEN_EXECUTABLE $ENV{DOXYGEN_ROOT_DIR}/bin/doxygen ) 
   ENDIF( DEFINED DOXYGEN_ROOT_DIR)
 ENDIF( DEFINED DOXYGEN_ROOT_DIR OR DEFINED ENV{DOXYGEN_ROOT_DIR})
  find_package (Doxygen)                                                                                    #
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, doxygen, graphviz, mscgen")
  set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, doxygen, graphviz, mscgen")
  if   (NOT DOXYGEN_FOUND)
    message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
  endif(NOT DOXYGEN_FOUND)
endif(SOLVERLAB_WITH_DOCUMENTATION)			

if (SOLVERLAB_WITH_TESTS)                                                                                  #
  enable_testing ()
 IF   ( DEFINED CPPUNIT_ROOT_DIR OR DEFINED ENV{CPPUNIT_ROOT_DIR})
   IF   ( DEFINED CPPUNIT_ROOT_DIR )
     set(CPPUNIT_ROOT_USER    ${CPPUNIT_ROOT_DIR} ) 
   ELSE ( DEFINED CPPUNIT_ROOT_DIR )
     set(CPPUNIT_ROOT_USER $ENV{CPPUNIT_ROOT_DIR} ) 
   ENDIF( DEFINED CPPUNIT_ROOT_DIR)
 ENDIF( DEFINED CPPUNIT_ROOT_DIR OR DEFINED ENV{CPPUNIT_ROOT_DIR})
  find_package (CPPUNIT REQUIRED)
endif(SOLVERLAB_WITH_TESTS)			
				    
# Enter subdirectories
set(CDMATH_WITH_POSTPRO ${SOLVERLAB_WITH_POSTPRO} )
set(CDMATH_WITH_PETSC ${SOLVERLAB_WITH_PETSC} )
set(CDMATH_WITH_DOCUMENTATION ${SOLVERLAB_WITH_DOCUMENTATION} )
set(CDMATH_WITH_PYTHON ${SOLVERLAB_WITH_PYTHON} )
set(CDMATH_WITH_PACKAGE ${SOLVERLAB_WITH_PACKAGE} )
set(CDMATH_WITH_TESTS ${SOLVERLAB_WITH_TESTS} )
add_subdirectory (CDMATH)
if(SOLVERLAB_WITH_COREFLOWS)
  set(COREFLOWS_WITH_DOCUMENTATION ${SOLVERLAB_WITH_DOCUMENTATION} )
  set(COREFLOWS_WITH_PYTHON ${SOLVERLAB_WITH_PYTHON} )
  set(COREFLOWS_WITH_GUI ${SOLVERLAB_WITH_GUI} )
  set(COREFLOWS_WITH_PACKAGE ${SOLVERLAB_WITH_PACKAGE} )
  set(COREFLOWS_WITH_TESTS ${SOLVERLAB_WITH_TESTS} )
  add_subdirectory (CoreFlows)
endif(SOLVERLAB_WITH_COREFLOWS)

# Configuration file
configure_file(
    ${PROJECT_SOURCE_DIR}/env_SOLVERLAB.sh
    ${CMAKE_INSTALL_PREFIX}/env_SOLVERLAB.sh
    @ONLY
)

