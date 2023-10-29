#!/usr/bin/env bash

# 211103
# this script is for command line solverlab native ubuntu20 compilation,
# needs salome 9.7 libraries as prerequisites (PETC MEDCOUPLING HDF5)
# and launch test to check install (or not),
# comprehension and use at our own risks


function svl_helper {
    echo "Bash script to install solverlab quickly from SOLVERLAB_SRC to install directory.

args (only one possible):    
    -h | --help      : this help
    
environment args :    
    svl_nproc=N      : lance make -jN where N is the number of core
    svl_set_env      : lanch env_solverlab.sh after installation complete
    svl_test         : set env_solverlab and launch ctest -I 1,5
    svl_test_verbose : set env_solverlab and launch ctest -v -I 1,5
    svl_clean        : rm SOURCE and download them from github
    
Example :
    >>svl_nproc=12 \\
      svl_test=1 \\
      svl_clean=1 \\
      ./compile_solverlab.sh    
    "
}

function in_ko () {
  # Red bold
  echo
  echo -e "\e[31m\e[1m"$@"\e[0m"
  echo
}

function in_ok () {
  # green bold
  echo
  echo -e "\e[32m\e[1m"$@"\e[0m"
  echo
}

function in_magenta () {
  # Magenta bold
  echo
  echo -e "\e[35m\e[1m"$@"\e[0m"
  echo
}



function set_env_svl {
    export ROOTDIR=/volatile/catB/${USER}
    export SOLVERLABDIR=${ROOTDIR}/SOLVERLAB
    export SALOMEDIR=${ROOTDIR}/salome
    export SALOMEVERSION=SALOME-9.7.0-native-UB20.04-SRC
    export SALOME_ROOT_DIR=${SALOMEDIR}/${SALOMEVERSION}
    export PETSC_DIR=${SALOME_ROOT_DIR}/INSTALL/petsc
    export PETSC_ARCH=""
    export MEDFILE_ROOT_DIR=${SALOME_ROOT_DIR}/INSTALL/medfile
    export MEDCOUPLING_ROOT_DIR=${SALOME_ROOT_DIR}/INSTALL/MEDCOUPLING
    export HDF5_ROOT=${SALOME_ROOT_DIR}/INSTALL/PRODUCTS
}

function compile_svl {
    in_magenta '*** Compilation Solverlab'
    mkdir build

    cd build

    cmake ../SOLVERLAB_SRC \
	    -DCMAKE_INSTALL_PREFIX=../install \
	    -DCMAKE_BUILD_TYPE=Release  \
	    -G"Eclipse CDT4 - Unix Makefiles"  \
	    -D_ECLIPSE_VERSION=4.3  \
	    -DSOLVERLAB_WITH_DOCUMENTATION=ON  \
	    -DPETSC_DIR=${PETSC_DIR}  \
	    -DPETSC_ARCH=${PETSC_ARCH}  \
	    -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR}  \
	    -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR}  \
	    -DSOLVERLAB_WITH_GUI=ON

    make -j$svl_nproc

    make install
}

############
# Main
###########
export svl_debug=${svl_debug:-0}

# echoes only for debug this bash script
if [ $svl_debug == "1" ]
then
  set -x
fi

if [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
  svl_helper
  exit 0
fi


export svl_test=${svl_test:-0}
export svl_clean=${svl_clean:-0}
export svl_test_verbose=${svl_test_verbose:-0}
export svl_set_env=${svl_set_env:-0}
export svl_nproc=${svl_nproc:-1}


in_magenta '*** environ variables'
env | grep "svl_" | sort
echo

if [ $svl_clean == "1" ];then
    in_magenta '*** git clone source'
    rm -rf SOLVERLAB_SRC
    git clone --branch master https://github.com/ndjinga/SOLVERLAB.git SOLVERLAB_SRC
fi
    
rm -rf build install

set_env_svl
compile_svl

if [ $svl_set_env == "1" ];then
    in_magenta '*** set solverlab env'
    source ../install/env_SOLVERLAB.sh
fi

if [ $svl_test == "1" ]; then
    in_magenta '*** test 1,5'
    source ../install/env_SOLVERLAB.sh
    ctest -I 1,5
fi

if [ $svl_test_verbose == "1" ]; then
    in_magenta '*** test 1,5 verbose'
    source ../install/env_SOLVERLAB.sh
    ctest -V -I 1,5
fi

set +x

exit 0
