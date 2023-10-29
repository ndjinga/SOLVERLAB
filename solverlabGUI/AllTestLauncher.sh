#!/bin/bash

# params are:
# AllTestLauncher.py --rootPath=$1 --pattern=$2 --type=$3
# example: AllTestLauncher.sh . 'test_???_*.py' std

if [ -z "$SOLVERLABGUI_ROOT_DIR" ]
then
  echo "ERROR: You need set env var SOLVERLABGUI_ROOT_DIR before, continue at your own risks..."
  SOLVERLABGUI_ROOT_DIR=$(pwd)
fi

if [ -z "$1" ]
then
  PAR1=${SOLVERLABGUI_ROOT_DIR}
  echo "set default --rootPath=${PAR1}"
else
  PAR1=$1
fi

if [ -z "$2" ]
then
  PAR2="test_???_*.py"
  echo "set default --pattern=${PAR2}"
else
  PAR2=$2
fi

if [ -z "$3" ]
then
  PAR3="std"
  echo "set default --type=${PAR3}"
else
  PAR3=$3
fi


# may be cache cause problems
find ${SOLVERLABGUI_ROOT_DIR} -name "__pycache__" -type d  -exec rm -rf {} \;

ATL="${SOLVERLABGUI_ROOT_DIR}/unittestpy/AllTestLauncher.py"

if [ -f ${ATL} ]

then   # case standalone solverlabgui installation
  echo "try 'AllTestLauncher.py --help' for explanations"
  echo ${ATL} "--rootPath=${PAR1} --pattern='"${PAR2}"' --type=${PAR3}"
  ${ATL} --rootPath=${PAR1} --pattern=${PAR2} --type=${PAR3}

else   # case salome matix sat solverlabgui installation
  ATL="${PACKAGESPY_ROOT_DIR}/packagespy/unittestpy/AllTestLauncher.py"
  echo "try 'AllTestLauncher.py --help' for explanations"

  echo ${ATL} "--rootPath=${PAR1} --pattern='"${PAR2}"' --type=${PAR3}"
  ${ATL} --rootPath=${PAR1} --pattern=${PAR2} --type=${PAR3}

  echo ${ATL} "--rootPath=${PACKAGESPY_ROOT_DIR}/packagespy --pattern='"${PAR2}"' --type=${PAR3}"
  ${ATL} --rootPath=${PACKAGESPY_ROOT_DIR}/packagespy --pattern=${PAR2} --type=${PAR3}
fi



echo "#### END of SOLVERLABGUI AllTestLauncher.sh"
