#!/bin/bash

# params are:
# AllTestLauncher.py --rootPath=$1 --pattern=$2 --type=$3
# example: AllTestLauncher.sh ./SOURCES/PACKAGESPY 'da*Test.py' std

if [ -z "$PACKAGESPY_ROOT_DIR" ]
then
  echo "You need command 'matix context' or 'salome context' before."
  echo "Or set environ variable PACKAGESPY_ROOT_DIR before."
  exit 1
fi

if [ -z "$1" ]
then
  PAR1=${PACKAGESPY_ROOT_DIR}/pythonAppliMatix
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

ATL="${PACKAGESPY_ROOT_DIR}/pythonAppliMatix/unittestpy/AllTestLauncher.py"
# HTMLOUT=/tmp/${USER}/AllTestLauncher_PACKAGESPY.html

echo "try 'AllTestLauncher.py --help' for explanations"
echo $ATL "--rootPath=${PAR1} --pattern='"${PAR2}"' --type=${PAR3} > ${HTMLOUT}"
# --debug
$ATL --rootPath=${PAR1} --pattern=${PAR2} --type=${PAR3} # > ${HTMLOUT}

# firefox $HTMLOUT &

echo "#### END of AllTestLauncherPACKAGESPY.sh"
