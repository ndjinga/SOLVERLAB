#!/bin/bash

export CDMATH_INSTALL=@CMAKE_INSTALL_PREFIX@
source @CMAKE_INSTALL_PREFIX@/env_SOLVERLAB.sh

export CoreFlows_INSTALL=@CMAKE_INSTALL_PREFIX@
export PETSC_DIR=@PETSC_DIR@
export PETSC_ARCH=@PETSC_ARCH@
export PETSC_INCLUDES=@PETSC_INCLUDES_INSTALL@

#------------------------------------------------------------------------------------------------------------------- 
export CoreFlows=$CoreFlows_INSTALL/bin/CoreFlowsMainExe
export LD_LIBRARY_PATH=$CoreFlows_INSTALL/lib:$CDMATH_INSTALL/lib:${PETSC_DIR}/${PETSC_ARCH}/lib:${MEDCOUPLING_LIBRARIES}:${MEDFILE_C_LIBRARIES}:${LD_LIBRARY_PATH}
export PYTHONPATH=$CoreFlows_INSTALL/lib/coreflows:$CoreFlows_INSTALL/bin/coreflows:$CDMATH_INSTALL/lib/cdmath:$CDMATH_INSTALL/bin/cdmath:$CDMATH_INSTALL/bin/cdmath/postprocessing:${PETSC_DIR}/${PETSC_ARCH}/lib:${MEDCOUPLING_LIBRARIES}:${MEDFILE_C_LIBRARIES}:${PYTHONPATH}

