### Compilation SOLVERLAB from github

See `Compile and install SOLVERLAB` at https://github.com/ndjinga/SOLVERLAB

#### is221713 ubuntu 20 example

Using some SOLVERLAB as SALOME libraries as prerequisites.

Set differents ${USER} (wambeke or yacine-ym268439)

```
function set_env_solverlab {
    export ROOTDIR=/volatile/catB/${USER}
    export SOLVERLABDIR=${ROOTDIR}/solverlab
    export SALOMEDIR=${ROOTDIR}/salome
    export SALOMEVERSION=SALOME-9.7.0-native-UB20.04-SRC
    export SALOME_ROOT_DIR=${SALOMEDIR}/${SALOMEVERSION}
}
```


##### install SALOME 970

Optionally

```
mkdir ${SALOMEDIR}
cd ${SALOMEDIR}
wget ftp://ftp.cea.fr/pub/salome/targz/SALOME-9.7.0/SALOME-9.7.0-native-UB20.04-SRC.tar.gz
tar -xf SALOME-9.7.0-native-UB20.04-SRC.tar.gz
cd ${SALOMEDIR}/SALOME-9.7.0-native-UB20.04-SRC
which python      # /usr/bin/python
python -V         # Python 3.8.10
./sat config SALOME-9.7.0-native --check_system 
pluma README      # see infos
# to compilation based on the binaries used as prerequisites
./install_bin.sh  # copies BINARIES-UB20.04 into INSTALL
# launch salome
${SALOMEDIR}/SALOME-9.7.0-native-UB20.04-SRC/mesa_salome -k
```

##### install SOLVERLAB

```
# SALOME previously installed
echo ${SALOME_ROOT_DIR}   # /volatile/catB/wambeke/salome/SALOME-9.7.0-native-UB20.04-SRC

mkdir ${SOLVERLABDIR}
cd ${SOLVERLABDIR}
git clone --branch master https://github.com/ndjinga/SOLVERLAB.git SOLVERLAB_SRC

# set salome libraries, differents ${USER} (wambeke or yacine-ym268439)
export PETSC_DIR=${SALOME_ROOT_DIR}/INSTALL/petsc
export PETSC_ARCH=""
export MEDFILE_ROOT_DIR=${SALOME_ROOT_DIR}/INSTALL/medfile
export MEDCOUPLING_ROOT_DIR=${SALOME_ROOT_DIR}/INSTALL/MEDCOUPLING
export HDF5_ROOT=${SALOME_ROOT_DIR}/INSTALL/PRODUCTS

# try using conda library Paraview
# Michael doit adapter le CMakeLists.txt pour chercher dans lib/paraview 
# au lieu de paraview/lib et dans include/paraview au lieu de paraview/include
# export PARAVIEW_ROOT_DIR=${ROOTDIR}/miniconda3/envs/py38s97/
# problem as no file paraview.h in install conda

rm -rf ${SOLVERLABDIR}/build
mkdir ${SOLVERLABDIR}/build
cd ${SOLVERLABDIR}/build
pwd  # /volatile/catB/wambeke/solverlab/build

cmake ${SOLVERLABDIR}/SOLVERLAB_SRC \
    -DCMAKE_INSTALL_PREFIX=${SOLVERLABDIR}/install \
    -DCMAKE_BUILD_TYPE=Release  \
    -G"Eclipse CDT4 - Unix Makefiles"  \
    -D_ECLIPSE_VERSION=4.3  \
    -DSOLVERLAB_WITH_DOCUMENTATION=ON  \
    -DPETSC_DIR=${PETSC_DIR}  \
    -DPETSC_ARCH=${PETSC_ARCH}  \
    -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR}  \
    -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR}  \
    -DSOLVERLAB_WITH_GUI=ON

make -j4
make install

# set environment SOLVERLAB
# pluma ${SOLVERLABDIR}/install/env_SOLVERLAB.sh
# bash  # optionally
source ${SOLVERLABDIR}/install/env_SOLVERLAB.sh 

# OLD solverlabGUI try CoreFlows_Standalone.py 
# absent in install 211029 because set obsolete by mickael
# launch directly in SOLVERLAB_SRC if cvw modifications 211029 OK
find ${SOLVERLABDIR} -name "CoreFlows*py"
find ${SOLVERLABDIR} -name "CoreFlows_Standalone.py"
find ${SOLVERLABDIR} -name "CoreFlows*"
tree ${SOLVERLABDIR}/SOLVERLAB_SRC/CoreFlows/gui
${SOLVERLABDIR}/SOLVERLAB_SRC/CoreFlows/gui/CoreFlows_Standalone.py


# test installation -V verbose -I interval de test
cd ${SOLVERLABDIR}/build
# -V is verbosity++ in cas of problem
ctest -V -I 1,1 # test installation -V verbose -I interval de test
ctest -I 1,3

# to check ParaView functionalities
# this is not working with 'ssh -Y'
ctest -V -I 4,5
ctest -V -I 5,5
```


##### Compilation with paraview of salome (NOT VALID YET 211029)

```
# Pour utiliser la lib Paraview de salome si pas dans usr/bin
export PARAVIEW_ROOT_DIR=${SALOME_ROOT_DIR}/BINARIES-UB20.04/ParaView

cmake ${SOLVERLABDIR}/SOLVERLAB_SRC \
    -DCMAKE_INSTALL_PREFIX=${SOLVERLABDIR}/install \
    -DCMAKE_BUILD_TYPE=Release  \
    -G"Eclipse CDT4 - Unix Makefiles"  \
    -D_ECLIPSE_VERSION=4.3  \
    -DSOLVERLAB_WITH_DOCUMENTATION=ON  \
    -DPETSC_DIR=${PETSC_DIR}  \
    -DPETSC_ARCH=${PETSC_ARCH}  \
    -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR}  \
    -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR}  \
    -DPARAVIEW_ROOT_DIR=${PARAVIEW_ROOT_DIR}  \
    -DSOLVERLAB_WITH_GUI=ON
```

#### Usage SOLVERLAB

You have to create a python script from existing examples in 
`${SOLVERLABDIR}/SOLVERLAB_SRC/CoreFlows/examples/Python`


### installation NEW solverlabGUI

Have to install PACKAGESPY, could be outside salome install, or not.

```
# inside salome (as matix do)
export PACKAGESPY_ROOT_DIR=${SALOME_ROOT_DIR}/INSTALL/PACKAGESPY
# outside salome (as standalone)
export PACKAGESPY_ROOT_DIR=${SOLVERLABDIR}/PACKAGESPY

# clone
echo ${PACKAGESPY_ROOT_DIR}
git clone --branch ym_M30 \
    https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git ${PACKAGESPY_ROOT_DIR}
# username wambeke on is240933, mail password
export PYTHONPATH=${PACKAGESPY_ROOT_DIR}/packagespy:${PYTHONPATH}
envs | grep ${SOLVERLABDIR}
```
