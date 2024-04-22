## compilation SOLVERLAB in SALOME 990 with sat

- From a salome.targz installation
- From SOLVERLAB github
- From PACKAGESPY tuleap SALOME-MATIX


### Compilation sat

```
export SROOT=/export/home/catA/wambeke/SALOME
export S990=SALOME-9.9.0-native-FD34-SRC
alias cd_salome="cd ${SROOT}/${S990}"

cd ${SROOT}
export SAT=${SROOT}/${S990}/sat/sat
$SAT config --list
export TRG=SALOME-9.9.0-native

cd_salome
sat/sat config $TRG --check_system
salome

# prepare sat compile
./install_bin.sh # prepare sat compile

# get references

cd ${SROOT}
git clone https://codev-tuleap.intra.cea.fr/plugins/git/matix/SAT_MATIX.git
cd SAT_MATIX
git checkout M30

cd ${SROOT}
git clone https://codev-tuleap.cea.fr/plugins/git/salome/sat_salome.git
cd sat_salome
git checkout cv_solverlab_M30

cd ${SROOT}
git clone https://codev-tuleap.cea.fr/plugins/git/salome/sat.git
cd sat
git checkout master


cd ${SROOT}
$SAT config $TRG -p SOLVERLAB

$SAT config $TRG -g .| grep "SOLVERLAB.pyconf"
find . -name "SOLVERLAB*"
  PROJECT/products/SOLVERLAB.pyconf
  PROJECT/products/env_scripts/SOLVERLAB.py
  PROJECT/products/compil_scripts/SOLVERLAB.sh

cd ${SROOT}
find . -name SOLVERLAB.pyconf
cp ./SALOME-9.9.0-native-FD34-SRC/PROJECT/products/SOLVERLAB.pyconf \
   ./SALOME-9.9.0-native-FD34-SRC/PROJECT/products/SOLVERLAB.pyconf_ORI
meld $(find . -name "SOLVERLAB.pyconf") &
meld $(find . -name "SOLVERLAB.py") &
meld $(find . -name "SOLVERLAB.sh") &

# avoid "SOURCES of CONFIGURATION not found!""
$SAT prepare $TRG -p CONFIGURATION

find . -name "*.pyconf" -exec fgrep -Hni SOLVERLAB {} \; | grep tag
  ...etc...
  ./SAT_MATIX/applications/RELEASED_MATIX_V30_970_CV_COS73.pyconf:170:   
            'SOLVERLAB' : {section : 'default_github', tag : 'yacine_gui'}

find . -name "*.pyconf" -exec fgrep -Hni SOLVERLAB {} \; | grep $TRG

cd SAT_MATIX
git remote -v
  origin https://codev-tuleap.intra.cea.fr/plugins/git/matix/SAT_MATIX.git (push)


pluma ~/.salomeTools/SAT.pyconf  # editor->pluma
$SAT config $TRG -e  # -> SOLVERLAB : {section : 'default_github', tag : 'yacine_gui'}

$SAT prepare $TRG -p PACKAGESPY
$SAT compile $TRG -p PACKAGESPY --clean_all

$SAT prepare $TRG -p SOLVERLAB
$SAT compile $TRG -p SOLVERLAB --clean_all


cp ./SAT_MATIX/products/PACKAGESPY.pyconf /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/PROJECT/products/.

# compilation log
$SAT log $TRG
cat ${SROOT}/${S990}/LOGS/SOLVERLAB/* | less -r
cat ${SROOT}/${S990}/LOGS/PACKAGESPY/* | less -r

mv mesa_salome mesa_salome_ORI
mv salome salome_ORI

$SAT launcher $TRG

# try salome with SOLVERLAB
$SROOT/$S990/salome
```

### Configuration SOLVERLAB sat

```
$SAT config $TRG -p SOLVERLAB
  SOLVERLAB is a product
    depends on = GUI,          KERNEL,     MEDCOUPLING,        ParaView,          Python,
                hdf5,      matplotlib,         medfile,           numpy,           petsc,

    build depend on = cmake,         cppunit,         doxygen,        graphviz,            swig,


  configuration:
    pyconf file path = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/PROJECT/products/SOLVERLAB.pyconf
    section = default

  prepare:
    get method = archive
    get from = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/ARCHIVES/SOLVERLAB.tgz

  compile:
    compilation method = script
    Compilation script = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/PROJECT/products/compil_scripts/SOLVERLAB.sh
    make -j = 1
    source dir = '/export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/SOURCES/SOLVERLAB' ** not found
    build dir = '/export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/BUILD/SOLVERLAB' ** not found
    install dir = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB
    debug  = no
    verbose  = no
    hpc  = no
    dev  = no

  environ :
    script = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/PROJECT/products/env_scripts/SOLVERLAB.py
    set          PYTHON_LIBDIR = lib/python${PYTHON_VERSION}/site-packages
    set          SOLVERLAB_ROOT_DIR = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB
    set          SOLVERLAB_SRC_DIR = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/SOURCES/SOLVERLAB
    prepend      PATH = ${SOLVERLAB_ROOT_DIR}/bin/salome:${PATH}
    prepend      LD_LIBRARY_PATH = ${SOLVERLAB_ROOT_DIR}/lib/salome:${LD_LIBRARY_PATH}
    prepend      PYTHONPATH = ${SOLVERLAB_ROOT_DIR}/bin/salome:${SOLVERLAB_ROOT_DIR}/lib/salome:${SOLVERLAB_ROOT_DIR}/${PYTHON_LIBDIR}/salome:${PYTHONPATH}
    append       SALOME_MODULES = ${SALOME_MODULES},SOLVERLAB
    append       SalomeAppConfig = ${SalomeAppConfig}:/export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB/share/salome/resources/solverlab
    set          CoreFlows_INSTALL = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB
    set          CoreFlows_ROOT_DIR = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB
    set          CoreFlows_ROOT = /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/SOLVERLAB
    set          CoreFlows_PYTHON = ON
    set          CoreFlows_DOC = ON
    set          CoreFlows_GUI = ON
    set          CoreFlows = ${CoreFlows_INSTALL}/bin/CoreFlowsMainExe
    set          CoreFlowsGUI = ${CoreFlows_INSTALL}/bin/CoreFlows_Standalone.py
    set          COREFLOWS_ROOT_DIR = ${CoreFlows_ROOT_DIR}
    prepend      PATH = ${CoreFlows_ROOT_DIR}/include:${PATH}
    prepend      LD_LIBRARY_PATH = ${CoreFlows_ROOT_DIR}/lib:${LD_LIBRARY_PATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/lib:${PYTHONPATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/lib/coreflows:${PYTHONPATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/bin/coreflows:${PYTHONPATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/lib/cdmath:${PYTHONPATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/bin/cdmath:${PYTHONPATH}
    prepend      PYTHONPATH = ${CoreFlows_ROOT_DIR}/bin/cdmath/postprocessing:${PYTHONPATH}
```
