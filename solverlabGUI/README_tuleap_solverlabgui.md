
## Install SOLVERLABUI from tuleap salome+matix

```
export MYROOT=${HOME}/PACKAGESPY    # as my choice
mkdir ${MYROOT}
cd ${MYROOT}
git clone --branch cv_M30_2023 ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/matix/packagespy.git

cd packagespy
export PACKAGESPY_ROOT_DIR=$(pwd)
${MYROOT}/packagespy/packagespy/scripts/AllTestLauncherPACKAGESPY.sh

cd ${MYROOT}
git clone --branch bsr_newgui ssh://gitolite@ssh-codev-tuleap.cea.fr:2044/salome/solverlab.git
cd solverlab
export SOLVERLAB_ROOT_DIR=$(pwd)

cd ${MYROOT}
# TODO rewrite part solverlabgui of env_SOLVERLAB.sh

export SOLVERLAB_ROOT_DIR=${SOLVERLAB_ROOT_DIR}/solverlabGUI


export PYTHONPATH=${PACKAGESPY_ROOT_DIR}/packagespy:${PYTHONPATH}
# export PATH=${PACKAGESPY_ROOT_DIR}/packagespy:${PATH}

${SOLVERLAB_ROOT_DIR}/solverlabGUI/solverlabGUI --gui
```


## Create packagespy.tgz

```
cd ${PACKAGESPY_ROOT_DIR}
pwd
./scripts/tgz_for_salome.sh
```


## Usage test SOLVERLABGUI with packagespy.tgz

```
export PACKAGESPY_ROOT_DIR=${MYROOT}/tmp/packagespy
export SOLVERLAB_ROOT_DIR=${MYROOT}/solverlab
export PYTHONPATH=${PACKAGESPY_ROOT_DIR}/packagespy:${PYTHONPATH}

${SOLVERLAB_ROOT_DIR}/solverlabGUI/solverlabGUI --gui
```


