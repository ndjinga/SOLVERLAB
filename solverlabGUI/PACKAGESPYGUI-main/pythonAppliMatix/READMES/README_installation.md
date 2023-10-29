### PACKAGESPY installation


#### Introduction

PACKAGESPY is all-python scripts. With some utilities as bash scripts.


#### Standalone mode

PACKAGESPY needs Python3 PyQt5 etc... environment,
for *example* use [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Download packagespy.git repository where you want

```bash
export P_ROOT_DIR=/export/home/catA/wambeke
cd ${P_ROOT_DIR}
git clone https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git PACKAGESPY
export PACKAGESPY_ROOT_DIR=${P_ROOT_DIR}/PACKAGESPY
```

User have to set PATH etc. for Python3 and PACKAGESPY (if not present)

```bash
export PATH=${PACKAGESPY_ROOT_DIR}/pythonAppliMatix/scripts:${PATH}
export PYTHONPATH=${PACKAGESPY_ROOT_DIR}/pythonAppliMatix:${PYTHONPATH}
```

As all-python scripts PACKAGESPY may be used directly, without any compilation.


#### SALOME mode

Compilation with SAT do the job (with CMake) to get PACKAGESPY in directory SALOME-xx/INSTALL

Example with a SALOME-9.9.0-native-FD34-SRC.tar.gz,
downloaded from [salome-platform.org](https://files.salome-platform.org/Salome)

```bash
export SROOT=/export/home/catA/wambeke/SALOME
export S990=SALOME-9.9.0-native-FD34-SRC

cd ${SROOT}
tar -xf ${S990}.tar.gz

alias ssat=${SROOT}/${S990}/sat/sat  # for example
alias cd_salome="cd ${SROOT}/${S990}"

cd_salome
ssat config --list
export TRG=SALOME-9.9.0-native
```

User have to set PACKAGESPY.pyconf etc. in sat configuration (if not present)

```bash
# TODO This is tricky operation
```

Integration of PACKAGESPY as SALOME module in SALOME-9.9.0 context

```bash
ssat config $TRG --check_system
salome
./install_bin.sh
ssat prepare $TRG -p CONFIGURATION
ssat prepare $TRG -p PACKAGESPY
ssat compile $TRG -p PACKAGESPY --clean_all
ssat launcher $TRG
```
