CDMATH-CoreFlows
================

CDMATH-CoreFlows is an open source C++/Python library intended at solving PDE systems
arising from the thermalhydraulics of two phase flows in power plant boilers. It
is a simple environment meant at students and researchers to test new numerical
methods on general geometries with unstructured meshes. It is developped by
CEA Saclay since 2014 and proposes a few
basic models and finite volume numerical methods. Some of the main objectives
are the study of

- Numerical schemes for compressible flows at low Mach numbers
- Well balanced schemes for stiff source terms (heat source, phase change, pressure losses)
- Flow inversion and counter-current two phase flows
- Schemes that preserve the phasic volume fraction α ∈ [0, 1]
- Convergence of finite volume methods
- New preconditioners for implicit methods for two phase flows
- The coupling of fluid models or multiphysics coupling (eg thermal hydraulics and neutronics or thermal hydraulics and solid thermics)

CDMATH-CoreFlows relies on the numerical toolbox [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) of the project [CDMATH](http://cdmath.jimdo.com) for the handling of meshes and fields, and on the library [PETSC](https://www.mcs.anl.gov/petsc/) for the handling of large sparse matrices.
You will need the packages 'doxygen' if you want to generate de documentation and 'swig' if you want to use python scripts.  The software is currently developed for linux distributions and is maintained on Ubuntu 14.04 LTS, 16.04 LTS and 18.04 LTS, as well as on Fedora 24, 26 and 29.

User guide
----------
The user guide is organized as follows :
- [The physical models](./Documentation/PhysicalModels.md)
    - [The linear scalar problems](./Documentation/PhysicalModels/ScalarModelsPage.ipynb)
        - [The transport equation](./Documentation/PhysicalModels/TransportEq.ipynb) for pure advection phenomena
        - [The diffusion equation](./Documentation/PhysicalModels/DiffusionEq.ipynb) for pure diffusion phenomena
    - [The compressible Navier-Stokes equations](./Documentation/PhysicalModels/NSModelsPage.ipynb)
    - [The two-phase flow models](./Documentation/PhysicalModels/TwoPhasePage.ipynb)
        - [The drift model](./Documentation/PhysicalModels/TwoPhase/DriftModelPage.ipynb) with two partial masses, one momentum and one energy equation
        - [The isothermal two-fluid model](./Documentation/PhysicalModels/TwoPhase/IsothermalPage.ipynb) with two partial masses and two momentum equations (no energy equation)
        - [The five equation two-fluid model](./Documentation/PhysicalModels/TwoPhase/FiveEqPage.ipynb) with two partial masses, two momentum equations and one energy equation
- [Software structure](Documentation/software.md)
- [The numerical methods](Documentation/numericalPage.ipynb)
- [Summary of  available functionalities](Documentation/functionalities.ipynb)
- [CDMATH-CoreFlows example scripts](Documentation/examples.md)

Download and compilation of CDMATH and PETSc
--------------------------------------------
[CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) can be downloaded and compiled together with [PETSC](https://www.mcs.anl.gov/petsc/) in a single process, thanks to the cmake option -DCDMATH_WITH_PETSC=ON.
We summarise the installation procedure that you can find in detailed form here
- https://github.com/ndjinga/CDMATH

In order to compile [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) you will need the packages 'cmake', 'gcc', 'gfortran', 'hdf5' plus 'numpy' and 'swig'.
In a linux terminal, first create and access a working directory :
- `mkdir -p ~/workspace/cdmath `
- `cd ~/workspace/cdmath `

Then create build and install repositories:
- `mkdir cdmath_build cdmath_install `

In order to download the approriate branch of [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) either unzip the following file to a directory cdmath-master
- `https://github.com/ndjinga/CDMATH/archive/master.zip`
or clone the git repository to a folder cdmath-master
- `git clone https://github.com/ndjinga/CDMATH.git cdmath-master`

This latter command results in the creation of a directory `~/workspace/cdmath/cdmath-master` containing the source files of [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH).

Go to the build directory
- `cd cdmath_build `

Then run the commands
- `cmake ../cdmath-master/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PYTHON=ON -DCDMATH_WITH_PETSC=ON`
- `make`
- `make install`

By default, [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) will compile a new sequential installation of [PETSC](https://www.mcs.anl.gov/petsc/). If an installation of [PETSC](https://www.mcs.anl.gov/petsc/) (version 3.4 or later) is already available in the system, it is possible to save time by first setting the environment variables PETSC_DIR and PETSC_ARCH to the appropriate values as can be found in petscconf.h, and then running the above cmake command.

Download and compilation of CoreFlows
---------------------------------------------
First create and access a working directory :
- `mkdir -p ~/workspace/CDMATH-CoreFlows `
- `cd ~/workspace/CDMATH-CoreFlows `
Now create build and install repositories:
- `mkdir CDMATH-CoreFlows_build CDMATH-CoreFlows_install `

In order to download CDMATH-CoreFlows either unzip the following file to a directory CDMATH-CoreFlows-master
- `https://github.com/ndjinga/CDMATH-CoreFlows/archive/master.zip`
or clone the git repository to a folder CDMATH-CoreFlows-master
- `git clone https://github.com/ndjinga/CDMATH-CoreFlows.git CDMATH-CoreFlows-master`
Either of these latter commands results in the creation of a directory `~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows-master`  containing the source files.

In the following steps we assume that [PETSC](https://www.mcs.anl.gov/petsc/) (version 3.4 or more recent) has been installed with CDMATH with the process described above.
You need to set the following variables 
- `CDMATH_INSTALL`, the path to your CDMATH installation, for example  `~/workspace/cdmath/cdmath_install//share/petsc-3.8.3 `
- `PETSC_DIR`, the path to your PETSc installation. If [PETSC](https://www.mcs.anl.gov/petsc/) was installed by CDMATH then [CDMATH-Toolbox](https://github.com/ndjinga/CDMATH) can be defined as `~/workspace/cdmath/cdmath_install`
- `PETSC_ARCH`, the type of installation used (usually arch-linux2-c-opt or linux-gnu-c-opt)

In order to do so, it is sufficient to source the 'CDMATH' environment file. Type in you linux terminal
- `source ~/workspace/cdmath/cdmath_install/env_CDMATH.sh`

then create build and install repositories next to CoreFlows-master :
- `mkdir CoreFlows_build CoreFlows_install`

Go to the build directory
- `cd CoreFlows_build `

Then run the command
- `../CDMATH-CoreFlows-master/configure  --prefix=../CDMATH-CoreFlows_install/ --with-petsc-dir=$PETSC_DIR --with-petsc-arch=$PETSC_ARCH --with-cdmath-dir=$CDMATH_INSTALL --with-python --with-doc`
- `make doc install`

You can add the following optional commands
- `--with-gui`, if you want to use CDMATH-CoreFlows as a Salomé module (you will need to use a Salomé shell)
- `--with-debug`, if you want to use CDMATH-CoreFlows in debug mode instead of the default optimised mode

Use of CDMATH-CoreFlows
-----------------------
First load CDMATH-CoreFlows environment from the CoreFlows-master directory
- `source ~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_install/env_CoreFlows.sh `

If you use C language: edit the file CoreFlows-master/CoreFlows_src/main.cxx then in a terminal type
- `cd ~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_build  `
- `make`
- `make install`
Then you can run the simulation in any directory with the command line
- `$CoreFlows `

If you use python language: edit your own python file `my_file.py` following for example the pattern of the file `CDMATH-CoreFlows-master/main.py`. Then in a terminal type
- `python my_file.py `

If you use the graphic interface, you need to run a [SALOME](https://www.salome-platform.org/) Unix shell 
- `./salome shell`
and type the command line
- `runSalome -mCOREFLOWS`
then click on new study to open CoreFlows interface

The complete documentation is available in the directory `~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_install/share/doc/`
