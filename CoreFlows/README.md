SOLVERLAB-CoreFlows
================

SOLVERLAB-CoreFlows is an open source C++/Python library intended at solving PDE systems
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

SOLVERLAB-CoreFlows relies on the numerical toolbox [SOLVERLAB-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) originating from the project [CDMATH](http://cdmath.jimdo.com) for the handling of meshes and fields, and on the library [PETSC](https://petsc.org/release/) for the handling of large sparse matrices.
You will need the packages 'doxygen' and 'sphinx' if you want to generate de documentation and 'swig' if you want to use python scripts.  
The library is currently maintained and distributed by the SALOME developpement team on various linux distributions (Ubuntu, CentOS, Fedora, Debian) and on Windows-10.

User guide of the CoreFlows module
----------------------------------
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
- [SOLVERLAB-CoreFlows example scripts](Documentation/examples.md)

Download and compilation of CDMATH and PETSc
--------------------------------------------
[CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) can be downloaded and compiled together with [PETSC](https://petsc.org/release/) in a single process, thanks to the cmake option -DCDMATH_WITH_PETSC=ON.
We summarise the installation procedure that you can find in detailed form here
- https://github.com/ndjinga/CDMATH

In order to compile [CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) you will need the packages 'cmake', 'gcc', 'gfortran', 'hdf5' plus 'numpy' and 'swig'.
In a linux terminal, first create and access a working directory :
- `mkdir -p ~/workspace/cdmath `
- `cd ~/workspace/cdmath `

Then create build and install repositories:
- `mkdir cdmath_build cdmath_install `

In order to download the approriate branch of [CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) either unzip the following file to a directory cdmath-master
- `https://github.com/ndjinga/CDMATH/archive/master.zip`
or clone the git repository to a folder cdmath-master
- `git clone https://github.com/ndjinga/CDMATH.git cdmath-master`

This latter command results in the creation of a directory `~/workspace/cdmath/cdmath-master` containing the source files of [CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH).

Go to the build directory
- `cd cdmath_build `

Then run the commands
- `cmake ../cdmath-master/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PYTHON=ON -DCDMATH_WITH_PETSC=ON`
- `make`
- `make install`

By default, [CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) will compile a new sequential installation of [PETSC](https://petsc.org/release/). If an installation of [PETSC](https://petsc.org/release/) (version 3.4 or later) is already available in the system, it is possible to save time by first setting the environment variables PETSC_DIR and PETSC_ARCH to the appropriate values as can be found in petscconf.h, and then running the above cmake command.

Download and compilation of CoreFlows
---------------------------------------------
First create and access a working directory :
- `mkdir -p ~/workspace/SOLVERLAB-CoreFlows `
- `cd ~/workspace/SOLVERLAB-CoreFlows `
Now create build and install repositories:
- `mkdir SOLVERLAB-CoreFlows_build SOLVERLAB-CoreFlows_install `

In order to download SOLVERLAB-CoreFlows either unzip the following file to a directory SOLVERLAB-CoreFlows-master
- `https://github.com/ndjinga/SOLVERLAB-CoreFlows/archive/master.zip`
or clone the git repository to a folder SOLVERLAB-CoreFlows-master
- `git clone https://github.com/ndjinga/SOLVERLAB/CoreFlows/CoreFlows.git SOLVERLAB-CoreFlows-master`
Either of these latter commands results in the creation of a directory `~/workspace/SOLVERLAB-CoreFlows/SOLVERLAB-CoreFlows-master`  containing the source files.

In the following steps we assume that [PETSC](https://petsc.org/release/) (version 3.4 or more recent) has been installed with CDMATH with the process described above.
You need to set the following variables 
- `CDMATH_INSTALL`, the path to your CDMATH installation, for example  `~/workspace/cdmath/cdmath_install/`
- `PETSC_DIR`, the path to your PETSc installation. If [PETSC](https://petsc.org/release/) was installed by CDMATH then [CDMATH-Toolbox](https://github.com/ndjinga/SOLVERLAB/tree/master/CDMATH) can be defined as `~/workspace/cdmath/cdmath_install`
- `PETSC_ARCH`, the type of installation used (usually arch-linux2-c-opt or linux-gnu-c-opt)

In order to do so, it is sufficient to source the 'CDMATH' environment file. Type in you linux terminal
- `source ~/workspace/cdmath/cdmath_install/env_CDMATH.sh`

then create build and install repositories next to CoreFlows-master :
- `mkdir CoreFlows_build CoreFlows_install`

Go to the build directory
- `cd CoreFlows_build `

Then run the command
- `../SOLVERLAB-CoreFlows-master/configure  --prefix=../SOLVERLAB-CoreFlows_install/ --with-petsc-dir=$PETSC_DIR --with-petsc-arch=$PETSC_ARCH --with-cdmath-dir=$CDMATH_INSTALL --with-python --with-doc`
- `make doc install`

You can add the following optional commands
- `--with-gui`, if you want to use SOLVERLAB-CoreFlows as a Salomé module (you will need to use a Salomé shell)
- `--with-debug`, if you want to use SOLVERLAB-CoreFlows in debug mode instead of the default optimised mode

Use of SOLVERLAB-CoreFlows
-----------------------
First load SOLVERLAB-CoreFlows environment from the CoreFlows-master directory
- `source ~/workspace/SOLVERLAB-CoreFlows/SOLVERLAB-CoreFlows_install/env_CoreFlows.sh `

If you use C language: edit the file CoreFlows-master/CoreFlows_src/main.cxx then in a terminal type
- `cd ~/workspace/SOLVERLAB-CoreFlows/SOLVERLAB-CoreFlows_build  `
- `make`
- `make install`
Then you can run the simulation in any directory with the command line
- `$CoreFlows `

If you use python language: edit your own python file `my_file.py` following for example the pattern of the file `SOLVERLAB-CoreFlows-master/main.py`. Then in a terminal type
- `python my_file.py `

If you use the graphic interface, you need to run a [SALOME](https://www.salome-platform.org/) Unix shell 
- `./salome shell`
and type the command line
- `runSalome -mCOREFLOWS`
then click on new study to open CoreFlows interface

The complete documentation is available in the directory `~/workspace/SOLVERLAB-CoreFlows/SOLVERLAB-CoreFlows_install/share/doc/`
