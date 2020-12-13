CDMATH
======

CDMATH is a geometrical and numerical toolbox designed for numerical analysts who work on the discretisation of partial differential equations on general shapes and meshes and would rather focus on high-level scripting. The library originates from [CDMATH](http://cdmath.jimdo.com), a collaborative workgroup with the same name. It is based on the [MEDcoupling](https://docs.salome-platform.org/latest/dev/MEDCoupling/tutorial/index.html) C++/python library of the [SALOME](http://www.salome-platform.org/) project for the handling of meshes and fields, and on the C++ library [PETSC](https://www.mcs.anl.gov/petsc/) for the handling of matrices and linear solvers. The library is currently developed for linux distributions and is maintained on Ubuntu 14.04 LTS, 16.04 LTS, 18.04 LTSand 20.04 LTS, as well as on Fedora 24, 26, 28, 30 and 32.

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
- [Summary of  available functionalities](Documentation/functionalities.ipynb)
- [Some example scripts](Documentation/examples.md)

Examples of use
---------------
- [Examples of stable numerical methods for the 1D linear transport equation](tests/doc/1DTransportEquation/RegularGrid/TransportEquation1D_RegularGrid.ipynb)
- [Example of stable numerical scheme for the 1D heat equation](tests/doc/1DHeatEquation/HeatEquation1D_RegularGrid.ipynb)
- [Examples of unstable numerical methods for the 1D linear transport equation](tests/doc/1DTransportEquation/UnstableSchemes/TransportEquation1D_UnstableSchemes.ipynb)
- [Shock formation and numerical capture issues for the 1D Burgers' equations](tests/doc/1DBurgersEquation/BurgersEquation1D.ipynb)
- [Numerical convergence and maximum principle analysis for the linear finite element method applied to the 2D Poisson equation](tests/doc/2DPoissonEF/Convergence_Poisson_FE_SQUARE.ipynb)
- [Numerical convergence analysis for the FV5 finite volume method applied to the 2D Poisson equation](tests/doc/2DPoissonVF/Convergence_Poisson_FV5_SQUARE.ipynb)
- [Numerical convergence analysis for the linear finite element and FV5 finite volume methods applied to the 2D Poisson equation with discontinuous boundary conditions](tests/doc/2DPoisson_StiffBC_DISK/2DPoisson_StiffBC_DISK.ipynb)
- [Numerical convergence analysis for the FV5 finite volume method applied to a 2D anisotropic diffusion equation](tests/doc/2DDiffusionVF/Convergence_Diffusion_FV5_SQUARE.ipynb)
- [Influence of the mesh on the convergence and low Mach precision for the UPWIND finite volume method applied to the 2D wave system](tests/doc/2DWaveSystemVF_stationary/Convergence_WaveSystem_Upwind_SQUARE.ipynb)
- [Influence of the mesh on the convergence and low Mach precision  for the CENTERED finite volume method applied to the 2D wave system](tests/doc/2DWaveSystemVF_stationary/Convergence_WaveSystem_Centered_SQUARE.ipynb)
- [Influence of the mesh on the convergence and low Mach precision  for the STAGGERED finite volume method applied to the 2D wave system](tests/doc/2DWaveSystemVF_stationary/Convergence_WaveSystem_Staggered_SQUARE_squares.ipynb)
- [Influence of the mesh on the convergence and low Mach precision  for the PSEUDO-STAGGERED (colocated) finite volume method applied to the 2D wave system](tests/doc/2DWaveSystemVF_stationary/Convergence_WaveSystem_PStag_SQUARE.ipynb)
- [Finite elements for the Poisson problem on a cube in 3D (by S. Kameni Ngwamou, PhD student)](tests/doc/3DPoissonEF/FiniteElements3DPoisson_CUBE.ipynb)
- [Finite elements for the stationary diffusion of the temperature in a 3D room. Influence of the radiator position (by S. Kameni Ngwamou, PhD student)](tests/doc/3DRoomCoolingEF/3DRoomCoolingEF.ipynb)
- [Surface Finite elements for the Poisson-Beltrami problem on a sphere in 3D (by M. Nguemfouo, PhD student)](tests/doc/3DPoissonSphereEF/SynthesisConvergenceFESphere.pdf)
- [Surface Finite elements for the Poisson-Beltrami problem on a torus in 3D (by M. Nguemfouo, PhD student)](tests/doc/3DPoissonTorusEF/SynthesisConvergenceFETorus.pdf)

Download CDMATH sources to compile
----------------------------------

Create your source directory. For instance:
* `mkdir ~/workspace/cdmath`
* `cd ~/workspace/cdmath`

Download from GitHub
* click on the following link : `https://github.com/ndjinga/CDMATH/archive/master.zip`, then unzip the file in a directory cdmath-master
* or type the following in a terminal : `wget https://github.com/ndjinga/CDMATH/archive/master.zip`, then unzip the file in a directory cdmath-master
* or clone the git repository to a folder cdmath-master:  `git clone https://github.com/ndjinga/CDMATH.git cdmath-master`


Set environment for the compilation of CDMATH
---------------------------------------------
Dependencies. The following package list is sufficient on Ubuntu 14.04, Ubuntu 16.04, Ubuntu 18.04 :

 - `cmake3` (mandatory)
 - `g++` or another C++ compiler (mandatory)
 - `python-dev`, `python-numpy` and `swig`, if you want to use CDMATH commands in Python scripts. Use the compilation option `-DCDMATH_WITH_PYTHON=ON`. (highly recommended)
 - `python-matplotlib` and `paraview` for postprocessing tools such as plotting curves (matplotlib) or generating 3D view images (paraview). Use the compilation option `-DCDMATH_WITH_POSTPRO=ON` (recommended).
 - `jupyter`, in order to generate and visualise nice reports from test case simulations (optional)
 - `doxygen`, `graphviz` and `mscgen`, if you want to generate a nice source code documentation in `~/workspace/cdmath/cdmath_install/doc/`. Use the compilation option `-DCDMATH_WITH_DOCUMENTATION=ON` (optional).
 - `libcppunit-dev`, if you want to generate unit tests. Use the compilation option `-DCDMATH_WITH_TESTS=ON` (optional).
 - `rpm`, if you want to generate RPM installation packages. Use the compilation option `-DCDMATH_WITH_PACKAGE=ON` (optional).

Directories. Create the suggested build and installation folders:
* `cd ~/workspace/cdmath`
* `mkdir cdmath_build`
* `mkdir cdmath_install`
* `cd cdmath_build`


Compile and install CDMATH
--------------------------
Simpler build for a minimum version:
* `cmake ../cdmath-master/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PETSC=ON -DCDMATH_WITH_PYTHON=ON `  
> This will download and build the following dependencies
> - PETSc from http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.2.tar.gz
> - SLEPc from https://slepc.upv.es/download/distrib/slepc-3.14.1.tar.gz
> - HDF5 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.3/src/hdf5-1.10.3.tar.gz
> - MEDFILE from http://files.salome-platform.org/Salome/other/med-4.1.0.tar.gz
> - MEDCOUPLING from http://files.salome-platform.org/Salome/other/medCoupling-9.6.0.tar.gz

Advanced build for an all-options version:
* `cmake ../cdmath-master -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DCDMATH_WITH_PETSC=ON -DCDMATH_WITH_PYTHON=ON  -DCDMATH_WITH_POSTPRO=ON -DCDMATH_WITH_TESTS=ON -DCDMATH_WITH_DOCUMENTATION=ON -DPETSC_DIR=${PETSC_DIR} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR}`  
> This assumes that you have an existing 
> - install of PETSc (with submodules SLEPC and HDF5) at the location given by the environment variable PETSC_DIR and the architecture variable PETSC_ARCH  
> See the instructions given in [the official documentation](http://www.mcs.anl.gov/petsc/documentation/installation.html)
> - install of MED                                    at the location given by the environment variable MEDFILE_ROOT_DIR
> - install of MEDCOUPLING                            at the location given by the environment variable MEDCOUPLING_ROOT_DIR

The 3 dependencies PETSC, MED and MEDCOUPLING should have been compiled with the same version of HDF5  
Warning : the linux package libhdf5-dev is generally not compatible with the libraries MED and MEDCoupling
Compile and install:
* `make`
* `make doc install`

Run unit and example tests:
* make example

Run validation tests:
* make validation

Use of CDMATH
-------------
We recommend using CDMATH in python scripts which avoids the hassle of having to edit Makefile files and compiling every C++ scripts.

To use CDMATH with your Python code, you can load the CDMATH environment in your terminal using the command
 * source `~/workspace/cdmath/cdmath_install/env_CDMATH.sh`
Then in your terminal simply type
- `python main.py `

If performance or parallelism is an issue for your simulations, you can use CDMATH librairies with your C++ code :
 * C++ libraries: `export LD_LIBRARY_PATH=~/workspace/cdmath/cdmath_install/lib`
 * To know how to include the right libraries for compilation, see the makefiles of the examples. They include the list ` -lmedC -lmedloader -lmedcoupling -lbase -lmesh -llinearsolver`.

The CDMATH environment variables consist in :
 * CDMATH C++ library path: `~/workspace/cdmath/cdmath_install/lib`
 * CDMATH Python library paths: `~/workspace/cdmath/cdmath_install/lib/cdmath:~/workspace/cdmath/cdmath_install/bin/cdmath`
 * PETSc, SLEPc and HDF5 library path: `${PETSC_DIR}/${PETSC_ARCH}/lib`
 * MED library path: `${MEDFILE_ROOT_DIR}/lib`
 * MEDCOUPLING library path: `${MEDCOUPLING_ROOT_DIR}/lib`

Create Linux installation packages for CDMATH
---------------------------------------------
After popular request, here is how you can create packages for Debian and Red Hat-based Linux distributions:

1. Download CDMATH as explained hereabove.
2. Set the environment as explained hereabove (in particular, make sure you have `rpm` installed).
3. Generate a makefile with `cmake -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PACKAGE=ON ../cdmath-master/` and eventually other options (documentation, tests, swig, etc).
4. Compile with `make package`.

You will then find a Debian package in the build directory; you may install it on Ubuntu 14.04 or Ubuntu 16.04. You will also find an RPM package, which you may install on Red Hat-based distributions. This way, the packages you generate may include all the compilation options you want.

Unfortunately, the Debian package may be said to be of “bad quality” for Debian standards as far as ownership is concerned. This is true and due to limitations in CMake/CPack. The package should still install nonetheless.
