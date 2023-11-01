Standalone compilation of SOLVERLAB from sources
================================================

Prerequisites for the compilation of SOLVERLAB
---------------------------------------------
The following package list is sufficient on Ubuntu 20.04 :

 - `cmake3` (mandatory)
 - `g++` or another C++ compiler (mandatory)
 - `python3-dev`, `python3-numpy` and `swig3`for python scripts (mandatory)
 - `pyqt5-dev-tools` to generate the Graphical User Interface (optional)
 - `python3-matplotlib`, `paraview-dev`, `libnetcdf-dev` (on Ubuntu 20.04) and `python3-paraview` for postprocessing tools such as plotting curves (matplotlib) or generating 3D view images (paraview) (optional)
 - `ffmpeg` and `ffmpeg-devel` to generate an animation from a set of curves (optional)
 - `libmpich-dev` or  `libopenmpi-dev` for mpi parallelisation (optional)
 - `libcppunit-dev`, if you want to generate unit tests. Use the compilation option `-DSOLVERLAB_WITH_TESTS=ON` (optional).
 - `rpm`, if you want to generate RPM installation packages. Use the compilation option `-DSOLVERLAB_WITH_PACKAGE=ON` (optional).


Download SOLVERLAB sources
--------------------------
Download SOLVERLAB sources from GitHub
* use the following url in a browser : `https://github.com/ndjinga/SOLVERLAB/archive/master.zip`, then download and unzip the file in a directory SOLVERLAB-master
* or type the following in a terminal : `wget https://github.com/ndjinga/SOLVERLAB/archive/master.zip`, then unzip the file in a directory SOLVERLAB-master
* or clone the git repository to a folder SOLVERLAB-master:  `git clone https://github.com/ndjinga/SOLVERLAB.git SOLVERLAB-master`


Compile and install SOLVERLAB from source files
-----------------------------------------------
First create a directory named 'build' where the compilation will take place and open a terminal in that directory.

**Simple configuration for a minimum version:**  
Type the following cmake instruction in a terminal for a minimal version
* `cmake /path/to/SOLVERLAB-master/ -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -DSOLVERLAB_WITH_DOCUMENTATION=ON `  
> This will download and build the following dependencies
> - PETSc from https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz
> - SLEPc from https://slepc.upv.es/download/distrib/slepc-3.17.0.tar.gz
> - F2CBLASLAPACK from https://ftp.mcs.anl.gov/pub/petsc/externalpackages/f2cblaslapack-3.4.2.q4.tar.gz
> - HDF5 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.3/src/hdf5-1.10.3.tar.gz
> - MEDFILE from http://files.salome-platform.org/Salome/other/med-4.1.1.tar.gz
> - MEDCOUPLING from http://files.salome-platform.org/Salome/other/medCoupling-9.8.0.tar.gz

**Advanced configuration for a complete version:**  
If you already have an installation of PETSC, MED and MEDCoupling, you may save computational time and memory by typing the following cmake instruction in a terminal :
* `cmake /path/to/SOLVERLAB-master -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DSOLVERLAB_WITH_DOCUMENTATION=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH=${PETSC_ARCH} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR} `  
> This assumes that you have an existing 
> - installation of PETSc (with submodule SLEPC) at the location given by the environment variable PETSC_DIR and the architecture variable PETSC_ARCH  
> See the instructions given in [the official documentation](https://petsc.org/release/install/)
> - installation of MED                                    at the location given by the environment variable MEDFILE_ROOT_DIR
> - installation of MEDCOUPLING                            at the location given by the environment variable MEDCOUPLING_ROOT_DIR

The 2 dependencies MED and MEDCOUPLING should have been compiled with the same version of HDF5  
Warning : the linux package libhdf5-dev is generally not compatible with the libraries MED and MEDCoupling  

**Compile and install:**
* `make`
* `make install`

Run tests and documentation
---------------------------
Run unit and example tests:
* `make examples`

Run validation tests:
* `make validation`

Generate html user guide 
* `make doc-user`

Generate html developer guide
* `make doc`

Parallel version
---------------------
To compile a parallel version of SOLVERLAB, you first need parallel C and C++ compilers and a parallel installation of the packages PETSc, MED and MEDCoupling.  

Then you need to add the following extra CMAKE parameters to the configuration step :  
* ` -DMPI_HOME=${MPI_HOME} -DCMAKE_C_COMPILER=${MPI_C_COMPILER} -DCMAKE_CXX_COMPILER=${MPI_CXX_COMPILER} `

where
- MPI_HOME stores the location of the mpi library
- MPI_C_COMPILER stores the location of the parallel C compiler
- MPI_CXX_COMPILER stores the location of the parallel C++ compiler

*/
