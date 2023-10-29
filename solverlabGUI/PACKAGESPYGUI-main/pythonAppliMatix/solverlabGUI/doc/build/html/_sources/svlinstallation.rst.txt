
.. include:: ../rst_prolog.rst


.. _svlInstallation:

**************
Installation
**************

In a normal use case Solverlab should be installed with SALOME (TODO link salome)
This section describe other method to get it.

Prerequisite
--------------

Solverlab in his core componment is a really small package and need other librairies to work properly.

Some dependencies are needed for the installation and need to be on your computer before trying to install Solverlab:

- cmake3
- g++ (or another C++ compiler)
- python3-dev, python3-numpy, swig3
- PyQt5_


Standalone installation
--------------------------

The most simple installation process and should suit to most users.
This way of installation will download the some dependencies, you need to have internet.

Start by creating two folders on the same level as the source.

.. code-block:: bash

    cd ~/workspace/SOLVERLAB    
    mkdir SOLVERLAB_build
    mkdir SOLVERLAB_install

Then we are going to compile the project with cmake in the build folder and install it.

.. code-block:: bash
    
    cd SOLVERLAB_build
    cmake ../SOLVERLAB-master/ -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -DSOLVERLAB_WITH_GUI=ON -DSOLVERLAB_WITH_DOCUMENTATION=ON
    make
    make install

To check if your installation use 

.. code-block:: bash

    make validation
    
If all test succed you can now go in the install folder.

Before any use of Solverlab don't forget to set its environement variable. 
You will find the script env_SOLVERLAB.sh in the install folder.

.. code-block:: bash
    
    source ~/workspace/SOLVERLAB/SOLVERLAB_install/env_SOLVERLAB.sh

You can now launch the GUI with

.. code-block:: bash

   $SOLVERLABGUI


Advanced Installation
------------------------

If you already have PETSC MED or MEDCoupling on your device you can specify their path and the installation process will not download them again.

This assumes that you have an existing install of:

- PETSC (with submodules SLEPC and HDF5) at the location given by the environment variable PETSC_DIR and the architecture variable PETSC_ARCH

- MED at the location given by the environment variable MEDFILE_ROOT_DIR

- MEDCoupling at the location given by the environment variable MEDCOUPLING_ROOT_DIR

The 3 dependencies PETSC, MED and MEDCOUPLING should have been compiled with the same version of HDF5

.. code-block:: bash 
    
    cmake ../SOLVERLAB-master -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DSOLVERLAB_WITH_DOCUMENTATION=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH=${PETSC_ARCH} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR} -DSOLVERLAB_WITH_GUI=ON


