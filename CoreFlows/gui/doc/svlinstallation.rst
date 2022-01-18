
.. include:: ./rst_prolog.rst


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


Standalone installation
--------------------------

The most simple installation process and should suit to most users.

To use the GUI you need to have installed Python3_ on your computer with the following dependency:

- PyQt5_

- TODO packagespy

.. warning:: For the moment Solverlab **doesn't** work with Conda TODO why ? can't find some librairy

When all the dependency are installed on your computer you can follow the `Solverlab install process <Solveralabinstall_>`_

When it's done you should be able to launch the GUI with the command

.. code-block:: bash

   python3 $SOLVERLABGUI


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


