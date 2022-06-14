
.. include:: ../rst_prolog.rst


.. _svlInstallation:

**************
Installation
**************

The easiest way to use SOLVERLAB_ is to download the SALOME_ targz file correponding to your operating system (OS), and to use its Python_ shell environment.
However SALOME_ is a much bigger project than SOLVERLAB_, binaries for your OS may not be avalaible, and one way which to build SOLVERLAB_ as a light standalone library.
This section describes the build of SOLVERLAB_ from its source files.

Prerequisite
--------------

Solverlab in his core componment is a really small package and need other librairies to work properly.

Some dependencies are needed for the installation and need to be on your computer before trying to install Solverlab:

- cmake3
- g++ (or another C++ compiler)
- python3-dev, python3-numpy, swig3, python3-pandas
- PyQt5_

.. _Simple-standalone-installation:

Simple standalone installation
------------------------------

The most simple installation process and should suit to most users.
This way of installation will download the some dependencies, you need to have an internet onnection.

Suppose the source files are located in `~/workspace/SOLVERLAB/SOLVERLAB-master`.  
Start by creating two folders on the same level as the source.

.. code-block:: bash

    cd ~/workspace/SOLVERLAB    
    mkdir SOLVERLAB_build  
    mkdir SOLVERLAB_install

Then compile the project with cmake in the build folder.

.. code-block:: bash
    
    cd SOLVERLAB_build
    cmake ../SOLVERLAB-master/ -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -DSOLVERLAB_WITH_DOCUMENTATION=ON -DSOLVERLAB_WITH_GUI=ON -DSOLVERLAB_WITH_PACKAGESPY=ON
    make

To generate the user GUI_ documentation use 

.. code-block:: bash

    make docGUI
    
To install the project and its GUI_ documentation use 

.. code-block:: bash

    make install
    

Before any use of Solverlab don't forget to set its environement variable. 
You will find the script env_SOLVERLAB.sh in the install folder.

.. code-block:: bash
    
    source ~/workspace/SOLVERLAB/SOLVERLAB_install/env_SOLVERLAB.sh

You can now launch the GUI with

.. code-block:: bash

   python3 $SOLVERLABGUI -g


Advanced Installation
------------------------

If you already have PETSC, MED and MEDCoupling installed on your system, you can specify their path and the installation process will save time and disk space by not downloading them again.

Let us assume that you have an existing installation of:

- PETSC (with submodules SLEPC and HDF5) at the location given by the environment variable PETSC_DIR and the architecture variable PETSC_ARCH

- MED at the location given by the environment variable MEDFILE_ROOT_DIR

- MEDCoupling at the location given by the environment variable MEDCOUPLING_ROOT_DIR

The 3 dependencies PETSC, MED and MEDCOUPLING should have been compiled with the same version of HDF5.

You also need a copy of the PACKAGESPY python library in a folder that we will call PACKAGESPY_ROOT_DIR.

Now you can use an advanced cmake configuration :

.. code-block:: bash 
    
    cmake ../SOLVERLAB-master -DCMAKE_INSTALL_PREFIX=../SOLVERLAB_install -DCMAKE_BUILD_TYPE=Release -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DSOLVERLAB_WITH_DOCUMENTATION=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH=${PETSC_ARCH} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ROOT_DIR=${MEDCOUPLING_ROOT_DIR} -DSOLVERLAB_WITH_GUI=ON -DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}

The next steps of the installation are similar as those of the :ref:Simple-standalone-installation above.

.. code-block:: bash

   make  
   make docGUI  
   make install

You can now launch the GUI with

.. code-block:: bash

   source ~/workspace/SOLVERLAB/SOLVERLAB_install/env_SOLVERLAB.sh  
   python3 $SOLVERLABGUI -g

