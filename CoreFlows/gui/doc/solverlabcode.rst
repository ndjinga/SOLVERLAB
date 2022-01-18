
.. include:: ./rst_prolog.rst

.. _iraSolverlabCode:

***************************
Solverlab code compilation
***************************

SolverlabGUI uses a specific version of Solverlab code, modified by J.P. Crocombette (cea),
which is tagged version 1.1.x for now.
This code comes from original version 1.0.8 by Christian Borschel.

Users find two current compiled executable files, which are used by solverlabGUI, located at
*.../solverlabGUI/solverlabCode*.

#. *solverlab_mingw64.exe*, compiled by MinGW_ gcc compiler, for Windows (64 bits).
#. *solverlab_linux64.exe*, compiled by GNU_ `gcc <GNUgcc_>`_ compiler, for Linux (64 bits).

A development installation of Solverlab code allows programmer improvments.
The following chapters explain Solverlab code compilation processes.


.. _iraSolverlabCodeSources:

Solverlab code sources
============================

With GPL licence, sources are available in solverlabGUI directory tree, located at
*.../solverlabGUI/solverlabCodes/solverlab_cea*.

User find also the useful Corteo_ data base, located at
*.../solverlabGUI/solverlabCodes/data_4bit*

.. code-block:: bash

    .../solverlabCodes > tree
    .
    ├── data_4bit
    │   ├── 10.asp
    │   ├── 11.asp
    etc.
    │   ├── corteo.mat
    │   └── erfinv.dat
    ├── doc
    │   ├── 20140804_solverlab_manual.pdf
    │   ├── Corteo20160816.pdf
    │   ├── solverlab-1-s2.0-S0168583X11006318-main.pdf
    │   └── Solverlab_tuto_installation.pdf
    ├── solverlab_cea
    │   ├── compileSolverlab.bat
    │   ├── fileio.c
    │   ├── fileio.h
    │   ├── fromcorteo.c
    │   ├── fromcorteo.h
    │   ├── geometry.c
    │   ├── geometry.h
    │   ├── indexvalues6bit.h
    │   ├── indexvalues.h
    │   ├── solverlab.c
    │   ├── solverlab.h
    │   ├── license.txt
    │   ├── makefile_cea
    │   ├── target.c
    │   ├── target.h
    │   ├── transport.c
    │   ├── transport.h
    │   ├── utils.c
    │   └── utils.h
    ├── compileSolverlab.lnk
    └── README.txt



.. _iraSolverlabCodeCompilation_linux:

Solverlab code compilation Linux
==================================

Example of compilation *(Linux-bash)*:

.. code-block:: bash

    # this is your location
    cd .../solverlabGUI/solverlabCodes/solverlab_cea
    # verifications
    cat README.txt
    # compilation
    make -f makefile_cea clean
    make -f makefile_cea solverlab
    make -f makefile_cea installGUI  # install executable in solverlabGUI/solverlabCode



.. _iraSolverlabCodeCompilation_windows:

Solverlab code compilation Windows7-10
=======================================

.. warning:: #. MinGW_ is supposed to be set
                and useful in environment path, at an usual location
                *C:\\MinGW* (for example).
             #. Git-windows_ is supposed to be set
                and useful in environment path, at an usual location
                *C:\\Program Files\\Git* (for example).
                In order to use like-Linux commands.


Example of compilation *(Windows7/10-cmd.exe shell)*:

.. code-block:: bat

    # this is mandatory location
    C:\
    cd C:\Users\Public\solverlab\solverlabGUI\solverlabCodes\solverlab_cea
    # verifications
    where make                       # --> C:\MinGW\bin\make.exe
    where gcc                        # --> C:\MinGW\bin\gcc.exe
    where uname                      # --> C:\Program Files\Git\usr\bin\uname.exe
    # compilation
    make -f makefile_cea clean
    make -f makefile_cea solverlab
    make -f makefile_cea installGUI  # install executable in solverlabGUI/solverlabCode


.. note:: To launch Solverlab code compilation, you may use Windows shortcut
          *C:\\User\\Public\\solverlab\\solverlabGUI\\solverlabCodes\\compileSolverlab(.lnk)*
          .




