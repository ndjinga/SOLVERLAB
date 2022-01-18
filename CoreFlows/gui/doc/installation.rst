
.. include:: ./rst_prolog.rst

.. _iraInstallation:

**************************
Development installations
**************************

A development installation of solverlabGUI allows programmer improvments.
It is a classical usage of Python_ packages.
It needs a directory for Python_ interpreter (usually named *miniconda3*),
and another directory for solverlabGUI scripts (usually named *solverlabGUI*),

*In fine*, user could find (and use) command *solverlabGUI* directly after a
*detar/unzip* installation. Or a *git clone*.

.. warning:: #. **Windows7-10** all-in-one installation is a *development installation* of solverlabGUI,
                which allows programmer improvments.
             #. Users find *miniconda3* and *solverlabGUI* directories in parent directory named
                *C:\\Users\\Public\\solverlab*.
             #. **Linux** all-in-one installation PyInstaller_ bundle is **NOT** like that.
             #. To get *development installation* **Linux** (freely located) of solverlabGUI,
                users have to follow the two next chapters.


.. _iraInstallation_pythonlinux:

Python installation Linux
============================

To install python3 (and its mandatory packages PyQt5 etc.) *locally*, we suggest to use miniconda_.
Note that *miniconda* is windows7-10 compliant.

.. note:: You may use this Python interpreter for another python scripting code than solverlabGUI.

For information:

#. https://conda.io/miniconda.html
#. https://conda.io/docs/index.html

Example of install *(Linux-bash)*:

.. code-block:: bash

    bash Miniconda3-latest-Linux-x86_64.sh
    # -> Miniconda3 will now be installed into this location:
    # -> /volatile/common/miniconda3 (for example. It is located as you want.)
    # -> Thank you for installing Miniconda3!

    export PATH=/volatile/common/miniconda3/bin:$PATH
    which conda
    # -> /volatile/common/miniconda3/bin/conda

    conda create --name py3qt5 python=3 \
       pip jupyter matplotlib numpy pandas pandas-datareader \
       pyqt=5 scipy sympy jsonschema pyyaml libxml2 paramiko
    # -> Solving environment: done
    # -> Proceed ([y]/n)? y
    # -> Downloading and Extracting Packages
    # -> To activate this environment, use:
    # -> source activate py3qt5

    conda info --envs
    # -> conda environments:
    # ->   base         /volatile/common/miniconda3
    # ->   py3qt5       /volatile/common/miniconda3/envs/py3qt5
    # ->   etc...

    source activate py3qt5
    which python
    # - > /volatile/common/miniconda3/envs/py3qt5/bin/python


.. _iraInstallation_linux:

solverlabGUI installation Linux
=================================


.. warning:: Python interpreter py3qt5 is supposed to be set
             and useful in environment path.
             Usually command *source activate py3qt5* assume that.

Example of install/launch *(Linux-bash)*:

.. code-block:: bash

    cd whereYouWant
    tar -xf .../solverlabGUI_xxxx.tgz
    cd solverlabGUI
    ls -l solverlabGUI        # the launch executable command (is a script python)
    which python            # --> py3qt5
    ./solverlabGUI -h         # on line help
    ./solverlabGUI -g -w ...  # launch GUI


.. _iraInstallation_windows:

solverlabGUI installation Windows7-10
======================================

.. warning:: Python interpreter py3qt5 is supposed to be set
             and useful in environment path, at **mandatory** usual location
             *C:\\Users\\Public\\solverlab\\miniconda3*.
             Usually command *conda activate py3qt5* assume that.


Example of install/launch *(Windows7/10-cmd.exe shell)*, using 7-zip_:

.. code-block:: bat

    C:\
    cd C:\Users\Public\solverlab   # this is mandatory location, useful for all users
    "C:\Program Files\7-Zip\7z.exe" x .../solverlabGUI_xxxx.7z
    cd C:\Users\Public\solverlab\solverlabGUI
    where python                 # --> py3qt5
    python solverlabGUI -h         # on line help
    python solverlabGUI -g -w ...  # launch GUI


.. note:: To launch solverlab GUI, you may use Windows shortcut
          *C:\\User\\Public\\solverlab\\solverlabGUI\\LaunchSolverlabGUI(.lnk)*
          .




