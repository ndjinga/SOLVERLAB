
.. include:: ./rst_prolog.rst


.. _iraPrerequisites:

***************
Prerequisites
***************

There are some definitions, and links.

#. Solverlab_ code, and its `manual <SolverlabManualPdf_>`_.
#. PYTHON_ 3.5, with packages PyQt5_, numpy, matplotlib_, pandas_, etc. (usually named *py3qt5*).
#. PyInstaller_ 3.4, free sofware to make solverlabGUI bundle (*only* valid for Linux).
#. 7-zip_, free sofware to compress/uncompress .7z files (for Windows installations).

Installation needs a PYTHON_ 3.6 interpreter, which is included in
*All-in-one* solverlab installations (Linux *and* Windows).
see :ref:`iraInstallation_pythonlinux`.

PyInstaller_ is a program that freezes Python programs
in *bundle*, which is **almost** a python package.
For more information about *bundle*, see `PyInstaller manual <PyInstallerManual_>`_.


.. _iraAllInOneInstallation:

*********************************
All-in-one installation
*********************************

These installations contain in **one** compressed file:

#. the solverlabGUI python scripts.
#. An interpreter PYTHON_.
#. The Solverlab code (GPL), its **source files**, and
   two executable files, one for Linux and one Windows.
#. The useful Corteo_ data base (4bits).


.. warning:: Corteo data base used version in solverlabGUI is **NOT** Version 20160816.


*********************************
All-in-one installation Linux
*********************************

.. warning:: This is a PyInstaller_ bundle,
             installed locally **where users want**.
             Python interpreter (named py3qt5) is simultaneously installed,


Source tar file *solverlabGUI_bundle_xxxx.tgz* is the **one** compressed archive file
of the PyInstaller bundle in *one folder* mode.
See `more information here <https://pyinstaller.readthedocs.io/en/stable/operating-mode.html#bundling-to-one-folder>`_.

There are two ways to install solverlabGUI:

#. Install solverlabGUI directly, using usual *file manager* functionalities:
   uncompress tar file in user's choice directory.
#. Install solverlabGUI typing bash command, in usual *terminal*:

    .. code-block:: bash

        # install
        cd yourChoiceDirectory       # which is really where you want
        tar -xf .../solverlabGUI_bundle_xxxx.tgz
        # launch
        cd solverlabGUI_bundle         # folder name as linux pyinstaller bundle
        ./solverlabGUI -h              # on line help
        ./solverlabGUI -g -w ...       # launch GUI


*********************************************
All-in-one installation Windows7-10
*********************************************

.. warning:: #. This is **not** a PyInstaller_ bundle.
             #. The **mandatory located** root directory is *C:\\Users\\Public\\solverlab*.
             #. The **mandatory located** solverlabGUI directory
                is *C:\\Users\\Public\\solverlab\\solverlabGUI*.
             #. The **mandatory located** Python interpreter py3qt5 directory
                is *C:\\Users\\Public\\solverlab\\miniconda3*.
             #. The sofware tool to uncompress .7z files is 7-zip_, which **has to be installed**.

Source .7z file *solverlabGUI_xxxx.7z* is the **one** compressed archive file.

There are two ways to install solverlabGUI:

#. Install solverlabGUI directly, using usual *file manager* functionalities:
   uncompress .7z file in mandatory *C:\\Users\\Public* directory.
#. Install solverlabGUI typing DOS command, in usual *(Windows7/10-cmd.exe shell)*:

    .. code-block:: bat

        rem install
        C:\
        cd C:\Users\Public   # this is mandatory location, useful for all users
        "C:\Program Files\7-Zip\7z.exe" x .../solverlab_xxxx.7z
        rem launch
        C:\User\Public\solverlab\solverlabGUI\LaunchSolverlabGUI.bat

.. note:: To launch solverlab GUI, you may use Windows shortcut
          *C:\\User\\Public\\solverlab\\solverlabGUI\\LaunchSolverlabGUI(.lnk)*

          .
