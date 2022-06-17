
.. include:: ../rst_prolog.rst


.. _generalUse:

=======
Usage
=======

Usage of the script
---------------------

SolverlabGUI usage is a Command Line Interface (CLI_), which is
Windows *and* Linux compatible (but only tested on Linux for the moment).

.. code-block:: bash

  solverlabGUI --[options]


Options of solverlabGUI
............................

Useful but *not exhaustive* generic options of *solverlabGUI* CLI.


Option *--help or -h*
............................

Get help as simple text.

.. code-block:: bash

    solverlabGUI --help          # get list of existing options


Option *--doc or -d*
............................

Get documentation as browser html.

.. code-block:: bash

    solverlabGUI --doc           # see html doc


Option *--verbose or -v*
............................

Change verbosity level (default is 'info').

.. code-block:: bash

    # execute solverlabGUI command in verbose debug mode
    solverlabGUI -v debug


Option *--workdir or -w*
............................

Change working directory (user data directory). Default is
../SOLVERLAB_WORKDIR

.. code-block:: bash

    # execute solverlab GUI in user choice working directory
    solverlabGUI -w .../MY_WORKDIR



Usage of the GUI
-----------------

Create a new *Solverlab tree* by clicking on the first item in the **TOOLBAR**

.. image:: images/guitips1.png
    :align: center

Right click on the *fileMed* Field to import your .med you want to study.

.. image:: images/guitips2.png
    :align: center
    
Then Right click on *Model* to choose the physical model you want to apply on your file.

.. image:: images/guitips3.png
    :align: center

.. note:: All model are not implemented yet in the GUI.

You can now modify all the data on your model.
             
All model use *computation_parameters* and *numerical_parameters*.

.. image:: images/guitips5.png
    :align: center
 
*file_name* is the name of the result file solverlab at the end of the simulation. 

To launch the simulation you have to go in the *Analysis* section of the **TREE VIEW**.
Before launching the simulation check *Analysis* to select how you want to launch the simulation.

.. image:: images/guitips4.png
    :align: center

*name* under *Analysis/dataInformations* is the directory name where the result file will be stored. 

Choose to launch in foreground (this can stop the gui for a long time if the simulation is big) or background and launch another simulation.
