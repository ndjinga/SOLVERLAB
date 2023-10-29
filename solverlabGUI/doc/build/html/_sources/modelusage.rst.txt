
.. include:: ../rst_prolog.rst


.. _modelUsage:

Model
---------

This section will explicit all options of the implemented model in the GUI

.. image:: images/emptyTreeview.png
    :align: center

The Solverlab Tree Object you can see in the **TREE VIEW** broke down as follows:

* GeometryMed

    * fileMed: the .med file you want to work on (you can import one by right-clicking and select "Browse file")
    * contents: read the selected file and display some important information about its content

* Model

    * List of Model (Right-click to add one)

* Analysis

    * datainformations
        * name: name of the directory where it will save your work
        * directory: location of the directory
        
    * caseSolverlab
        * launchInBackground: 
        * Equation: Choose which simulation to launch from the List of Model
        * NumberOfProcessors: If in background try to launch solverlab in multicore

Diffusion Equation
+++++++++++++++++++

You can see the documentation `here <SolverlabDiffusion_>`_

.. image:: images/diffusionequationtree.png
    :align: center

A field need to be present in your mesh file to be visible in the GUI. All "field_option" are here for advanced user, it is recommanded to leave them with default value.

Some value can be a field present in the mesh file or a scalar and the GUI let you choose between those two options. 

The boundary condition are created dynamically by reading in the .med file.
