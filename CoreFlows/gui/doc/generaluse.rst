
.. include:: ./rst_prolog.rst


.. _generalUse:

Usage of the GUI
-----------------

Create a new *Solverlab tree* by clicking on the first item in the **TOOLBAR**

Right click on the *fileMed* Field to import your .med you want to study.
Then Right click on *Model* to choose the physical model you want to apply on your file.

.. note:: All model are not implemented yet in the GUI.

You can now modify all the data on your model.

Before launching the simulation check *Analysis* to select how you want to launch the simulation.

.. warning:: file_name under *Model/computation_parameters* is the name of the result file.
             name under *Analysis/dataInformations* is the directory name where the result file will be stored. 
             
All model use *computation_parameters*.
(TODO explicit paramters with pictures)
 
To launch the simulation you have to go in the *Analysis* section of the **TREE VIEW**.
(picture of analysis field)
You can select where to save your file.
Choose to launch in foreground (will stop the gui until the solverlab is done) or background and launch another simulation.
