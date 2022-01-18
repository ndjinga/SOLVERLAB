
.. include:: ./rst_prolog.rst

.. _iraExploreDirWidget:

Explore Dir widget
--------------------------

This widget displays the contents of user solverlabGUI working directory.
This directory is usually referenced as *SOLVERLABGUI_WORKDIR*.

Its usage is like a **simple** file explorer.

.. note:: Theses result files of Solverlab code are located in sub-directories named *output*.


.. image:: images/solverlabExploreDirWidget1.png
   :scale: 80 %
   :align: center

There are three main widgets (from left to right):

#. Directories names widget
#. Files names widget
#. File contents widget

.. note:: The lower *selected files* widget is for future improvments, no usage *yet*.


There are some contextual menus to explore directories, and to display input/output text files.


Directories names widget menu
::::::::::::::::::::::::::::::

.. image:: images/solverlabExploreDirectoriesMenu1.png
   :scale: 80 %
   :align: center

This menu contains some elementary actions to navigate in **all** disk directories.


.. _iraFilesNamesMenu:

Files names widget menu
:::::::::::::::::::::::::

.. image:: images/solverlabExploreFilesMenu1.png
   :scale: 80 %
   :align: center


This menu contains some actions to apply on selected file.


Files contents widget menu
::::::::::::::::::::::::::::::::::

.. image:: images/solverlabExploreDisplayMenu1.png
   :scale: 80 %
   :align: center

This menu contains some elementary actions to apply on displayed file.
The files are *syntax highlighted* if possible, using highlight_ tool, only on Linux distributions for now.

.. warning:: The lower *highlight theme* action is displayed only for Linux distribution.


