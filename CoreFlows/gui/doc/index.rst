
.. include:: ./rst_prolog.rst

.. empty first toctree is used for pdf contents maxdepth, see sphinx/builders/latex/__init__.py, toctrees[0].get('maxdepth')

.. toctree::
   :maxdepth: 2


.. _iraHome:

**************
Solverlab GUI
**************

.. image:: images/icon_about_SolverlabGUI.png
   :align: center

.. note:: Find *Solverlab code* documentation `here <SolverlabManualPdf_>`_ [1]_.

The Solverlab GUI code is a GUI_ (Graphical User Interface)
used to perform simulations with `Solverlab physical models <SolverlabPresentationPdf_>`_.

The Solverlab GUI is designed with (Packagespy based on) PyQt5_ and Python3_ scripts files.



User's manual
==============

.. toctree::
   :maxdepth: 2

   svlinstallation
   iramainwidget
   generaluse
   modelusage



Programmer's guide
====================

.. toctree::
   :maxdepth: 1

   background

Packagespy 
============   
(TODO packagespy doc)

.. toctree::
   :maxdepth: 1
   
   packagespy


Release Notes
=============

.. toctree::
   :maxdepth: 1

   Release Notes 9.8.0 <release_notes/release_notes_9.8.0>

.. [1] solverlabGUI/doc/src/solverlabDocuments/solverlab_manual.pdf
