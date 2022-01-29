
.. include:: ../rst_prolog.rst

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

The Solverlab GUI is designed with the library PACKAGESPY_, which is based on PyQt5_ and Python3_ scripts files.



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
   packagespy


Release Notes
=============

.. toctree::
   :maxdepth: 1

   Release Notes 1.0.0 <release_notes/release_notes_1.0.0>

.. [1] solverlabGUI/doc/src/solverlabDocuments/20140804_solverlab_manual.pdf
