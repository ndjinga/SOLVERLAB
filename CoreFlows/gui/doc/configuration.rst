
.. include:: ./rst_prolog.rst

.. _iraConfiguration:

****************************
SolverlabGUI configuration
****************************

SolverlabGUI uses files to store its configuration parameters.
It uses ConfigParser_ package, from The Python Standard Library.

Two configuration files are created or used at solverlabGUI launch,
and located at *SOLVERLABGUI_WORKDIR* directory.

#. file *.../SOLVERLABGUI_WORKDIR/solverlabGUI_user.cfg*
#. file *.../SOLVERLABGUI_WORKDIR/solverlabGUI_default.cfg*

.. note:: **User may edit/modify** file *solverlabGUI_user.cfg*


Syntax
========

See https://docs.python.org/3/library/configparser.html


.. _iraConfiguration_description:

Description
=============

The effective configuration **is a merge** of these two previous files.
Parameters in *solverlabGUI_user.cfg* override parameters in *solverlabGUI_default.cfg*.

User will find **all allowed parameters** in systematically up-to-dated
*solverlabGUI_default.cfg*.









