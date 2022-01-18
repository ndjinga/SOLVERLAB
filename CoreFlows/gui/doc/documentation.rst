
.. include:: ./rst_prolog.rst

.. _iraDocumentation:

************************************
Documentation
************************************

.. _iraDocumentation_consultation:


Doc consultation
=================

To display solverlabGUI html documentation in your web browser *firefox*, or *else*.
The initial entry file is located at *solverlabGUI/doc/build/html/index.html*.

.. code-block:: bash

  # Linux bash, as an example
  cd .../solverlabGUI
  firefox doc/build/html/index.html &
  # or as CLI_
  solverlabGUI --doc
  
.. _iraDocumentation_modification:

Doc modification
=================

To modify solverlabGUI documentation with simple editor *pluma*, or else.

Read the manual, see http://www.sphinx-doc.org/en/stable/tutorial.html,
or may be copy/paste from 'Show Source' item.

.. code-block:: bash

  # Linux bash, as an example
  cd ...solverlabGUI/
  tmp=$(find doc -name "*.rst")
  pluma $tmp &


.. _iraDocumentation_compilation:

Doc compilation Linux
======================

On a Linux system, to compile solverlabGUI html documentation,
programmers use installed GNU_ *make*, and SPHINX_.

.. warning:: To make documentation **pdf** programmers needs
             installed *texlive* package (preferably up to date version).
             See: https://www.tug.org/texlive/quickinstall.html


.. code-block:: text

  cd ...solverlabGUI/doc
  cat README  # read some environment setup information
  #  ... and read it
  make
    Please use `make <target>' where <target> is one of
    html       to make standalone HTML files
    dirhtml    to make HTML files named index.html in directories
    singlehtml to make a single large HTML file
    pickle     to make pickle files
    json       to make JSON files
    htmlhelp   to make HTML files and a HTML help project
    qthelp     to make HTML files and a qthelp project
    devhelp    to make HTML files and a Devhelp project
    epub       to make an epub
    latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter
    latexpdf   to make LaTeX files and run them through pdflatex
    text       to make text files
    man        to make manual pages
    changes    to make an overview of all changed/added/deprecated items
    linkcheck  to check all external links for integrity
    doctest    to run all doctests embedded in the documentation (if enabled)
  
  # and then
  make html      # make html doc
  make latexpdf  # make pdf doc
  



