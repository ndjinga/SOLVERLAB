.. _packagespyPPY:

Présentation de PACKAGESPY
==================================

PACKAGESPY est un ensemble de programmes écrit en Python.

* Il utilise entre autres la bibliothèque PyQt5.
* Il n'est *pas encore* un `package python <https://packaging.python.org>`_.

.. _prerequis:

Prérequis
----------
Il faut vérifier la présence des librairies python.

Les principaux prérequis conseillés validés optimaux sont ceux tirés de SALOME_,
en fonction de la version de SALOME

* Linux OS ()
* python 3
* PyQt 5
* matplotlib
* numpy
* pandas
* `sphinx <http://sphinx-doc.org>`_

Base git initiale
------------------------------------------

Pour récupérer tout PACKAGESPY *sur les ordinateurs CEA/DEN/DM2S/LGLS*.

.. code-block:: bash

  echo MYDIR="/home/MyChoice/tmp"
  mkdir $MYDIR
  cd $MYDIR
  ## recuperation package entier
  git clone /home/matix/GitRepo/packagespy.git

La documentation
------------------------------------------

| Elle est générée par `SPHINX <http://sphinx-doc.org>`_.
| C'est **préférentiellement** une doc html.
| Elle est dans le répertoire .../PACKAGESPY/doc/html.

Elle s'affiche avec votre browser préféré:

.. code-block:: bash

  matix context
  firefox $PACKAGESPY_ROOT_DIR/doc/html/index.html

Elle peut s'imprimer, presque bien (si la generation latexpdf est faite):

.. code-block:: bash

  evince $PACKAGESPY_ROOT_DIR/doc/packagespy.pdf

Elle est le résultat d'une compilation SPHINX sur les fichiers packagespy/doc/xxxx.rst:

.. code-block:: bash

  ## la documentation
  cd $PACKAGESPY_SRC_DIR/doc
  ls *.rst
  cd $PACKAGESPY_SRC_DIR/../../BUILD/PACKAGESPY/doc
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
    texinfo    to make Texinfo files
    info       to make Texinfo files and run them through makeinfo
    gettext    to make PO message catalogs
    changes    to make an overview of all changed/added/deprecated items
    linkcheck  to check all external links for integrity
    doctest    to run all doctests embedded in the documentation (if enabled)
  make html  # creation _build/html/index.html
  mv ./_build/latex/packagespy.aux ./_build/latex/packagespy.aux_en # solve french problem
  make latexpdf  # creation _build/latex/packagespy.pdf
  make singlehtml  # creation _build/singlehtml/index.html
