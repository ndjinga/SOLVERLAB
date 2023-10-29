.. _essaisphinx:

exemples Sphinx
==================================

allez voir `doc sphinx <http://openmdao.org/dev_docs/documenting/sphinx.html>`_

Du code python:

.. code-block:: python
   :linenos:
   :emphasize-lines: 3,5

   """:emphasize-lines:: marche pas sans Pygments!"""
   def some_function_code-block():
       interesting = False
       print 'This line is highlighted.'
       print 'This one is not...'
       print '...but this one is.'

.. code-block:: python

   """.. highlight:: python marche pas sans Pygments!"""
   def some_function_highlight():
       interesting = False
       print 'This line is highlighted.'
       print 'This one is not...'
       print '...but this one is.'

une liste.

* liste1
* liste2
* liste3

une liste numérotée.

#. listea
#. listeb
#. listeb

aller a hello_ dans cette page.

aller a `un texte lien namefile.html (home page) <index.html>`_ hors de la page.

aller a :ref:`un autre texte lien namePageRst (home page) <entrydocpackagespy>` hors de la page.

::

 plein de texte *tel que sans formattage mais dans une zone speciale*
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte
 plein de texte

.. _hello:

Hello c'est ici!!! dans la meme page

| plein de texte *tel que* mais avec du formattage
| plein de texte
| plein de texte
| plein de texte
| plein de texte

| plein de texte
| plein de texte
| plein de texte
| plein de texte
| plein de texte
| plein de texte
