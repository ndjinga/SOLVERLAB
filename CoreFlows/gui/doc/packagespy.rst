
.. include:: ../rst_prolog.rst

===============================================
How to add a new physical model to SolverlabGUI
===============================================

.. _packagespy:


=============
Packagespy
=============

Packagespy est une bibliothèque python développée au CEA permettant la création simplifiée d'interfaces graphique (IHMs). 
Elle est basé sur le *design-pattern* Modèle-Vue-Controleur (MVC).

Packagespy a d'abord été conçu pour une disposition graphique spécifique basée sur PyQt V5.


Vue
--------------------

Un arbre (TREEVIEW) dans le dock de gauche, un barre d'action (TOOLBAR) dans le dock du haut et une fenetre centrale pour afficher du contenu.
La classe TreeXmlXyz peut être directement instancié et utilisera donc des réglages par defaut.
On peut aussi la dériver et créer un affichage spécifique pour une application. 

.. code-block:: python

    class TreeViewSvl(TreeXmlXyz):
      class COLS:
      labels = ['Name', 'Value', 'Attributes']
      Tag = 0
      Text = 1
      Attributes = 2

      def __init__(self, parent=None):
      super(TreeViewSvl, self).__init__(parent)

      self.setHeaderLabels(self.COLS.labels)
      self.setAlternatingRowColors(True)
      self.formats_treeview = FORMATS_TREEVIEW



Modèle
-----------

Le Modèle est une structure bien organisée de données de type Arbre où les feuilles contiennent 
les valeurs élémentaires (modifiables).

Le Modèle est la structure de donnée qui va enregistrer les données entrées par l'utilisateur. 
Packagespy à été concu pour fonctionner avec un arbre de donnée. 
Chaque noeud comportant une clé, son nom et une valeur, à chaque noeud, une classe.

Toutes les classes utilisées en tant que noeud héritent de quelques classes de base de packagespy.
Des classes abstraites (reconnaissables lorsque commençant par "_" ).
Ainsi que des classes génériques tres courantes
(flottant, entier, string, entier positif, etc...).
permettant de répondre immediatement à une grande partie des définitions élémentaires
d'une structure de données arborescente qui représente les données d'entrée d'une applications particulière.


Une classe typique est composé de 3 éléments:

- _attributesList : une liste contenant les fils du noeud actuel.
- _helpDict : un dictionnaire contenant des tooltips a afficher a l'utilisateur pour chaque fils du noeud.
- __init__ : la fonction qui va definir le comportement a l'initialisation.

.. code-block:: python

    import xyzpy.classFactoryXyz as CLFX
    from xyzpy.intFloatListXyz import StrInListXyz #only need to import class we want to derivate.
    from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz
    
    class AnimalList(StrInListXyz):
      _allowedList = ["None", "Cat", "Dog", "Other"]

    class NodeExample(_XyzContrainBase):
      _attributesList = [
      ("Name","StringXyz"),
      ("Room","IntPosXyz"),
      ("Animal","AnimalList")
      ]
      
      _helpDict = {
      "Name": ("Name of the customer",""),
      "Room": ("Room of the customer",""),
      "Animal":("Which animal is with him ?",""),
      }
      
       def __init__(self):
         super(NodeExample, self).__init__()
         self.setIsCast(True)
         self._setAllAttributesList()
      
    class ListExample(ListOfBaseXy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z)
      _allowedClasses = ["NodeExample"]

    class MyModel(_XyzConstrainBase):
      
      _attributesList = [
      ("Customers","ListExample"),
      ]
      
      _helpDict = {
      "Customers": ("List of all actual customer","")
      }
      
      def __init__(self):
        super(MyModel, self).__init__()
        self.setIsCast(True)
        self._defautNameAsRoot = "Hotel"
        self._setAllAttributesList()
      
    CLFX.appendAllXyzClasses([AnimalList, NodeExample, ListExample, MyModel]) 
      
CLFX.appendAllXyzClasses() est une méthode qui permet d'informer n'importe quel partie du code de la présence des classes ajoutées en parametres. Ca permet au code d'intancier un classe uniquement en connaissant son nom.  

Controller
------------

Le Controller est la partie du code qui va gérer les interactions entre le Model en mémoire et les actions de l'utilisateur sur la fenetre ainsi que celle avec le code sur lequel la GUI s'appuie.

.. code-block:: python

    
    
    
    
Ajouter un model dans SolverlabGui
------------------------------------

L'ajout d'un nouveau modele utilisable dans l'interface devrait normalement se limiter à la creation d'un model spécifique pour l'équation choisi et du script faisant le lien entre celui ci et solverlab.

EquationSvl est une classe qui ne doit pas etre instanciée mais qui contient tous les paramètres commun à chaque modèle ainsi que ceux necessaires au fonctionnement de solverlab.
Il faut donc créer une nouvelle classe pour acceuillir les données specifiques au modele que l'on veut implémenter.

































