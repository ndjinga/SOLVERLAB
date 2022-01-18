
.. include:: ./rst_prolog.rst


.. _packagespy:


=============
Packagespy
=============

Packagespy est une bibliothèque python développé au CEA permettant la création d'interface graphique. Elle est basé sur le design du Model-Vue-Controller.
Le Model est une structure de donnée de type Arbre ou les feuilles contiennent les valeurs.
(Chaque partie est normalement indépendante et ne peux pas agir directement sur l'autre ainsi elle communique la l'aide de request au coeur de l'api)

Vue
--------------------

Packagespy a d'abord été pensé pour une disposition graphique spécifique. Un arbre (TREEVIEW) dans le dock de gauche, un barre d'action (TOOLBAR) dans le dock du haut et une fenetre central pour afficher du contenu.


Model
-----------

Le Model est la structure de donnée qui va enregistrer les données entrées par l'utilisateur. Packagespy à été concu pour fonctionner avec un arbre de donnée. Chaque noeud comportant une Key, son nom et une Value, sa classe.

Toute les classes utilisées en tant que noeud doivent hérité de xyz.BaseFreeXyz
Packagespy fournit des classes abstraites (reconnaissable car commançant par "_" ) ainsi que des classes basique permettant de répondre a une grande parties des besoins de la plupart des structures de données necessaire aux lancements des applications.

Une classe typique est composé de 3 éléments:

- _attributesList : une liste contenant les fils du noeud actuel.
- _helpDict : un dictionnaire contenant des tooltips a afficher a l'utilisateur pour chaque fils du noeud.
- __init__ : la fonction qui va definir le comportement a l'initialisation.

.. code-block:: python

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
      
    class ListExample(ListOfBaseXyz)
      _allowedClasses = ["NodeExample"]

    class MyModel(_XyzConstrainBase):
      
      _atributesList = [
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
      
      
Controller
------------

Le Controller est la partie du code qui va gérer les interactions entre le Model en mémoire et les actions de l'utilisateur sur la fenetre ainsi que celle avec le code sur lequel la GUI s'appuie.

