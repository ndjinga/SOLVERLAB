
.. include:: ../rst_prolog.rst


.. _packagespy:


=============
Packagespy
=============

Packagespy est une bibliothèque python développée au CEA permettant la création d'interfaces graphiques. Elle est basé sur le design du *Model-View-Controller*.
Le *Model* est une structure de donnée de type Arbre où les feuilles contiennent les valeurs.
(Chaque partie est normalement indépendante et ne peut pas agir directement sur l'autre. Ainsi, elle communique à l'aide de request au coeur de l'api)

*View*
--------------------

Packagespy a d'abord été pensé pour une disposition graphique spécifique et est basé sur PyQt5. Un arbre (TREEVIEW) dans le dock de gauche, un barre d'action (TOOLBAR) dans le dock du haut et une fenêtre centrale pour afficher du contenu.
La classe TreeXmlXyz peut être directement instanciée et utilisera donc des réglages par defaut.
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



*Model*
-----------

Le *Model* est la structure de donnée qui va enregistrer les données entrées par l'utilisateur. Packagespy à été concu pour fonctionner avec un arbre de donnée. Chaque noeud comportant une *Key*, son nom et une *Value*, sa classe.

Toute les classes utilisées en tant que noeud doivent hériter de xyz.BaseFreeXyz
Packagespy fournit des classes abstraites (reconnaissables car commançant par "_" ) ainsi que des classes basiques permettant de répondre a une grande parties des besoins de la plupart des structures de données necessaires aux lancements des applications.

Une classe typique est composé de 3 éléments:

- _attributesList : une liste contenant les fils du noeud actuel.
- _helpDict : un dictionnaire contenant des tooltips à afficher à l'utilisateur pour chaque fils du noeud.
- __init__ : la fonction qui va definir le comportement à l'initialisation.

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
      
CLFX.appendAllXyzClasses() est une méthode qui permet d'informer n'importe quelle partie du code de la présence des classes ajoutées en parametres. Ca permet au code d'intancier un classe uniquement en connaissant son nom.  

*Controller*
------------

Le *Controller* est la partie du code qui va gérer les interactions entre le *Model* en mémoire et les actions de l'utilisateur sur la fenêtre ainsi que celles avec le code sur lequel la GUI s'appuie.
Il faut donc créer des signaux PyQt qui vont ensuite pouvoir être interceptés lorsque l'utilisateur va faire des actions sur la GUI pour pouvoir répondre en conséquence.

.. code-block:: python

    
    
    
    
Ajouter un modèle physique dans SolverlabGui
--------------------------------------------

L'ajout d'un nouveau modèle physique ie d'une nouvelle équation mathématique utilisable dans l'interface devrait normalement se limiter à la creation d'un *Model* spécifique pour l'équation choisie et du script faisant le lien entre celui ci et solverlab.

EquationSvl est une classe qui ne doit pas etre instanciée, elle contient tous les paramètres commun aux modèles ainsi que ceux necessaires au fonctionnement de solverlab.
Il faut donc créer une nouvelle classe pour accueillir les données specifiques au modèle physique que l'on veut implémenter, en initialisant les données avec les valeurs par defaut de solverlab.

Pour exemple:

.. code-block:: python

    class TransportEq(_XyzConstraintBase):
        _attributesList = [
          ("field_name", "SelectField"),
          ("field_option", "FieldOptionSvl"),
          ("SpecificValue","NewValueClass"),
          ("SpecificBoundaryCondition","NewBoundaryClass"),
          ("numerical_parameters", "NumericalSvl"),
          ("computation_parameters", "ComputationSvl"),
          ]
          
          
Il faut ensuite créer un script qui va faire la liaison entre les données du *Model* et solverlab.
C'est un script qui va parcourir le *Model* tout en transmettant les données à l'api solverlab.
Le script déjà implémenté se nomme launchDiffusionEQ.py.

Enfin il faut modifier la classe ListOfEquation dans equationsvl.py afin d'y autoriser la nouvelle classe,
ainsi que permettre à RunSolverlab.py de lancer le script de calcul de la nouvelle classe.

































