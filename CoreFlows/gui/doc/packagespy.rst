
.. include:: ../rst_prolog.rst

===============================================
How to add a new physical model to SolverlabGUI
===============================================

PACKAGESPY_ est une bibliothèque python développée au CEA_ permettant la création d'interfaces graphiques. Elle est basé sur le design du *Model-View-Controller* (voir ModelViewController_).
Le *Model* est une structure de donnée de type Arbre où les feuilles contiennent les valeurs.
(Chaque partie est normalement indépendante et ne peut pas agir directement sur l'autre. Ainsi, elle communique à l'aide de request au coeur de l'api_).
    
L'ajout d'un nouveau modèle physique ie d'une nouvelle équation mathématique utilisable dans l'interface devrait normalement se limiter à la creation d'un *Model* spécifique pour l'équation choisie et du script faisant le lien entre celui ci et SOLVERLAB_.

EquationSvl est une classe qui ne doit pas etre instanciée, elle contient tous les paramètres commun aux modèles ainsi que ceux necessaires au fonctionnement de SOLVERLAB_.
Il faut donc créer une nouvelle classe pour accueillir les données specifiques au modèle physique que l'on veut implémenter, en initialisant les données avec les valeurs par defaut de SOLVERLAB_.

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
          
          
Il faut ensuite créer un script qui va faire la liaison entre les données du *Model* et SOLVERLAB_.
C'est un script qui va parcourir le *Model* tout en transmettant les données à l'api_ de SOLVERLAB_.
Le script déjà implémenté se nomme launchDiffusionEQ.py.

Enfin il faut modifier la classe ListOfEquation dans equationsvl.py afin d'y autoriser la nouvelle classe,
ainsi que permettre à RunSolverlab.py de lancer le script de calcul de la nouvelle classe.

































