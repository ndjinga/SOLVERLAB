.. _doepyPPY:

Présentation de doepy
==================================

| doepy est un ensemble de programmes écrit en Python.
| C'est un utilitaire générique pour gérer des `plans d'expériences <https://fr.wikipedia.org/wiki/Plan_d%27exp%C3%A9riences>`_ (Design Of Experience: DOE_).
| Un plan d'expérience est une suite ordonnée d'essais d'une expérimentation. Un essai correspond ici au lancement unitaire d'un code de modélisation quelconque.
| doepy est prévu pour permettre de lancer (de nombreuses fois) de manière générique *tous* les codes de la plateforme MATIX_, et en première instance le code NUMODIS_.

.. _prerequis:


Prérequis
----------

| doepy est integré à la plateforme MATIX_. Il utilise certains prérequis de la plateforme.
| En particulier URANIE_, (prérequis principal de doepy).


Utilisation
------------

Ce package ne dispose pas de GUI il est utilisé en mode *ligne de commande*, dans une console bash.

Il faut préalablement activer l'environnement MATIX_.

.. code-block:: bash
  
  bash
  .../matix context
 
 
Premier exemple d'utilisation de URANIE_
--------------------------------------------

| Un exemple DOE *simple* est present dans le package doepy pour démonstration.
| Il consiste a lancer un pseudo-code python *runMyCode.py*.
| *runMyCode.py* lit dans le fichier *inMyCode.dat* un paramètre du DOE nommé *LongueurDislocation*. Il le multiplie par 2. Il range le résultat dans le fichier *outMyCode.dat* dans un paramètre de sortie du DOE nommé *resultat*.
| Il y a un tracé (sur l'écran) des résultats du DOE.
| Les résultats des essais du DOE sont rapatriés/concaténés dans le fichier *Uranie_MyCode1_exportData.dat* que l'on va retrouver dans le répertoire courant.

son lancement peut se résumer à cela (sous bash):

.. code-block:: bash
  
  bash
  .../matix context
  cd $PACKAGESPY_ROOT_DIR/pythonAppliMatix/doepy/exampleDoeUranie
  root -l Uranie_MyCode1.C #avec visu plot des résultats
  #n.b. tapez '.q' pour quitter interpréteur root URANIE
  gedit Uranie_MyCode1_exportData.dat #visu array des résultats
  


Création-lancement d'un DOE_ NUMODIS_
--------------------------------------------

Etapes pour créer son propre DOE: 

#. Copier *en local* doepy/doeNumodis/exampleForeman (ou un autre exemple qui marche).
#. Lancer le DOE *en local*, vérifier que tout marche toujours.
#. Modifier/renommer les fichiers pour adapter le DOE à *son propre* problème.
#. Lancer le DOE (a travers URANIE_).
#. Analyser les résultats, tracer des courbes, en particulier avec les utilitaires URANIE_.


Ce qui peut se résumer à cela (sous bash):

.. code-block:: bash
  
  bash
  .../matix context
  export MYDIR=/tmp/doe_myProblem #choose your name as your problem
  cd /tmp && rm -rf $MYDIR        #clean previous
  cp -r $PACKAGESPY_ROOT_DIR/pythonAppliMatix/doepy/doeNumodis/exampleForeman $MYDIR
  cd $MYDIR
  
  #test all right
  rm -rf /tmp/Doe_Foreman_uranie  #clean previous
  root -l -q Uranie_MyCode_Foreman.C
  
  #modify/rename and adapt... your problem
  gedit doeForeman.py Uranie_MyCode_Foreman.C sample/*.xml
  
  #launch your problem DOE
  root -l -q Uranie_MyCode_MyProblem.C #launch your problem DOE
  
  #supposed copied/modified/renamed as doeMyProblem.py Uranie_MyCode_MyProblem.C by user
  #clean initial files
  rm doeForeman.py Uranie_MyCode_Foreman.C
  

Contenu initial de exampleForeman (pour information):

.. code-block:: bash

  cd $PACKAGESPY_ROOT_DIR/pythonAppliMatix/doepy/doeNumodis/exampleForeman
  tree
  .
  ├── doeForeman.py
  ├── inMyCode.dat
  ├── rootlogon.C
  ├── sample
  │   ├── data_files.xml
  │   ├── donnees.xml
  │   ├── Fe.xml
  │   ├── graph.xml
  │   ├── micro.xml
  │   └── topo.xml
  └── Uranie_MyCode_Foreman.C


Options de mise en oeuvre du DOE NUMODIS
--------------------------------------------

Il est possible de lancer un DOE sans disposer de URANIE_, on est alors limité à un *simple* DOE_ de type Plan d'expérience complet. 

Nous conseillons fortement à l'utilisateur de lire la documentation URANIE_.

remarque:

Un plan d'expérience complet génère au minimum 2 puissance k essais pour k paramètres, et plus souvent n1*n2*n3*...*nk essais pour k paramètres, ou ni est le nombre de discrétisations du parametre i (souvent >> 2). Ce mode devient vite couteux en ressouces informatiques lorsque k (ou ni) est grand.

Exemple de Plan d'expérience complet sans disposer de URANIE_ (sous bash):

.. code-block:: bash
  
  bash
  .../matix context
  export MYDIR=/tmp/doe_myProblem #choose your name as your problem
  cd /tmp
  rm -rf $MYDIR
  cp -r $PACKAGESPY_ROOT_DIR/pythonAppliMatix/doepy/doeNumodis/exampleFrankRead $MYDIR
  cd $MYDIR
  $PACKAGESPY_ROOT_DIR/pythonAppliMatix/doepy/doeNumodis/runDoeNumodis.py \
    --doe ./doeFrankRead.py \
    --benchmark ./sample/data_files.xml \
    --workdir /tmp/Doe_FrankRead \
    --verbose
  


