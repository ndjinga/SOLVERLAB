## TODO's

- **POLICY:** Les plus vieux en bas.
- **POLICY:** Les TODO's résolus, ajouter 'fixed 220610' pour résolu le 10/6/2022 (par exemple).


### 2022/06

- respecter la chronologie dans ce fichier, en attendant d'utiliser le bug tracker tuleap ou gitlab (par exemple).


### 2021/12

- extraire la pression autour de l'aile
- renomer la librairie medcoupling en medcouplingcpp
- créer une variable SOLVERLAB_WITH_KERNEL/SOLVERLAB_WITH_SALOME, charger le kernel salomé à un seul endroit
(actuellement cela est fait dans examples/python/CMakeLists.txt et gui/salome/CMakeLists.txt).
Echec fin Avril 2021.
Seule la gui a été chargée à un seul endroit.
Le kernel et la configuration sont chargés à deux endroits.
- Débugger Euler 2D
- En 3D MPI, résoudre le problème d'allocation de la matrice PETSC
Débuger l'erreur mémoire du test PoissionFE_3DTorus
- gestion des fichiers .ipynb depuis le makefile ou depuis chaque test (import nbconvert)
- comprendre pourquoi les tests wavesystem avec source ne convergent pas
- Finaliser l'élasticité linéaire
- Finaliser loi d'état
- Finaliser IJKmesh
- Débugger GUI->yacine
- modéliser L(p) ou L(T)
- TranportEquation :  créer un exemple miltiD lisant un maillage med
- Dans CoreFlows : ne plus stocker en mémoire le maillage après la première itération,
mais seulement les valeurs des champs dans des dataarray
