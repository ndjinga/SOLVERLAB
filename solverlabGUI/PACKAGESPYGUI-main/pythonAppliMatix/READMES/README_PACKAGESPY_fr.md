
## PACKAGESPY introduction

**TODO** améliorer ce fichier do in english

PACKAGESPY

- Est une bibliothèque python développée au [CEA](https://www.cea.fr)
  permettant la création d'interfaces graphiques.

- Est conçu pour une disposition graphique spécifique et utilise
  [PyQt5](https://pypi.org/project/PyQt5). Ici simplification et généricité commandent.

- Est basé, avec quelques libertés, sur le *design pattern*
  [Modele-Vue-Controleur (Model-View-Controller)](https://en.wikipedia.org/wiki/Model%E2%80%93view%E2%80%93controller).
  Le *Modele* est normalement indépendant de la *Vue* ils ne peuvent pas agir directement l'un sur l'autre.

- Le *Modele* et la *Vue* communiquent à l'aide de requêtes au coeur de
  l'[API](https://en.wikipedia.org/wiki/API) par l'intermédiaire du *Controleur*


### Modele

- Le *Modele* de PACKAGESPY est un arbre de données.
- Cette structure de données possede des noeuds/branches qui
  portent d'autres branches et/ou des feuilles.
- Seules les feuilles contiennent les valeurs.
- Il va contenir ce qu'on appelle les *données d'entrées*, modifiables par l'utilisateur.


- PACKAGESPY assure génériquement la sérialisation/désérialisation du Modele
( format de fichier [XML](https://fr.wikipedia.org/wiki/Extensible_Markup_Language) ).
Chaque noeud/branche a un nom et ses branches/feuilles *filles* sont définies
dans la déclaration de la classe python associée.


- Toutes les classes utilisées en tant que noeud/branche héritent de `BaseFreeXyz`.
PACKAGESPY fournit des classes abstraites (reconnaissables car commençant par `_`
ainsi que des classes basiques permettant de répondre a une grande parties des besoins
des structures de données nécessaires aux lancements des applications.


Une classe typique noeud/branche est définie a minima par 3 attributs & méthodes:
- `_attributesList` : attribut liste contenant les fils du noeud actuel.
- `_helpDict` : attribut dictionnaire contenant des *tooltips* à présenter à l'utilisateur
  pour chaque fils du noeud.
- `__init__()` : méthode définissant le comportement du noeud à l'usage.


```python
import xyzpy.classFactoryXyz as CLFX
from xyzpy.intFloatListXyz import StrInListXyz #only need to import class we want to derivate.
from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz

class AnimalList(StrInListXyz):
  _allowedList = ["Cat", "Dog", "Other"]

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

# make known new classes in factory
CLFX.appendAllXyzClasses([AnimalList, NodeExample, ListExample, MyModel])
```

`CLFX.appendAllXyzClasses(...)` est un *design pattern*
[*fabrique (factory)*](https://fr.wikipedia.org/wiki/Fabrique_(patron_de_conception)) qui permet d'informer
n'importe quelle partie du code de la présence des classes ajoutées en parametres.
Cela permet au code d'intancier une classe uniquement en connaissant son nom.


### Vue

PACKAGESPY a comme pérequis
[PyQt5](https://pypi.org/project/PyQt5), et utilise:

- Le widget *QMainWindows* ou son equivalent SALOME-Desktop.
- Un arbre *QTreeView* dans le dock de gauche,
- Une barre d'action *QToolBar* dans le dock du haut.
- Une fenêtre centrale pour afficher des widgets utiles complémentaires
  (un explorateur de fichier, une fenêtre plot matplotlib etc.).
- La classe *TreeXmlXyz* peut être directement instanciée et utilisera donc des réglages par défaut.
  On peut aussi la dériver et créer un affichage spécifique pour chaque application
  (couleurs et titres etc.).

```python
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
```


### Controleur

Le *Controleur* est la partie du code qui va gérer les interactions
entre le *Modele* qu'il protège (en mémoire)
et les actions de l'utilisateur sur le
[GUI](https://en.wikipedia.org/wiki/Graphical_user_interface).

Il gère les signaux [PyQt5](https://pypi.org/project/PyQt5)
qui sont lancés lorsque l'utilisateur clique sur les items du *TreeView*.

```python
```
