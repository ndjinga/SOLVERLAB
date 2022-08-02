
![logo](../images/logoSOLVERLAB.png)

## SOLVERLAB_GUI

### Introduction

- Is a new Graphical User Interface (GUI-IHM) developped for SOLVERLAB,
  this GUI uses the library PACKAGESPY.

- To use the GUI, use the cmake option `-DSOLVERLAB_WITH_PACKAGESPY=ON`.
  That will download the library [PACKAGESPYGUI](https://github.com/ndjinga/PACKAGESPYGUI) from github.

- **OR** if you have a local copy of the library
  [PACKAGESPY](https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git)
  (as in SALOME tuleap contex), give the link to solverlab using the cmake option
  `-DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}`.

- Finally, to launch the Graphical User Interface of SOLVERLAB,
  you have to load the SOLVERLAB environment *in your terminal*, type:

  ```bash
  source .../SOLVERLAB_install/env_SOLVERLAB.sh
  $SOLVERLABGUI -h  # help CLI
  $SOLVERLABGUI -g  # lanch GUI
  ```


### PACKAGESPY (Tuleap) *outside* SOLVERLAB_GUI prerequisite

- PACKAGESPY is the *main* mandatory prerequisite for SOLVERLAB_GUI.

- The *reference* (SALOME-MATIX) base git for PACKAGESPY is
  `https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git`

- NOTE: This web site is not reacheable for everybody, unhappily.



### PACKAGESPYGUI (Github) *inside* SOLVERLAB_GUI prerequisite

- To clone-get *one version* of PACKAGESPY in SOLVERLAB_GUI, see
  [PACKAGESPYGUI](https://github.com/ndjinga/PACKAGESPYGUI) from github.

- WARNING: it is a *lightened* version of reference base git tuleap-SALOME-MATIX packagespy

- If this script launched with success, *authorized user*
  may set the new PACKAGESPY-SOLVERLAB_GUI
  prerequisite reference in base git SOLVERLAB.

- Commiting his new tag (...or his new branch) in
  https://github.com/ndjinga/SOLVERLAB.git (for example),
  *authorized user* allows *every users* to get and use SOLVERLAB_GUI.

- WARNING: This is a **FORK!** (alas unhappily!) (TODO: *avoid fork!, ever!*)
