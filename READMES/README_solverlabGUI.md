
![logo](../images/logoSOLVERLAB.png)

SOLVERLAB_GUI
==============

Introduction
------------

- Is a new Graphical User Interface (GUI-IHM) based on the library PACKAGESPY, 
which is being developped for SOLVERLAB. 

- In order to use the new GUI use the cmake option `-DSOLVERLAB_WITH_PACKAGESPY=ON`.
That will download the library [PACKAGESPYGUI](https://github.com/ndjinga/PACKAGESPYGUI) from github.

- If you have a local copy of the library PACKAGESPY, 
give the link to solverlab using the cmake option `-DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}`.

- To use the Graphical User Interface of SOLVERLAB, 
you can load the SOLVERLAB environment **in your terminal**, type:

    ```
    source .../SOLVERLAB_install/env_SOLVERLAB.sh`
    $SOLVERLABGUI -h  # to get CLI help
    $SOLVERLABGUI -g  # to get GUI
    ```

PACKAGESPY prerequisite
-------------------------

- PACKAGESPY is the *main* mandatory prerequisite for SOLVERLAB_**GUI**.

- The **reference** base git for PACKAGESPY is 
https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git

- NOTE: This web site is not reacheable for everybody, unhappily.


PACKAGESPY_*in*_SOLVERLAB_GUI prerequisite
-----------------------------------------

- To clone (get) **one version** of PACKAGESPY in SOLVERLAB_GUI,
as prerequisite **customized and integrated**
in base git SOLVERLAB, **advanced authorized user** may use 
bash script `...SOLVERLAB/scripts/get_packagespy.sh`.

- WARNING: it is a **lightened** version of reference base git matix-tuleap-packagespy

- If this script launched with success, **authorized user** 
may set the new PACKAGESPY-SOLVERLAB_**GUI** 
prerequisite reference in base git SOLVERLAB. 

- Commiting his new tag (...or his new branch) in 
https://github.com/ndjinga/SOLVERLAB.git (for example), 
**authorized user** allows *every users* to get and use SOLVERLAB_**GUI**.

- WARNING: This is a **FORK!** (alas unhappily!) (TODO: *avoid fork!, ever!*)




