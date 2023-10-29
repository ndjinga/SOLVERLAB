
![logo](../images/logoSOLVERLAB.png)

## SOLVERLAB_GUI

### Introduction

- Is a new Graphical User Interface (GUI/IHM) developped for SOLVERLAB.  

- This GUI uses the library [PACKAGESPY](https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git)
  contained in the platform SALOME (variable `PACKAGESPY_ROOT_DIR`).  

- To compile the GUI, use the cmake option `-DSOLVERLAB_WITH_PACKAGESPY=ON`, 
  and provide the location of PACKAGESPY with the cmake option
  `-DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}`.  

- Finally, to launch the Graphical User Interface of SOLVERLAB,
  you have to load the SOLVERLAB environment.  
  *In your terminal*, type:

  ```bash
  source .../SOLVERLAB_install/env_SOLVERLAB.sh
  ${SOLVERLAB_ROOT_DIR}/solverlabGUI --help # help CLI
  ${SOLVERLAB_ROOT_DIR}/solverlabGUI --gui  # lanch GUI
  ```


### PACKAGESPY

- PACKAGESPY is the *main* mandatory prerequisite for SOLVERLAB_GUI.

- The *reference* base git for PACKAGESPY is
  `https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git`
  (requires an account to log in)  

- PACKAGESPY files can also be found in SALOME binary files (location given by the variable `PACKAGESPY_ROOT_DIR`).


