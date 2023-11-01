
![logo](../images/logoSOLVERLAB.png)

## SOLVERLAB_GUI

### Introduction

- Is a new Graphical User Interface (GUI/IHM) developped for SOLVERLAB.  

- This GUI uses the library [PACKAGESPY](https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git)
  contained in the platform SALOME (variable `PACKAGESPY_ROOT_DIR`).  

- To compile the GUI, provide the location of the following libraries :  
  - PACKAGESPY (cmake option  `-DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}`)  
  - SALOME configuration (cmake option `-DCONFIGURATION_ROOT_DIR=${-DCONFIGURATION_ROOT_DIR}$`)  
  - SALOME KERNEL (cmake option `-DKERNEL_ROOT_DIR=${-DKERNEL_ROOT_DIR}$`).  

- Finally, to launch the Graphical User Interface of SOLVERLAB,
  you have to load the SOLVERLAB environment.  
  *In your terminal*, type:

  ```bash
  source .../path/to/SALOME/env_launch.sh
  python3 $SOLVERLABGUI -g
  ```


### PACKAGESPY

- PACKAGESPY is the *main* mandatory prerequisite for SOLVERLAB_GUI.

- The *reference* base git for PACKAGESPY is
  `https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git`
  (requires an account to log in)  

- PACKAGESPY files can also be found in SALOME binary files (location given by the variable `PACKAGESPY_ROOT_DIR`).


