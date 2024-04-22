
## QUICK TESTS


Easy tests *a mano*.


### quick test python-SALOME

With a SALOME or MATIX configuration/installation.

```bash
matix context   # or ...
salome context

python

import numpy as np
np.__file__

# CoreFlows solverlab
import CoreFlows as cf
cf.__file__

# cdmath solverlab
import cdmath as cm
cm.__file__

# GUI solverlab
import solverlabpy as SVL
SVL.__file__

# PACKAGESPY pythonAppliMatix
import xyzpy as XYZ
XYZ.__file__
help(XYZ)
```


### tricks


To explore and debug, if import problems.

```bash
export ADIR=/volatile2/wambeke/TULEAP_MATIX3/MATIX_30-CO7

${ADIR}/matix context
which solverlabGUI
solverlabGUI -g

# or...
# export PYTHONPATH=/volatile2/wambeke/TULEAP_MATIX3/MATIX_30-CO7/INSTALL/SOLVERLAB/bin:${PYTHONPATH} # TODO miss 2206 fixed now

export SLVGUI=$(find ${ADIR}/SOURCES/SOLVERLAB -name "solverlabGUI")
$SLVGUI -g

tree ${ADIR}/SOURCES/PACKAGESPY/pythonAppliMatix
find ${ADIR}/SOURCES/PACKAGESPY -name "helppy"

cd ${ADIR}
${ADIR}/matix
```
