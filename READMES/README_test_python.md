
QUICK TESTS
============


Eazy tests *a mano*.


quick test python 1
--------------------

With a SALOME or MATIX configuration/installation


```
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


tricks
=======

To explore and debug if import problems

```
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

/*** http://patorjk.com/software/taag/#p=display&v=1&c=c&f=Big&t=SolverlabGUI%0A%0A%0A%0A
 *       _____       _                _       _      _____ _    _ _____ 
 *      / ____|     | |              | |     | |    / ____| |  | |_   _|
 *     | (___   ___ | |_   _____ _ __| | __ _| |__ | |  __| |  | | | |  
 *      \___ \ / _ \| \ \ / / _ \ '__| |/ _` | '_ \| | |_ | |  | | | |  
 *      ____) | (_) | |\ V /  __/ |  | | (_| | |_) | |__| | |__| |_| |_ 
 *     |_____/ \___/|_| \_/ \___|_|  |_|\__,_|_.__/ \_____|\____/|_____|
 *                                                                      
 *                                                                      
 */
```
