#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
make known as python packages all subdirectories
"""

import os
aDir = os.path.split(__file__)[0]
all_= [o for o in os.listdir(_aDir) if os.path.isdir(os.path.join(_aDir, o))]
del _aDir

#__all__ = ['xxxpy']
