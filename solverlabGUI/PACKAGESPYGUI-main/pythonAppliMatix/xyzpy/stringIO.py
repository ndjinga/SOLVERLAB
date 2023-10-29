#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""\
fixing some StringIO Python 2-3 code but may be cause problem on unicode
"""

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

