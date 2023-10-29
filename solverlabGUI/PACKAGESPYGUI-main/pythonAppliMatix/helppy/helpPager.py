#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
utilities for pydoc help(...) to file directly

| example:
| > import helppy.helpPager as HP
| > import numpy as np
| > HP.helpToFile(np, "/tmp/${USER}/helpPager.tmp", )
"""

import os
import sys

# to get help() output without standby (esc-:-q as 'more' or 'less')
os.environ["PAGER"] = "cat"

def helpToFile(whatHelp, File="/tmp/${USER}/helpPager.tmp" ):
  out_tmp = sys.stdout
  name = os.path.expandvars(File)
  adir, abasename = os.path.split(name)
  if not os.path.exists(adir):
    os.mkdir(adir)
  sys.stdout = open(name, "w")
  help(whatHelp)
  sys.stdout.flush()
  sys.stdout.close()
  sys.stdout = out_tmp
  return name

